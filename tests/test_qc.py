"""Execute QC rules against synthetic BAMs and lightweight external-tool stubs."""

import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from picard_launcher_fixture import PICARD_LAUNCHER


RULE = Path(__file__).resolve().parents[1] / "workflow" / "rules" / "bam_qc.smk"
SAMTOOLS = shutil.which("samtools")
SNAKEMAKE_PYTHON = os.environ.get("SNAKEMAKE_PYTHON", sys.executable)
HAS_SNAKEMAKE = importlib.util.find_spec("snakemake") is not None


@unittest.skipUnless(SAMTOOLS and HAS_SNAKEMAKE, "samtools and the Snakemake Python environment are required")
class QcRuleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="trios-qc-")
        self.work = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)
        self.bindir = self.work / "bin"
        self.bindir.mkdir()
        (self.work / "tmp").mkdir()
        self.bam = self.work / "input reads.bam"
        self.ref = self.work / "reference genome.fa"
        self.ref.write_text(">chr1\n" + "A" * 1000 + "\n")
        self.trace = self.work / "tool_trace.jsonl"
        self.env = dict(
            os.environ,
            PATH=str(self.bindir) + os.pathsep + os.environ["PATH"],
            TOOL_TRACE=str(self.trace),
            SAMTOOLS_REAL=SAMTOOLS,
            XDG_CACHE_HOME=str(self.work / "cache"),
            TMPDIR=str(self.work / "tmp"),
        )
        self.install("samtools", '''import os
import json
import sys
with open(os.environ["TOOL_TRACE"], "a") as trace:
    trace.write(json.dumps({"tool": "samtools", "args": sys.argv[1:]}) + "\\n")
os.execv(os.environ["SAMTOOLS_REAL"], [os.environ["SAMTOOLS_REAL"], *sys.argv[1:]])
''')

    def install(self, name, source):
        path = self.bindir / name
        path.write_text("#!" + sys.executable + "\n" + source)
        path.chmod(0o700)
        return path

    def make_bam(self, records=()):
        header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        lines = []
        for number, (name, flag, mapq, cigar) in enumerate(records, start=1):
            lines.append(f"{name}\t{flag}\tchr1\t{number * 20}\t{mapq}\t{cigar}\t=\t900\t100\tAAAAAAAAAA\tIIIIIIIIII\n")
        subprocess.run([SAMTOOLS, "view", "-b", "-o", str(self.bam), "-"], input=header + "".join(lines), text=True, check=True, capture_output=True)
        subprocess.run([SAMTOOLS, "index", str(self.bam)], check=True, capture_output=True)

    def run_rules(self, qc, threads=None, check=True):
        config = {"bam_qc": qc, "threads": threads or {}, "parabricks": {"singularity_cmd": ""}}
        snakefile = self.work / "Snakefile"
        snakefile.write_text(
            "import os\n"
            + "config = " + repr(config) + "\n"
            + "output_dir = 'QC output with spaces'\n"
            + "genome = " + repr(str(self.ref)) + "\n"
            + "subjs = ['sample']\noptional_output = []\nqc_output = []\n"
            + "def get_bam_by_subj(wildcards):\n    return " + repr(str(self.bam)) + "\n"
            + "include: " + repr(str(RULE.with_name("runtime.smk"))) + "\n"
            + "include: " + repr(str(RULE)) + "\n"
            + "rule all:\n    default_target: True\n    input: optional_output\n"
        )
        result = subprocess.run(
            [SNAKEMAKE_PYTHON, "-m", "snakemake", "--snakefile", str(snakefile),
             "--directory", str(self.work), "--cores", "2", "--jobs", "1", "--nolock", "--nocolor", "all"],
            cwd=self.work, env=self.env, text=True, capture_output=True, check=False,
        )
        if check:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def traces(self, tool):
        return [row for row in map(json.loads, self.trace.read_text().splitlines()) if row["tool"] == tool]

    def output(self, relative):
        return self.work / "QC output with spaces" / relative

    def read_names(self, path):
        data = subprocess.run([SAMTOOLS, "view", str(path)], text=True, check=True, capture_output=True).stdout
        return [line.split("\t")[0] for line in data.splitlines()]

    def test_bam_filter_matches_reference_and_preserves_source(self):
        records = [
            ("pass_m", 99, 60, "10M"), ("duplicate", 99 + 0x400, 60, "10M"),
            ("secondary", 99 + 0x100, 60, "10M"), ("low_mapq", 99, 59, "10M"),
            ("not_proper_pair", 97, 60, "10M"), ("soft_clip", 99, 60, "1S9M"),
            ("hard_clip", 99, 60, "1H10M"), ("padding", 99, 60, "5M1P5M"),
            ("ref_skip", 99, 60, "5M1N5M"), ("eq_x_only", 99, 60, "5=5X"),
            ("insertion", 99, 60, "5M1I4M"), ("deletion", 99, 60, "5M1D5M"),
            ("supplementary", 99 + 0x800, 60, "10M"),
        ]
        self.make_bam(records)
        before = self.bam.read_bytes()
        self.run_rules({"bam2bam": {"enable": True, "threads": 2, "mem_mb": 1024}, "flagstat": {"enable": True}})
        redup = self.output("bam2bam/sample.triocompass.addRG.redup.bam")
        hmq = self.output("bam2bam/sample.triocompass.addRG.redup.HMQ.sorted.bam")
        self.assertEqual(self.bam.read_bytes(), before)
        self.assertEqual(set(self.read_names(redup)), {row[0] for row in records} - {"duplicate"})
        self.assertEqual(self.read_names(hmq), ["pass_m", "insertion", "deletion", "supplementary"])
        original_view = subprocess.run([SAMTOOLS, "view", "-f", "0x2", "-h", "-q", "60", "-F", "0x500", str(redup)], text=True, capture_output=True, check=True).stdout
        original_awk = subprocess.run(["awk", '$1~/^@/||($6!~/[SHPN]/ && $6~/M/)'], input=original_view, text=True, capture_output=True, check=True).stdout
        expected = [line.split("\t")[0] for line in original_awk.splitlines() if not line.startswith("@")]
        self.assertEqual(self.read_names(hmq), expected)
        for bam in (redup, hmq):
            self.assertTrue(Path(str(bam) + ".bai").exists())
            subprocess.run([SAMTOOLS, "quickcheck", str(bam)], check=True, capture_output=True)
            subprocess.run([SAMTOOLS, "idxstats", str(bam)], check=True, capture_output=True)
        calls = self.traces("samtools")
        for row in calls:
            args = row["args"]
            if "-@" in args:
                self.assertEqual(args[args.index("-@") + 1], "1")
        sort_args = next(row["args"] for row in calls if row["args"][0] == "sort")
        self.assertEqual(sort_args[sort_args.index("-m") + 1], "358M")
        self.assertIn("idxstats", {row["args"][0] for row in calls})
        self.assertTrue(self.output("flagstat/sample.idxstat").is_file())
        self.assertEqual(list(self.output("bam2bam/work/sample").glob("attempt.*")), [])

    def test_empty_bam_produces_valid_empty_indexed_outputs(self):
        self.make_bam()
        self.run_rules({"bam2bam": {"enable": True, "threads": 1, "mem_mb": 1024}})
        for suffix in ("redup.bam", "redup.HMQ.sorted.bam"):
            bam = self.output("bam2bam/sample.triocompass.addRG." + suffix)
            self.assertEqual(self.read_names(bam), [])
            self.assertTrue(Path(str(bam) + ".bai").exists())
            subprocess.run([SAMTOOLS, "quickcheck", str(bam)], check=True, capture_output=True)
            subprocess.run([SAMTOOLS, "idxstats", str(bam)], check=True, capture_output=True)
        for row in self.traces("samtools"):
            args = row["args"]
            if "-@" in args:
                self.assertEqual(args[args.index("-@") + 1], "0")

    def install_telomerecat(self):
        return self.install("telomerecat custom", '''import json
import os
from pathlib import Path
import sys
args = sys.argv[1:]
with open(os.environ["TOOL_TRACE"], "a") as trace:
    trace.write(json.dumps({"tool": "telomerecat", "args": args, "cwd": os.getcwd(), "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")}) + "\\n")
output = Path(args[args.index("--output") + 1])
output.write_text("" if os.environ.get("MOCK_MODE") == "empty" else "Sample,Length\\nsample,1234\\n")
print("telomerecat diagnostic")
if os.environ.get("MOCK_MODE") == "fail":
    sys.exit(27)
''')

    def test_telomerecat_uses_allocated_threads_and_staged_absolute_paths(self):
        self.make_bam()
        executable = self.install_telomerecat()
        self.run_rules({"telomerecat": {"enable": True, "threads": 4, "executable": str(executable)}})
        call, = self.traces("telomerecat")
        args = call["args"]
        self.assertEqual(args[:3], ["bam2length", "-p", "2"])
        self.assertEqual(args[-1], str(self.bam))
        self.assertIn("-v2", args)
        self.assertEqual(call["OMP_NUM_THREADS"], "1")
        self.assertEqual(Path(args[args.index("--output") + 1]).parent, Path(call["cwd"]))
        self.assertIn("attempt.", Path(call["cwd"]).name)
        self.assertEqual(self.output("telomerecat/sample.csv").read_text(), "Sample,Length\nsample,1234\n")
        self.assertIn("telomerecat diagnostic", self.output("logs/telomerecat/sample.log").read_text())

    def test_telomerecat_failure_does_not_publish_final_csv(self):
        self.make_bam()
        executable = self.install_telomerecat()
        self.env["MOCK_MODE"] = "fail"
        result = self.run_rules({"telomerecat": {"enable": True, "executable": str(executable)}}, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output("telomerecat/sample.csv").exists())
        self.assertIn("telomerecat diagnostic", self.output("logs/telomerecat/sample.log").read_text())

    def test_telomerecat_empty_csv_is_rejected(self):
        self.make_bam()
        executable = self.install_telomerecat()
        self.env["MOCK_MODE"] = "empty"
        result = self.run_rules({"telomerecat": {"enable": True, "executable": str(executable)}}, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output("telomerecat/sample.csv").exists())

    def install_picard(self):
        self.install("picard", PICARD_LAUNCHER)

    def test_picard_receives_intact_arguments_and_bounded_heap(self):
        self.make_bam()
        self.install_picard()
        self.run_rules({"collectwgsmetrics": {"enable": True, "mem_mb": 1024}})
        call, = self.traces("picard")
        self.assertEqual(call["program_args"], ["CollectWgsMetrics", "I=" + str(self.bam), "O=QC output with spaces/collectwgsmetrics/sample.collect_wgs_metrics.txt", "R=" + str(self.ref)])
        self.assertEqual(call["max_heap_bytes"], 819 * 1024 ** 2)
        self.assertEqual(call["initial_heap_bytes"], 1024 ** 2)
        self.assertEqual(call["java_options"], ["-Xms1m", "-Xmx819m", "-XX:+UseSerialGC"])
        self.assertEqual(call["injected"], dict.fromkeys(("JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS")))
        self.assertEqual(self.output("collectwgsmetrics/sample.collect_wgs_metrics.txt").read_text(), "metrics\n")

    def test_picard_failure_removes_incomplete_metrics(self):
        self.make_bam()
        self.install_picard()
        self.env["MOCK_MODE"] = "fail"
        result = self.run_rules({"collectwgsmetrics": {"enable": True, "mem_mb": 1024}}, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output("collectwgsmetrics/sample.collect_wgs_metrics.txt").exists())


if __name__ == "__main__":
    unittest.main()
