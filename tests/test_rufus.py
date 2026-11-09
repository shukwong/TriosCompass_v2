"""RUFUS integration contract tests: synthetic BAM headers and a fake caller."""
import gzip
import importlib.util
import json
from pathlib import Path
import struct
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/run_rufus.py"
SPEC = importlib.util.spec_from_file_location("run_rufus", SCRIPT)
rufus = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rufus)
HEADER = ("##fileformat=VCFv4.2\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t./child.bam\tfather.bam\tmother.bam\n")
RECORDS = (
    "chr1\t10\tsnv\tA\tG\t30\tPASS\t.\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t20\tins\tA\tAT\t30\tPASS\tSVTYPE=INS\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t30\tsv\tN\t<DEL>\t30\tPASS\tSVTYPE=DEL;END=130\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t40\tbnd\tN\tN]chr2:20]\t.\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\t0/0\n"
)


class RufusTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ped = self.root / "trio.ped"
        self.ped.write_text("F mother 0 0 2 1\nF child father mother 1 2\nF father 0 0 1 1\n")
        self.samples = ["mother", "child", "father"]
        # Deliberately use identical basenames and spaces in the original paths.
        self.bams = [self.bam(sample) for sample in self.samples]
        self.reference = self.root / "reference input.fa"
        self.reference.write_text(">chr1\nACGT\n")
        for suffix in (".fai",) + rufus.BWA_SUFFIXES:
            Path(str(self.reference) + suffix).write_text("index\n")
        self.output = self.root / "output" / "F.rufus.vcf"
        self.fake = self.root / "fake caller.py"
        self.fake.write_text(
            "#!/usr/bin/env python3\n"
            "import gzip, json, pathlib, sys\n"
            "pathlib.Path('argv.json').write_text(json.dumps(sys.argv[1:]))\n"
            "assert all(pathlib.Path(p).is_file() for p in "
            "['child.bam', 'father.bam', 'mother.bam', 'child.bam.bai', 'father.bam.bai', "
            "'mother.bam.bai', 'reference.fa', 'reference.fa.fai', 'reference.fa.bwt'])\n"
            "with gzip.open(" + repr(rufus.FINAL_NAME + ".gz") + ", 'wt') as f:\n"
            "    f.write(" + repr(HEADER + RECORDS) + ")\n"
        )
        self.fake.chmod(0o700)

    def bam(self, sample, filename="identical.bam"):
        path = self.root / (sample + " input") / filename
        path.parent.mkdir(exist_ok=True)
        header = ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:4\n"
                  "@RG\tID:reads\tSM:" + sample + "\n").encode()
        with gzip.open(path, "wb") as handle:
            handle.write(b"BAM\x01" + struct.pack("<i", len(header)) + header + struct.pack("<i", 0))
        Path(str(path) + ".bai").write_bytes(b"BAI\x01" + struct.pack("<i", 0))
        return str(path)

    def args(self, **settings):
        config = {"executable": str(self.fake)}
        config.update(settings)
        return SimpleNamespace(ped=str(self.ped), family="F", samples=self.samples,
                               bams=self.bams, reference=str(self.reference), threads=3,
                               settings_json=json.dumps(config), work_dir=str(self.root / "work"),
                               output=str(self.output))

    def test_fake_caller_orders_trio_preserves_records_and_keeps_attempts(self):
        self.assertEqual(rufus.run(self.args()), 0)
        self.assertEqual(rufus.run(self.args()), 0)
        attempts = sorted((self.root / "work/F").iterdir())
        self.assertEqual(len(attempts), 2)
        for attempt in attempts:
            provenance = json.loads((attempt / "run.json").read_text())
            self.assertEqual(provenance["status"], "complete")
            self.assertEqual([m["sample"] for m in provenance["trio"]], ["child", "father", "mother"])
            argv = json.loads((attempt / "argv.json").read_text())
            self.assertEqual(argv, ["-s", "child.bam", "-c", "father.bam", "-c", "mother.bam",
                                    "-r", "reference.fa", "-t", "3", "-k", "25"])
            self.assertTrue((attempt / (rufus.FINAL_NAME + ".gz")).is_file())
            self.assertEqual((attempt / "mother.bam").resolve(), Path(self.bams[0]).resolve())
        expected_header = HEADER.replace("./child.bam\tfather.bam\tmother.bam", "child\tfather\tmother")
        self.assertEqual(self.output.read_text(), expected_header + RECORDS)

    def test_summary_counts_records_and_classes_without_claiming_all_are_svs(self):
        source = self.root / "mixed.vcf"
        source.write_text(HEADER + RECORDS + "chr1\t50\tmixed\tA\tG,AT\t.\tPASS\t.\tGT\t1/2\t0/0\t0/0\n")
        summary = self.root / "summary.tsv"
        self.assertEqual(rufus.summary(SimpleNamespace(families=["F"], vcfs=[str(source)], output=str(summary))), 0)
        columns, row = [line.split("\t") for line in summary.read_text().splitlines()]
        record = dict(zip(columns, row))
        self.assertEqual(record["TotalRecords"], "5")
        for key in ("SNV", "INDEL", "SYMBOLIC", "BND", "MIXED"):
            self.assertEqual(record[key], "1")
        self.assertEqual(json.loads(record["SVTYPE_counts"]), {"INS": 1, "DEL": 1, "BND": 1})
        self.assertNotIn("dnSV", summary.read_text())

    def test_header_only_final_vcf_is_valid(self):
        source = self.root / "empty.vcf"
        source.write_text(HEADER)
        self.assertEqual(rufus.validate_vcf(source)["total_records"], 0)

    def test_missing_parent_and_duplicate_bam_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "both distinct parents"):
            rufus.trio_from_ped(self.ped, "F", self.samples[:2], self.bams[:2])
        with self.assertRaisesRegex(ValueError, "distinct BAM"):
            rufus.trio_from_ped(self.ped, "F", self.samples, [self.bams[0]] * 3)

    def test_mismatched_bam_sample_fails_before_caller(self):
        args = self.args()
        args.bams = [self.bams[1], self.bams[0], self.bams[2]]
        with mock.patch.object(rufus.subprocess, "run") as caller:
            with self.assertRaisesRegex(ValueError, "BAM sample mismatch"):
                rufus.run(args)
        caller.assert_not_called()
        manifest = next((self.root / "work/F").glob("*/run.json"))
        self.assertEqual(json.loads(manifest.read_text())["status"], "failed")
        self.assertFalse(self.output.exists())

    def test_stages_index_adjacent_to_bam_symlink(self):
        alias = self.root / "child-alias.bam"
        alias.symlink_to(self.bams[1])
        original_index = Path(self.bams[1] + ".bai")
        alias_index = Path(str(alias) + ".bai")
        original_index.rename(alias_index)
        args = self.args()
        args.bams = [self.bams[0], str(alias), self.bams[2]]
        self.assertEqual(rufus.run(args), 0)

    def test_index_symlink_keeps_bai_suffix_when_target_has_different_name(self):
        index = Path(self.bams[1] + ".bai")
        target = index.parent / "stored-index"
        index.rename(target)
        index.symlink_to(target)
        alias = self.root / "child-without-adjacent-index.bam"
        alias.symlink_to(self.bams[1])
        args = self.args()
        args.bams = [self.bams[0], str(alias), self.bams[2]]
        self.assertEqual(rufus.run(args), 0)
        attempt = next((self.root / "work/F").iterdir())
        self.assertEqual((attempt / "child.bam.bai").resolve(), target)

    def test_malformed_or_wrong_sample_vcf_never_replaces_existing_output(self):
        self.output.parent.mkdir()
        self.output.write_text("previous result\n")
        source = self.root / "invalid.vcf"
        trio = [{"staged_bam": role + ".bam", "sample": role} for role in ("child", "father", "mother")]
        for content in ("not a VCF\n", HEADER + "chr1\tbroken\n",
                        HEADER.replace("./child.bam", "unrelated.bam") + RECORDS,
                        HEADER + RECORDS.replace("chr1\t10", "chr1\t0", 1)):
            with self.subTest(content=content[:40]):
                source.write_text(content)
                with self.assertRaises(ValueError):
                    rufus.publish_vcf(source, self.output, trio)
                self.assertEqual(self.output.read_text(), "previous result\n")
        self.assertEqual(list(self.output.parent.glob(".F.rufus.vcf.*")), [])

    def test_caller_failure_keeps_logs_and_does_not_publish(self):
        self.fake.write_text("#!/usr/bin/env python3\nimport sys\nprint('caller failed')\nsys.exit(7)\n")
        with self.assertRaises(rufus.subprocess.CalledProcessError):
            rufus.run(self.args())
        attempt = next((self.root / "work/F").iterdir())
        self.assertIn("caller failed", (attempt / "caller.log").read_text())
        self.assertEqual(json.loads((attempt / "run.json").read_text())["status"], "failed")
        self.assertFalse(self.output.exists())

    def test_missing_final_vcf_is_failure(self):
        self.fake.write_text("#!/usr/bin/env python3\nprint('no final output')\n")
        with self.assertRaisesRegex(ValueError, "did not create"):
            rufus.run(self.args())
        self.assertFalse(self.output.exists())

    def test_container_command_binds_inputs_and_preserves_argument_boundaries(self):
        runtime = self.root / "fake-apptainer"
        runtime.write_text("#!/bin/sh\nexit 0\n")
        runtime.chmod(0o700)
        container = self.root / "rufus image.sif"
        container.touch()
        work = self.root / "work dir"
        work.mkdir()
        command = rufus.caller_command(
            {"container": str(container), "runtime": str(runtime), "runtime_args": ["--cleanenv"],
             "executable": "/opt/rufus/runRufus.sh", "extra_args": ["-m", "5"]},
            3, work, [self.reference])
        self.assertEqual(command[:3], [str(runtime), "exec", "--cleanenv"])
        self.assertIn(str(work) + ":" + str(work) + ":rw", command)
        self.assertIn(str(self.root) + ":" + str(self.root) + ":ro", command)
        self.assertEqual(command[command.index(str(container)) + 1], "/opt/rufus/runRufus.sh")
        self.assertEqual(command[-2:], ["-m", "5"])

    def test_core_arguments_and_minimum_threads_cannot_be_silently_overridden(self):
        for extra in (["--subject=other.bam"], ["-CLEAN"], ["--threads", "100"]):
            with self.assertRaisesRegex(ValueError, "cannot override"):
                rufus.caller_command({"extra_args": extra}, 3, self.root, [])
        with self.assertRaisesRegex(ValueError, ">= 3"):
            rufus.caller_command({}, 2, self.root, [])
        with self.assertRaisesRegex(ValueError, "list of strings"):
            rufus.caller_command({"extra_args": "-m 5"}, 3, self.root, [])


if __name__ == "__main__":
    unittest.main()
