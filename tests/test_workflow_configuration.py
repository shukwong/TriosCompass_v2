"""Dry-run the complete workflow using synthetic pedigree and input manifests."""

import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
REQUIRED = ("snakemake", "peds", "networkx", "peppy", "eido", "yaml")
HAS_DEPENDENCIES = all(importlib.util.find_spec(module) is not None for module in REQUIRED)


@unittest.skipUnless(HAS_DEPENDENCIES, "Snakemake and the workflow metadata dependencies are required")
class WorkflowConfigurationTests(unittest.TestCase):
    def setUp(self):
        import yaml
        self.yaml = yaml
        self.temporary = tempfile.TemporaryDirectory(prefix="trios-workflow-")
        self.work = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)
        (self.work / "tmp").mkdir()
        self.ref = self.work / "ref.fa"
        self.ref.write_text(">chr1\n" + "A" * 1000 + "\n")
        Path(str(self.ref) + ".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
        self.ref.with_suffix(".dict").write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n")
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            Path(str(self.ref) + suffix).touch()
        self.interval = self.work / "interval.bed"
        self.interval.write_text("chr1\t0\t1000\n")
        self.exclude = self.work / "legacy_exclude.bed"
        self.exclude.write_text("chr1\t500\t600\n")
        self.ped = self.work / "ped"
        self.ped.mkdir()
        (self.ped / "FAM1.ped").write_text("FAM1 dad 0 0 1 1\nFAM1 mom 0 0 2 1\nFAM1 kid dad mom 1 2\n")
        self.manifest = self.work / "samples.csv"
        with self.manifest.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample_name", "SAMPLE_ID", "BAM"])
            for sample in ("dad", "mom", "kid"):
                bam = self.work / (sample + ".bam")
                bam.touch()
                writer.writerow([sample, sample, str(bam)])
        self.pep = self.work / "pep.yaml"
        self.write_pep()
        self.config = {
            "name": "synthetic_trio", "pepfile": str(self.pep), "output_dir": "output",
            "ref": {"sequence": str(self.ref)}, "ped_dir": str(self.ped),
            "call_dnm": {"interval": str(self.interval)},
            "bam_qc": {"flagstat": {"enable": False}, "collectwgsmetrics": {"enable": False}, "collectmultiplemetrics": {"enable": False}},
            "fastq_input": {"enable": False, "fastqc": {"enable": False}},
            "multiqc": {"enable": False},
        }
        self.env = dict(os.environ, XDG_CACHE_HOME=str(self.work / "cache"), TMPDIR=str(self.work / "tmp"))
        self.env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + self.env["PATH"]

    def write_pep(self):
        self.pep.write_text(self.yaml.safe_dump({"pep_version": "2.0.0", "sample_table": str(self.manifest)}))

    def dryrun(self, targets=(), profile=False, expect_success=True):
        configfile = self.work / "config.yaml"
        configfile.write_text(self.yaml.safe_dump(self.config))
        command = [sys.executable, "-m", "snakemake", "--snakefile", str(REPO / "workflow" / "Snakefile"),
                   "--directory", str(self.work), "--configfile", str(configfile),
                   "--cores", "2", "--dry-run", "--printshellcmds", "--nolock", "--nocolor"]
        if profile:
            self.prepare_profile_runtime()
            command.extend(["--profile", str(REPO / "workflow" / "profiles" / "biowulf")])
        command.extend(targets)
        result = subprocess.run(command, cwd=self.work, env=self.env, text=True, capture_output=True, check=False, timeout=60)
        output = result.stdout + result.stderr
        if expect_success:
            self.assertEqual(result.returncode, 0, output)
            self.assertIn("This was a dry-run", output)
            self.assertFalse((self.work / "output" / "dnm_vcf_summary" / "DNM_summary.txt").exists())
        return result, output

    def test_legacy_config_omits_new_options_and_builds_complete_dag(self):
        _, output = self.dryrun(targets=("all",))
        self.assertIn("rule all:", output)
        self.assertIn("rule gatkhc_pb:", output)
        self.assertIn("rule deepvariant_pb:", output)
        self.assertIn("rule dnm_vcf_summary:", output)
        self.assertIn("bcftools filter -i 'QD", output)
        for name in ("run_rufus", "telomerecat", "bam2bam", "igv_reports_dnsv"):
            self.assertNotIn("rule " + name + ":", output)

    def test_default_target_remains_all_without_multiqc(self):
        _, output = self.dryrun()
        self.assertIn("rule all:", output)
        self.assertIn("output/dnm_vcf_summary/DNM_summary.txt", output)
        self.assertNotIn("rule multiqc:", output)

    def test_gatk_without_hard_filter_renders_the_input_gvcf(self):
        self.config["gatk_hc"] = {"gatk_hard_filter": {"enable": False}}
        _, output = self.dryrun()
        self.assertIn("rule gatk_cgp:", output)
        self.assertIn("output/gatk_genotype_gvcf_pb/FAM1.genotype_gvcf.g.vcf.gz | bcftools view", output)
        self.assertNotIn("bcftools filter -i 'QD", output)

    def test_all_new_modules_build_into_the_full_dag(self):
        self.config["rufus"] = {"enable": True, "threads": 2}
        self.config["igv_reports"] = {"enable": True, "source": "rufus"}
        self.config["bam_qc"].update({"telomerecat": {"enable": True, "threads": 2}, "bam2bam": {"enable": True, "threads": 2}})
        _, output = self.dryrun()
        for name in ("all", "run_rufus", "rufus_summary", "telomerecat", "bam2bam", "igv_reports_dnsv", "dnm_vcf_summary"):
            self.assertIn("rule " + name + ":", output)
        self.assertIn("output/rufus/FAM1.rufus.vcf", output)
        self.assertNotIn("--in-bam output/bam2bam", output)

    def test_igv_sources_connect_to_their_callers(self):
        for source, caller, vcf in (
            ("joint", "joint_dnSV", "output/joint_dnSV/FAM1.dnSV.vcf"),
            ("manta", "dnSV_manta", "output/dnSV_manta/FAM1.dnSV.manta.vcf"),
            ("rufus", "run_rufus", "output/rufus/FAM1.rufus.vcf"),
        ):
            with self.subTest(source=source):
                self.config["igv_reports"] = {"enable": True, "source": source}
                self.config["rufus"] = {"enable": source == "rufus"}
                self.config["dnSV"] = {"enable": source != "rufus", "exclude_bed": str(self.exclude),
                    "GraphTyper2+smoove": {"enable": source == "joint"}, "Manta": {"enable": source == "manta"}}
                _, output = self.dryrun()
                self.assertIn("rule " + caller + ":", output)
                self.assertIn("rule igv_reports_dnsv:", output)
                self.assertIn("--vcf " + vcf, output)

    def set_fastq_manifest(self, invalid_lane=False):
        self.config["fastq_input"]["enable"] = True
        self.config["pepschema"] = str(REPO / "workflow" / "schemas" / "rerf_fastq_schema.yaml")
        with self.manifest.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample_name", "SAMPLE_ID", "FLOWCELL", "LANE", "INDEX", "R1", "R2"])
            for sample in ("dad", "mom", "kid"):
                lanes = range(1, 10 if invalid_lane else 9) if sample == "dad" else (1,)
                for lane in lanes:
                    reads = []
                    for read in (1, 2):
                        path = self.work / f"{sample}_FC1_L{lane}.R{read}.fastq.gz"
                        path.touch()
                        reads.append(str(path))
                    writer.writerow([f"{sample}_FC1_L{lane}", sample, "FC1", lane, "IDX001", *reads])
        self.write_pep()

    def test_fastq_mode_accepts_all_eight_lanes(self):
        self.set_fastq_manifest()
        _, output = self.dryrun()
        self.assertIn("rule fq2bam:", output)
        self.assertIn("dad_FC1_L8", output)
        self.assertIn("rule dnm_vcf_summary:", output)
        self.assertIn("--num-gpus 1 --low-memory", output)

    def test_fastq_schema_rejects_lane_nine(self):
        self.set_fastq_manifest(invalid_lane=True)
        result, output = self.dryrun(expect_success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("EidoValidationError", output)
        self.assertIn("Validation unsuccessful", output)

    def test_null_top_level_exclusion_keeps_the_legacy_path(self):
        self.config["dnSV"] = {"enable": True, "exclude_bed": None,
            "GraphTyper2+smoove": {"enable": True, "exclude_bed": str(self.exclude)}, "Manta": {"enable": False}}
        _, output = self.dryrun()
        self.assertIn("rule smoove:", output)
        self.assertIn("--exclude " + str(self.exclude), output)

    def prepare_profile_runtime(self):
        # Snakemake probes deployment tools even in dry-run mode. These metadata-
        # only stubs reject any pull, environment creation, or scheduler action.
        bindir = self.work / "probe_tools"
        bindir.mkdir()
        source = r'''import json
import os
from pathlib import Path
import sys
name = Path(sys.argv[0]).name
args = sys.argv[1:]
if name == "singularity" and args == ["--version"]:
    print("singularity version 3.11.4")
elif name in ("conda", "mamba") and args == ["--version"]:
    print(name + " 24.1.0")
elif name == "conda" and args == ["info", "--json"]:
    print(json.dumps({"conda_prefix": os.environ["MOCK_CONDA_PREFIX"], "platform": "linux-64"}))
elif name == "conda" and args == ["config", "--get", "channel_priority", "--json"]:
    print(json.dumps({"get": {"channel_priority": "strict"}}))
else:
    sys.exit("Unexpected execution in workflow dry-run: " + name + " " + " ".join(args))
'''
        for name in ("singularity", "conda", "mamba", "sbatch", "scancel"):
            executable = bindir / name
            executable.write_text("#!" + sys.executable + "\n" + source)
            executable.chmod(0o700)
        self.env["PATH"] = str(bindir) + os.pathsep + self.env["PATH"]
        self.env["MOCK_CONDA_PREFIX"] = str(self.work / "mock_conda")

    def test_biowulf_profile_resolves_for_a_dry_run(self):
        from snakemake import get_profile_file
        profile = REPO / "workflow" / "profiles" / "biowulf"
        settings = self.yaml.safe_load((profile / "config.yaml").read_text())
        for key, filename in (("cluster", "slurm-submit.py"), ("cluster-status", "slurm-status.py"), ("jobscript", "slurm-jobscript.sh")):
            resolved = Path(get_profile_file(str(profile), settings[key])).resolve()
            self.assertEqual(resolved, (profile.parent / "slurm" / filename).resolve())
            self.assertTrue(resolved.is_file())
            self.assertTrue(os.access(resolved, os.X_OK))
        _, output = self.dryrun(profile=True)
        self.assertIn("rule all:", output)
        self.assertIn("partition=gpu", output)
        self.assertIn("mem_mb=200000", output)


if __name__ == "__main__":
    unittest.main()
