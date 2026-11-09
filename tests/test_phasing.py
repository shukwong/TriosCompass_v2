"""Exercise the actual phasing shell rules, including large and empty scatters."""

import ast
import gzip
import os
from pathlib import Path
import re
import resource
import shlex
import shutil
import subprocess
import tempfile
import unittest


RULE = Path(__file__).resolve().parents[1] / "workflow/rules/phasing.smk"
HEADER = ("VariantID\tParentalOrigin\tPhase\tHaplotypeBlockSize\t InformativeSiteCnt\t"
          "FM_count\tMF_count\tFM_status\tParentOriginChange\n")


class ShellPaths:
    """Format Snakemake's input/output path lists for the isolated shell tests."""

    def __init__(self, paths):
        self.paths = [str(path) for path in paths]

    def __format__(self, specification):
        if specification not in ("", "q"):
            raise ValueError(specification)
        quote = shlex.quote if specification == "q" else str
        return " ".join(quote(path) for path in self.paths)


def rule_shell(name):
    source = RULE.read_text()
    block = re.search(r"^(?:rule|checkpoint) " + re.escape(name)
                      + r":\n(.*?)(?=^(?:rule|checkpoint|def) |\Z)",
                      source, re.MULTILINE | re.DOTALL).group(1)
    literal = re.search(r'    shell:\s*(""".*?""")', block, re.DOTALL).group(1)
    return ast.literal_eval(literal)


@unittest.skipUnless(shutil.which("bash"), "bash is required")
class PhasingRuleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="trios-phasing-")
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.vcf = self.work / "input calls.vcf.gz"
        self.variants = self.work / "scattered variants"
        self.output = self.work / "parental origin.tab"

    def write_vcf(self, records=()):
        with gzip.open(self.vcf, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, pos, ref, alt in records:
                handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

    def execute(self, name, inputs, output, stdin="", descriptor_limit=None):
        command = ["bash", "--noprofile", "--norc", "-euo", "pipefail", "-c",
                   rule_shell(name).format(input=ShellPaths(inputs), output=ShellPaths([output]))]
        # Optional container prefix keeps these same checks usable with the
        # workflow's exact utility image, without requiring it for unit tests.
        prefix = shlex.split(os.environ.get("TRIOS_PHASING_TEST_PREFIX", ""))
        limit = None
        if descriptor_limit is not None:
            def limit():
                _, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
                resource.setrlimit(resource.RLIMIT_NOFILE, (min(descriptor_limit, hard), hard))
        return subprocess.run(prefix + command, cwd=self.work, input=stdin,
                              text=True, capture_output=True, timeout=30,
                              preexec_fn=limit)

    def check_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_scatter_closes_files_and_retains_duplicate_locus_records(self):
        records = [("chr1", pos, "A", "G") for pos in range(1, 2049)]
        records += [("chr1", 1, "A", "T"), ("chr2", 10, "AT", "A")]
        self.write_vcf(records)
        before = self.vcf.read_bytes()
        self.check_success(self.execute("scatter_dnms", [self.vcf], self.variants, descriptor_limit=64))
        self.assertEqual(len(list(self.variants.glob("*.dnm"))), 2049)
        self.assertEqual((self.variants / "chr1_1.dnm").read_text(), "chr1:1:A:G\nchr1:1:A:T\n")
        self.assertEqual((self.variants / "chr1_2048.dnm").read_text(), "chr1:2048:A:G\n")
        self.assertEqual((self.variants / "chr2_10.dnm").read_text(), "chr2:10:AT:A\n")
        self.assertEqual(sum(len(path.read_text().splitlines()) for path in self.variants.glob("*.dnm")), len(records))
        # Re-execution truncates each locus once rather than appending another copy.
        self.check_success(self.execute("scatter_dnms", [self.vcf], self.variants, descriptor_limit=64))
        self.assertEqual((self.variants / "chr1_1.dnm").read_text(), "chr1:1:A:G\nchr1:1:A:T\n")
        self.assertEqual(self.vcf.read_bytes(), before)

    def test_header_only_vcf_scatter_and_aggregate_are_successful(self):
        self.write_vcf()
        self.check_success(self.execute("scatter_dnms", [self.vcf], self.variants))
        self.assertTrue(self.variants.is_dir())
        self.assertEqual(list(self.variants.iterdir()), [])
        self.check_success(self.execute("aggregate_phase", [], self.output, stdin="must not read stdin\n"))
        self.assertEqual(self.output.read_text(), HEADER)

    def test_aggregate_preserves_header_and_natural_variant_order(self):
        first = self.work / "first parental origin.txt"
        second = self.work / "second parental origin.txt"
        first.write_text("chr10:9:A:G\tF\nchr2:100:A:G\tM\n")
        second.write_text("chr2:9:A:G\tF\nchr1:20:A:G\tM\n")
        self.check_success(self.execute("aggregate_phase", [first, second], self.output))
        self.assertEqual(self.output.read_text(), HEADER + "chr1:20:A:G\tM\nchr2:9:A:G\tF\nchr2:100:A:G\tM\nchr10:9:A:G\tF\n")

    def test_invalid_compressed_input_still_fails(self):
        self.vcf.write_bytes(b"not a gzip stream\n")
        result = self.execute("scatter_dnms", [self.vcf], self.variants)
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
