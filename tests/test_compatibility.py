"""Behavior checks for workflow compatibility fixes without running bioinformatics tools."""

import ast
import os
from pathlib import Path
import re
import shlex
import string
import subprocess
import tempfile
from types import SimpleNamespace
import unittest


RULES = Path(__file__).resolve().parents[1] / "workflow" / "rules"


class ShellFormatter(string.Formatter):
    def format_field(self, value, format_spec):
        if format_spec == "q":
            return shlex.quote(str(value))
        return super().format_field(value, format_spec)


def rule_shell(source, name):
    block = re.split(r"(?m)^rule " + re.escape(name) + r":\s*$", source, maxsplit=1)[1]
    literal = re.search(r'(?ms)^    shell:\s*(""".*?""")', block).group(1)
    return ast.literal_eval(literal)


class GatkMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        scope = {}
        source = (RULES / "gatk_hc.smk").read_text().split("\nrule gatkhc_pb:", 1)[0]
        exec(compile(source, "gatk_hc.smk", "exec"), scope)
        cls.heap = staticmethod(scope["_gatk_java_heap_mb"])

    def test_heap_fits_small_and_large_allocations(self):
        for mem_mb in (512, 1024, 4096, 16000, 64000):
            with self.subTest(mem_mb=mem_mb):
                heap_mb = self.heap(None, SimpleNamespace(mem_mb=mem_mb))
                self.assertGreater(heap_mb, 0)
                self.assertLess(heap_mb, mem_mb)
                self.assertGreaterEqual(mem_mb - heap_mb, 256)

    def test_missing_memory_has_a_bounded_default(self):
        heap_mb = self.heap(None, SimpleNamespace())
        self.assertGreater(heap_mb, 0)
        self.assertLess(heap_mb, 16000)

    def test_insufficient_memory_fails_clearly(self):
        for mem_mb in (-1, 0, 256, 511):
            with self.subTest(mem_mb=mem_mb):
                with self.assertRaisesRegex(ValueError, "at least 512 MB"):
                    self.heap(None, SimpleNamespace(mem_mb=mem_mb))


class HipstrStagingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = (RULES / "hipstr.smk").read_text()
        cls.shells = {name: rule_shell(source, name) for name in ("hipstr", "hipstr_recall")}

    def run_hipstr(self, rule, failure):
        with tempfile.TemporaryDirectory(prefix="hipstr-test-") as tmp:
            work = Path(tmp)
            bindir = work / "bin"
            bindir.mkdir()
            hipstr = bindir / "HipSTR"
            hipstr.write_text('''#!/usr/bin/env python3
import os
from pathlib import Path
import sys
args = sys.argv[1:]
def option(flag):
    return args[args.index(flag) + 1]
Path(option("--log")).write_text("HipSTR diagnostic log\\n")
Path(option("--str-vcf")).write_text("VCF data\\n")
Path(option("--viz-out")).write_text("alignment data\\n")
if os.environ.get("FAIL_TOOL") == "hipstr":
    sys.exit(23)
''')
            tabix = bindir / "tabix"
            tabix.write_text('''#!/usr/bin/env python3
import os
from pathlib import Path
import sys
if os.environ.get("FAIL_TOOL") == "tabix":
    sys.exit(29)
Path(sys.argv[-1] + ".tbi").write_text("index data\\n")
''')
            for executable in (hipstr, tabix):
                executable.chmod(0o700)
            output_dir = work / "configured output with spaces"
            output = SimpleNamespace(
                vcf=str(output_dir / "calls.vcf.gz"),
                vcf_index=str(output_dir / "calls.vcf.gz.tbi"),
                viz=str(output_dir / "calls.aln.viz.gz"),
                viz_index=str(output_dir / "calls.aln.viz.gz.tbi"),
            )
            log = str(work / "logs with spaces" / "hipstr.log")
            command = ShellFormatter().format(
                self.shells[rule], output=output, log=log,
                input=SimpleNamespace(ref="reference with spaces.fa", reg="regions.bed"),
                params=SimpleNamespace(bams="father.bam,mother.bam,child.bam"),
            )
            env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ["PATH"], FAIL_TOOL=failure)
            result = subprocess.run(
                ["bash", "-euo", "pipefail", "-c", command], cwd=work, env=env,
                text=True, capture_output=True, check=False,
            )
            expected = {"": 0, "hipstr": 23, "tabix": 29}[failure]
            self.assertEqual(result.returncode, expected, result.stderr)
            self.assertEqual(Path(log).read_text(), "HipSTR diagnostic log\n")
            self.assertEqual(list(output_dir.glob(".tmp.*")), [])
            for target in vars(output).values():
                self.assertEqual(Path(target).exists(), not failure, target)
            if not failure:
                for target in vars(output).values():
                    self.assertGreater(Path(target).stat().st_size, 0)

    def test_success_publishes_all_outputs_at_configured_path(self):
        for rule in self.shells:
            with self.subTest(rule=rule):
                self.run_hipstr(rule, "")

    def test_hipstr_failure_preserves_status_and_log_without_publishing(self):
        for rule in self.shells:
            with self.subTest(rule=rule):
                self.run_hipstr(rule, "hipstr")

    def test_index_failure_does_not_publish_partial_outputs(self):
        for rule in self.shells:
            with self.subTest(rule=rule):
                self.run_hipstr(rule, "tabix")


if __name__ == "__main__":
    unittest.main()
