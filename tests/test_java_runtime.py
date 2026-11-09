"""Execute the Java rule shells against runtime and launcher boundary fixtures."""

import ast
import json
import os
from pathlib import Path
import re
import shlex
import string
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

from picard_launcher_fixture import PICARD_LAUNCHER


RULES = Path(__file__).resolve().parents[1] / "workflow" / "rules"
JAVA_OPTIONS = ("JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS")


class ShellFormatter(string.Formatter):
    def format_field(self, value, format_spec):
        if format_spec == "q":
            return shlex.quote(str(value))
        return super().format_field(value, format_spec)


def rule_block(filename, name):
    source = (RULES / filename).read_text()
    match = re.search(r"(?m)^( *)rule " + re.escape(name) + r":\s*$", source)
    indent = match.group(1)
    block = source[match.end():]
    following = re.search(r"(?m)^" + indent + r"rule \w+:\s*$", block)
    return block[:following.start()] if following else block


def shell_from(block):
    return ast.literal_eval(re.search(r'(?ms)^\s*shell:\s*(""".*?""")', block).group(1))


class JavaRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="trios-java-")
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.bindir = self.work / "bin"
        self.bindir.mkdir()
        self.trace = self.work / "trace.jsonl"
        self.output = self.work / "result.txt"
        self.env = dict(os.environ, PATH=str(self.bindir) + os.pathsep + os.environ["PATH"],
                        TOOL_TRACE=str(self.trace))
        for key in ("CONDA_PREFIX", "JAVA_HOME", "MOCK_MODE"):
            self.env.pop(key, None)
        self.env.update({key: "-Xmx999g -Dunwanted=true" for key in JAVA_OPTIONS})

    def install(self, directory, name, source):
        directory.mkdir(parents=True, exist_ok=True)
        executable = directory / name
        executable.write_text("#!" + sys.executable + "\n" + source)
        executable.chmod(0o700)
        return executable

    def install_gatk(self, directory):
        return self.install(directory, "gatk", r'''import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
args = sys.argv[1:]
java = str(Path(os.environ["JAVA_HOME"]) / "bin" / "java") if os.environ.get("JAVA_HOME") else shutil.which("java")
with open(os.environ["TOOL_TRACE"], "a") as handle:
    handle.write(json.dumps({"tool": "gatk", "executable": sys.argv[0], "args": args,
                           "JAVA_HOME": os.environ.get("JAVA_HOME"), "java": java,
                           "injected": {key: os.environ.get(key) for key in
                                        ("JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS")}}) + "\n")
subprocess.run([java], check=True)
Path(args[args.index("-O") + 1]).write_text("combined\n")
''')

    def install_java(self, directory):
        return self.install(directory, "java", r'''import json
import os
import sys
with open(os.environ["TOOL_TRACE"], "a") as handle:
    handle.write(json.dumps({"tool": "java", "executable": sys.argv[0]}) + "\n")
''')

    def run_gatk(self, use_conda):
        block = rule_block("gatk_hc.smk", "gatk_combine_gvcf")
        # Evaluate the rule's actual deployment selector, not a test-only switch.
        expression = re.search(r"(?m)^\s*use_conda=(.*)$", block).group(1)
        selector = eval(expression, {"workflow": SimpleNamespace(use_conda=use_conda)})
        command = ShellFormatter().format(
            shell_from(block), params=SimpleNamespace(use_conda=selector,
                                                     java_heap_mb=12800, v="input.g.vcf.gz"),
            input=SimpleNamespace(ref="ref.fa"), output=str(self.output),
        )
        return subprocess.run(["bash", "-euo", "pipefail", "-c", command], cwd=self.work,
                              env=self.env, text=True, capture_output=True, check=False)

    def assert_gatk(self, result, gatk, java, java_home=None):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        calls = [json.loads(row) for row in self.trace.read_text().splitlines()]
        self.assertEqual([call["tool"] for call in calls], ["gatk", "java"])
        self.assertEqual(calls[0]["executable"], str(gatk))
        self.assertEqual(calls[0]["java"], str(java))
        self.assertEqual(calls[1]["executable"], str(java))
        self.assertEqual(calls[0]["JAVA_HOME"], str(java_home) if java_home else None)
        self.assertEqual(calls[0]["injected"], dict.fromkeys(JAVA_OPTIONS))
        self.assertEqual(calls[0]["args"][:3], ["--java-options", "-Xmx12800m", "CombineGVCFs"])
        self.assertEqual(self.output.read_text(), "combined\n")

    def test_path_only_runtime_needs_no_conda_prefix_or_java_home(self):
        gatk = self.install_gatk(self.bindir)
        java = self.install_java(self.bindir)
        self.assert_gatk(self.run_gatk(False), gatk, java)

    def test_native_runtime_preserves_its_configured_java_home(self):
        gatk = self.install_gatk(self.bindir)
        native_java = self.work / "native java"
        java = self.install_java(native_java / "bin")
        self.env["JAVA_HOME"] = str(native_java)
        self.assert_gatk(self.run_gatk(False), gatk, java, native_java)

    def test_non_conda_runtime_ignores_an_unrelated_controller_conda_prefix(self):
        gatk = self.install_gatk(self.bindir)
        java = self.install_java(self.bindir)
        self.env["CONDA_PREFIX"] = str(self.work / "unrelated controller conda")
        self.assert_gatk(self.run_gatk(False), gatk, java)

    def test_conda_runtime_uses_its_own_gatk_and_java_despite_host_environment(self):
        prefix = self.work / "selected conda"
        gatk = self.install_gatk(prefix / "bin")
        java = self.install_java(prefix / "bin")
        self.install(self.bindir, "gatk", "raise SystemExit('wrong host GATK')\n")
        self.install(self.bindir, "java", "raise SystemExit('wrong host Java')\n")
        self.env.update(CONDA_PREFIX=str(prefix), JAVA_HOME=str(self.work / "wrong host java"))
        self.assert_gatk(self.run_gatk(True), gatk, java, prefix)

    def test_conda_runtime_requires_activation(self):
        self.install_gatk(self.bindir)
        result = self.run_gatk(True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("The GATK Conda environment must be activated", result.stderr)
        self.assertFalse(self.output.exists())

    def test_conda_runtime_cannot_fall_back_to_host_gatk(self):
        self.install_gatk(self.bindir)
        self.env["CONDA_PREFIX"] = str(self.work / "missing conda")
        result = self.run_gatk(True)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.trace.exists())
        self.assertFalse(self.output.exists())

    def test_picard_rules_override_launcher_defaults_including_small_initial_heap(self):
        self.install(self.bindir, "picard", PICARD_LAUNCHER)
        for filename, name, program in (
            ("bam_qc.smk", "collectwgsmetrics", "CollectWgsMetrics"),
            ("bam_input.smk", "replace_rg", "AddOrReplaceReadGroups"),
        ):
            for allocation in (256, 1024, 16000):
                with self.subTest(rule=name, mem_mb=allocation):
                    block = rule_block(filename, name)
                    expression = re.search(r"(?m)\bheap_mb=(lambda .*)$", block).group(1)
                    heap = eval(expression)(None, SimpleNamespace(mem_mb=allocation))
                    command = ShellFormatter().format(
                        shell_from(block), params=SimpleNamespace(heap_mb=heap),
                        input=(SimpleNamespace(bam="input reads.bam", ref="reference genome.fa")
                               if name == "collectwgsmetrics" else "input reads.bam"),
                        output=str(self.output), wildcards=SimpleNamespace(subj="sample"),
                    )
                    result = subprocess.run(["bash", "-euo", "pipefail", "-c", command],
                                            cwd=self.work, env=self.env, text=True,
                                            capture_output=True, check=False)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    call = json.loads(self.trace.read_text().splitlines()[-1])
                    self.assertEqual(call["program_args"][0], program)
                    self.assertIn("I=input reads.bam", call["program_args"])
                    self.assertEqual(call["max_heap_bytes"], heap * 1024 ** 2)
                    self.assertLess(call["max_heap_bytes"], allocation * 1024 ** 2)
                    self.assertEqual(call["initial_heap_bytes"], 1024 ** 2)
                    self.assertLessEqual(call["initial_heap_bytes"], call["max_heap_bytes"])
                    self.assertEqual(call["injected"], dict.fromkeys(JAVA_OPTIONS))
                    self.assertIn("-XX:+UseSerialGC", call["java_options"])
                    self.assertEqual(self.output.read_text(), "metrics\n")


if __name__ == "__main__":
    unittest.main()
