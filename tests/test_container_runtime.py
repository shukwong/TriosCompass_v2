"""Fail-closed behavior of optional scientific container execution."""
import importlib.util
from pathlib import Path
import shlex
from types import SimpleNamespace
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/container_runtime.py"
spec = importlib.util.spec_from_file_location("container_runtime", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ContainerRuntime = module.ContainerRuntime


class ContainerRuntimeTests(unittest.TestCase):
    def strict(self, config=None, **kwargs):
        settings = {"runtime": {"container_only": True},
                    "containers": {"tabix": "/images/hts tools.sif"},
                    "parabricks": {"singularity_cmd": "apptainer exec --cleanenv --nv /images/parabricks.sif"}}
        if config:
            settings.update(config)
        flags = {"use_singularity": True, "use_conda": False, "singularity_args": " --cleanenv"}
        flags.update(kwargs)
        return ContainerRuntime(settings, **flags)

    def test_legacy_mode_retains_conda_and_native_auxiliary_tools(self):
        runtime = ContainerRuntime({}, use_conda=True)
        self.assertIsNone(runtime.container("bcftools"))
        self.assertEqual(runtime.exec_command("tabix", "bgzip"), "bgzip")
        self.assertEqual(runtime.container("samtools", "docker://samtools:test"), "docker://samtools:test")

    def test_strict_mode_rejects_host_fallback_or_conda_activation(self):
        for flag, value, message in (("use_singularity", False, "use-singularity"),
                                      ("use_conda", True, "use-conda"),
                                      ("singularity_args", "", "cleanenv")):
            with self.subTest(flag=flag), self.assertRaisesRegex(ValueError, message):
                self.strict(**{flag: value})
        with self.assertRaisesRegex(ValueError, "runtime.args"):
            self.strict({"runtime": {"container_only": True, "args": []}})

    def test_profile_gate_enforces_strict_mode_without_runtime_section(self):
        with self.assertRaisesRegex(ValueError, "use-singularity"):
            ContainerRuntime({"container_only": True})
        runtime = self.strict({"container_only": True, "runtime": {"container_only": False}})
        self.assertTrue(runtime.strict)

    def test_missing_or_invalid_enabled_image_fails(self):
        with self.assertRaisesRegex(ValueError, "containers.bcftools"):
            self.strict().container("bcftools")
        with self.assertRaisesRegex(ValueError, "containers.tabix"):
            self.strict({"containers": {"tabix": ["wrong"]}}).container("tabix")

    def test_explicit_exec_preserves_paths_as_single_arguments(self):
        runtime = self.strict({"runtime": {"container_only": True, "executable": "/apps/apptainer",
                                           "args": ["--cleanenv", "--bind", "/data/path with spaces:/data/path with spaces"]}})
        tokens = shlex.split(runtime.exec_command("tabix", "bgzip"))
        self.assertEqual(tokens, ["/apps/apptainer", "exec", "--cleanenv", "--bind",
                                 "/data/path with spaces:/data/path with spaces", "/images/hts tools.sif", "bgzip"])

    def test_manual_gpu_launcher_must_be_container_exec(self):
        for command in ("pbrun", "bash run_parabricks.sh", "apptainer exec --nv image.sif",
                        "apptainer exec --cleanenv image.sif ; pbrun"):
            with self.subTest(command=command), self.assertRaisesRegex(ValueError, "parabricks.singularity_cmd"):
                self.strict({"parabricks": {"singularity_cmd": command}}).parabricks_command()
        command = "singularity exec --cleanenv --nv '/images/gpu tools.sif'"
        self.assertEqual(shlex.split(self.strict({"parabricks": {"singularity_cmd": command}}).parabricks_command()), shlex.split(command))

    @staticmethod
    def rule(name, *, image=None, norun=False, shell="echo test", script=None, wrapper=None):
        return SimpleNamespace(name=name, container_img=image, norun=norun, shellcmd=shell,
                               script=script, wrapper=wrapper)

    def test_audit_rejects_uncontained_new_rule_and_python_run_block(self):
        runtime = self.strict()
        with self.assertRaisesRegex(ValueError, "new_tool"):
            runtime.validate_rules([self.rule("new_tool")])
        with self.assertRaisesRegex(ValueError, "shell/script"):
            runtime.validate_rules([self.rule("python_tool", image="python.sif", shell=None)])

    def test_audit_allows_aggregate_and_manual_gpu_rules_without_nesting(self):
        runtime = self.strict()
        runtime.validate_rules([self.rule("all", norun=True), self.rule("deepvariant_pb"),
                                self.rule("contained_tool", image="tool.sif")])
        with self.assertRaisesRegex(ValueError, "nest runtimes"):
            runtime.validate_rules([self.rule("deepvariant_pb", image="wrong.sif")])


if __name__ == "__main__":
    unittest.main()
