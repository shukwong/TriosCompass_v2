# Controller metadata stays in Snakemake; scientific commands use selected images.
import sys

sys.path.insert(0, os.path.dirname(workflow.source_path("../scripts/container_runtime.py")))
from container_runtime import ContainerRuntime

_tc_runtime = ContainerRuntime(
    config, use_singularity=workflow.use_singularity,
    use_conda=workflow.use_conda, singularity_args=workflow.singularity_args,
)
TC_CONTAINER_ONLY = _tc_runtime.strict
tc_container = _tc_runtime.container
tc_exec = _tc_runtime.exec_command
tc_parabricks = _tc_runtime.parabricks_command
tc_validate_rules = _tc_runtime.validate_rules
