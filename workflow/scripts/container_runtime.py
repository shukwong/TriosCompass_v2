"""Container selection and fail-closed execution checks for TriosCompass."""

import os
import shlex


class ContainerRuntime:
    """Keep legacy Conda optional; require explicit containers in strict mode."""

    MANUAL_RULES = frozenset({
        "fq2bam", "deepvariant_pb", "gatkhc_pb", "gatk_genotype_gvcf_pb",
        "collectmultiplemetrics",
    })

    def __init__(self, config, *, use_singularity=False, use_conda=False,
                 singularity_args=""):
        self.config = config
        self.settings = config.get("runtime", {})
        self.images = config.get("containers", {})
        if not isinstance(self.settings, dict) or not isinstance(self.images, dict):
            raise ValueError("runtime and containers must be configuration mappings")
        # Profiles set the top-level gate so legacy configurations cannot
        # silently use host tools when the container-only profile is selected.
        gates = (self.settings.get("container_only", False), config.get("container_only", False))
        if not all(isinstance(gate, bool) for gate in gates):
            raise ValueError("container_only and runtime.container_only must be true or false")
        self.strict = any(gates)
        self.executable = self.settings.get("executable", "apptainer")
        if not isinstance(self.executable, str) or os.path.basename(self.executable) not in ("apptainer", "singularity"):
            raise ValueError("runtime.executable must name apptainer or singularity")
        self.args = self._tokens(self.settings.get("args", ["--cleanenv"]), "runtime.args")
        if self.strict:
            if not use_singularity:
                raise ValueError("runtime.container_only requires --use-singularity")
            if use_conda:
                raise ValueError("runtime.container_only requires use-conda: false (omit --use-conda)")
            if "--cleanenv" not in self._tokens(singularity_args or "", "singularity_args"):
                raise ValueError("runtime.container_only requires --singularity-args ' --cleanenv'")
            if "--cleanenv" not in self.args:
                raise ValueError("runtime.container_only requires --cleanenv in runtime.args")

    @staticmethod
    def _tokens(value, label):
        if isinstance(value, str):
            return shlex.split(value)
        if isinstance(value, list) and all(isinstance(item, str) for item in value):
            return value
        raise ValueError(label + " must be a string or a list of strings")

    def container(self, name, default=None):
        image = self.images.get(name, default)
        if image is None or image == "":
            if self.strict:
                raise ValueError("runtime.container_only requires containers." + name)
            return None
        if not isinstance(image, str) or not image.strip():
            raise ValueError("containers." + name + " must be a nonempty image path or URI")
        return image

    def exec_command(self, name, executable):
        """Run auxiliary tools beside a manual GPU launch, without nesting runtimes."""
        image = self.container(name)
        if image is None:
            return shlex.quote(executable)
        return shlex.join([self.executable, "exec", *self.args, image, executable])

    def parabricks_command(self):
        command = self.config.get("parabricks", {}).get("singularity_cmd", "")
        if not isinstance(command, str):
            raise ValueError("parabricks.singularity_cmd must be a command string")
        if self.strict:
            tokens = shlex.split(command)
            if (len(tokens) < 4 or os.path.basename(tokens[0]) not in ("apptainer", "singularity")
                    or tokens[1] != "exec" or "--cleanenv" not in tokens[2:]):
                raise ValueError("runtime.container_only requires parabricks.singularity_cmd to start "
                                 "with apptainer/singularity exec and include --cleanenv")
            # This is a command prefix, not an arbitrary shell program. Quoting is
            # retained below so spaces in image/bind paths remain one argument.
            if any(c in command for c in "\n\r") or any(t in (";", "&&", "||", "|", ">", "<") for t in tokens):
                raise ValueError("parabricks.singularity_cmd must be one container exec command")
            return shlex.join(tokens)
        return command

    def validate_rules(self, rules):
        """Reject newly introduced host commands when container-only mode is selected."""
        if not self.strict:
            return
        missing = []
        for rule in rules:
            if rule.norun:
                continue
            if rule.name in self.MANUAL_RULES:
                if rule.container_img:
                    raise ValueError(rule.name + " manually launches Parabricks; a rule container would nest runtimes")
                self.parabricks_command()
                continue
            if not rule.container_img:
                missing.append(rule.name)
            elif not (rule.shellcmd or rule.script or rule.wrapper):
                # Snakemake run: executes controller Python and cannot guarantee
                # that scientific Python is isolated by a container annotation.
                raise ValueError(rule.name + " must use shell/script instead of run in container-only mode")
        if missing:
            raise ValueError("runtime.container_only found executable rules without containers: "
                             + ", ".join(sorted(missing)))
