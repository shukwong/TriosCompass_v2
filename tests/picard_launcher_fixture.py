"""Lightweight model of the pinned Bioconda Picard launcher's JVM argument rules.

The picard:2.27.3--hdfd78af_0 launcher consumes -Xm*, -XX*, and -D* options
before calling Java. If no -Xm* option or _JAVA_OPTIONS is present it appends
-Xms512m -Xmx2g, which take precedence over JAVA_TOOL_OPTIONS. Production image
probes remain separate from this portable command-boundary regression fixture.
"""

PICARD_LAUNCHER = r'''import json
import os
from pathlib import Path
import shlex
import sys

args = sys.argv[1:]
memory = [arg for arg in args if arg.startswith("-Xm")]
properties = [arg for arg in args if arg.startswith(("-XX", "-D"))]
program = [arg for arg in args if not arg.startswith(("-Xm", "-XX", "-D"))]
if not memory and "_JAVA_OPTIONS" not in os.environ:
    memory = ["-Xms512m", "-Xmx2g"]
options = (shlex.split(os.environ.get("JAVA_TOOL_OPTIONS", ""))
           + shlex.split(os.environ.get("JDK_JAVA_OPTIONS", ""))
           + memory + properties
           + shlex.split(os.environ.get("_JAVA_OPTIONS", "")))

def size_bytes(prefix):
    values = [option[len(prefix):] for option in options if option.startswith(prefix)]
    if not values:
        return None
    value = values[-1].lower()
    unit = {"k": 1024, "m": 1024 ** 2, "g": 1024 ** 3}
    return int(value[:-1]) * unit[value[-1]] if value[-1] in unit else int(value)

trace = {"tool": "picard", "args": args, "program_args": program,
         "java_options": options, "initial_heap_bytes": size_bytes("-Xms"),
         "max_heap_bytes": size_bytes("-Xmx"),
         "injected": {key: os.environ.get(key) for key in
                      ("JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS")}}
with open(os.environ["TOOL_TRACE"], "a") as handle:
    handle.write(json.dumps(trace) + "\n")
if trace["initial_heap_bytes"] > trace["max_heap_bytes"]:
    sys.exit("Initial heap size set to a larger value than the maximum heap size")
Path(next(arg[2:] for arg in program if arg.startswith("O="))).write_text("metrics\n")
if os.environ.get("MOCK_MODE") == "fail":
    sys.exit(31)
'''
