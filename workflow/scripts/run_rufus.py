#!/usr/bin/env python3
"""Run one pedigree-defined RUFUS trio and validate/publish its final VCF.

Only the Python standard library is required. Caller intermediates and the raw
VCF are retained in a new work directory on every invocation.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import tempfile


CLASSES = ("SNV", "MNV", "INDEL", "SYMBOLIC", "BND", "MIXED", "OTHER")
BWA_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
FINAL_NAME = "child.bam.generator.V2.overlap.hashcount.fastq.bam.FINAL.vcf"
VCF_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


def trio_from_ped(ped, family, samples, bams):
    """Return child/father/mother using PED identifiers, never iteration order."""
    if len(samples) != len(bams) or len(set(samples)) != len(samples):
        raise ValueError("Sample/BAM lists must have equal lengths and unique sample IDs")
    mapping = dict(zip(samples, bams))
    members = {}
    with open(ped) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                raise ValueError("Malformed PED row {}".format(number))
            if fields[0] == family:
                if fields[1] in members:
                    raise ValueError("Duplicate sample in PED family: " + fields[1])
                members[fields[1]] = fields
    children = [row for row in members.values() if row[2] != "0" and row[3] != "0"]
    if len(children) != 1:
        raise ValueError("RUFUS requires exactly one child with both parents in family " + family)
    child = children[0]
    ids = (child[1], child[2], child[3])
    if len(set(ids)) != 3 or any(sample not in members or sample not in mapping for sample in ids):
        raise ValueError("The child and both distinct parents must be present in PED and BAM mapping")
    trio = []
    for role, sample in zip(("child", "father", "mother"), ids):
        bam = Path(mapping[sample]).expanduser().absolute()
        if bam.suffix.lower() != ".bam" or not bam.is_file():
            raise ValueError("RUFUS integration requires a BAM file for " + sample)
        trio.append({"role": role, "sample": sample, "bam": str(bam)})
    if len({Path(member["bam"]).resolve() for member in trio}) != 3:
        raise ValueError("The three trio members must use distinct BAM files")
    return trio


def validate_bam_sample(path, expected):
    """Read only the BGZF/gzip BAM text header to check its read-group SM IDs."""
    with gzip.open(path, "rb") as handle:
        if handle.read(4) != b"BAM\x01":
            raise ValueError("Not a BAM file: " + str(path))
        length_bytes = handle.read(4)
        if len(length_bytes) != 4:
            raise ValueError("Truncated BAM header: " + str(path))
        length = struct.unpack("<i", length_bytes)[0]
        if not 0 <= length <= 64 * 1024 * 1024:
            raise ValueError("Invalid BAM text header length: " + str(path))
        raw = handle.read(length)
        if len(raw) != length:
            raise ValueError("Truncated BAM text header: " + str(path))
    header = raw.decode("utf-8")
    sample_ids = {field[3:] for line in header.splitlines() if line.startswith("@RG\t")
                  for field in line.split("\t")[1:] if field.startswith("SM:")}
    if sample_ids != {expected}:
        raise ValueError("BAM sample mismatch for {}: expected {!r}, found {}".format(
            path, expected, sorted(sample_ids)))


def find_bam_index(bam):
    for location in (bam, bam.resolve()):
        candidates = (Path(str(location) + ".bai"), location.with_suffix(".bai"),
                      Path(str(location) + ".csi"), location.with_suffix(".csi"))
        for path in candidates:
            if path.is_file():
                return path.absolute()
    raise ValueError("BAM index (.bai or .csi) is missing for " + str(bam))


def stage_inputs(trio, reference, work):
    """Use short distinct names because upstream builds commands from basenames."""
    sources = []
    for member in trio:
        bam = Path(member["bam"])
        validate_bam_sample(bam, member["sample"])
        index = find_bam_index(bam)
        staged = work / (member["role"] + ".bam")
        staged.symlink_to(bam.resolve(strict=True))
        Path(str(staged) + index.suffix).symlink_to(index.resolve())
        member["staged_bam"] = staged.name
        member["index"] = str(index)
        sources.extend((bam.resolve(), index.resolve()))
    reference = Path(reference).expanduser().absolute()
    (work / "reference.fa").symlink_to(reference.resolve(strict=True))
    sources.append(reference.resolve())
    for suffix in (".fai",) + BWA_SUFFIXES:
        source = Path(str(reference) + suffix)
        # Accommodate a classic BWA index built with the extension-free prefix.
        if suffix != ".fai" and not source.is_file():
            source = reference.with_suffix(suffix)
        if not source.is_file():
            raise ValueError("Missing reference index: " + str(reference) + suffix)
        source = source.resolve()
        (work / ("reference.fa" + suffix)).symlink_to(source)
        sources.append(source)
    return sources


def string_list(settings, key):
    value = settings.get(key, [])
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise ValueError("rufus.{} must be a list of strings".format(key))
    return value


def positive_integer(value, label, minimum=1):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError("{} must be an integer >= {}".format(label, minimum))
    return value


def executable_path(executable):
    if not isinstance(executable, str) or not executable:
        raise ValueError("Executable must be a nonempty string")
    resolved = shutil.which(os.path.expanduser(executable))
    if resolved is None:
        raise ValueError("Executable is not available: " + executable)
    return str(Path(resolved).resolve())


def caller_command(settings, threads, work, sources):
    threads = positive_integer(threads, "threads", minimum=3)
    kmer = positive_integer(settings.get("kmer_size", 25), "rufus.kmer_size")
    extra = string_list(settings, "extra_args")
    # Core arguments define the biological comparison; do not allow duplicates
    # in the escape hatch to silently replace them or trigger upstream cleanup.
    reserved = {"-s", "--subject", "-c", "--controls", "-r", "--ref", "-t",
                "--threads", "-k", "--kmersize", "-CLEAN"}
    if any(arg.split("=", 1)[0] in reserved for arg in extra):
        raise ValueError("rufus.extra_args cannot override subject, controls, reference, threads, kmer size or cleanup")
    executable = settings.get("executable", "runRufus.sh")
    if not isinstance(executable, str) or not executable:
        raise ValueError("rufus.executable must be a nonempty string")
    container = settings.get("container")
    if container:
        if not isinstance(container, str):
            raise ValueError("rufus.container must be a path or container URI")
        runtime = executable_path(settings.get("runtime", "apptainer"))
        if "://" not in container:
            container = str(Path(container).expanduser().resolve(strict=True))
        command = [runtime, "exec"] + string_list(settings, "runtime_args")
        for parent in sorted({source.parent for source in sources}):
            if any(character in str(parent) for character in ",:\n"):
                raise ValueError("Container bind paths cannot contain comma, colon or newline: " + str(parent))
            command.extend(("--bind", str(parent) + ":" + str(parent) + ":ro"))
        if any(character in str(work) for character in ",:\n"):
            raise ValueError("Container work path cannot contain comma, colon or newline")
        command.extend(("--bind", str(work) + ":" + str(work) + ":rw",
                        "--pwd", str(work), container, executable))
    else:
        command = [executable_path(executable)]
    return command + ["-s", "child.bam", "-c", "father.bam", "-c", "mother.bam",
                      "-r", "reference.fa", "-t", str(threads), "-k", str(kmer)] + extra


def allele_class(ref, alt):
    if "[" in alt or "]" in alt or (len(alt) > 1 and (alt.startswith(".") or alt.endswith("."))):
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return "SYMBOLIC"
    if not re.fullmatch("[ACGTNacgtn]+", alt):
        return "OTHER"
    if len(ref) != len(alt):
        return "INDEL"
    return "SNV" if len(ref) == 1 else "MNV"


def validate_vcf(path, output_handle=None, trio=None):
    """Stream-check structure and counts, optionally copying with PED labels.

    This validates VCF structure and key fields, not caller biological validity.
    Empty (header-only) final VCFs are valid completed outputs.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    header = None
    counts = Counter({key: 0 for key in CLASSES})
    svtypes = Counter()
    total = 0
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if number == 1 and not re.fullmatch(r"##fileformat=VCFv4\.\d+", line.rstrip("\r\n")):
                raise ValueError("Final VCF must begin with a VCFv4 fileformat header")
            if line.startswith("##"):
                if header is not None:
                    raise ValueError("VCF metadata after column header")
            elif line.startswith("#CHROM\t"):
                if header is not None or fields[:8] != VCF_COLUMNS:
                    raise ValueError("Invalid or duplicated VCF column header")
                if len(fields) == 9 or (len(fields) > 8 and fields[8] != "FORMAT"):
                    raise ValueError("Invalid VCF FORMAT/sample columns")
                if len(set(fields[9:])) != len(fields[9:]) or any(not sample for sample in fields[9:]):
                    raise ValueError("Empty or duplicate VCF sample IDs")
                header = fields
                if trio is not None:
                    expected = [member["staged_bam"] for member in trio]
                    # RUFUS emits './child.bam' but unprefixed control names.
                    observed = [name[2:] if name.startswith("./") else name for name in fields[9:]]
                    if observed != expected:
                        raise ValueError("Unexpected RUFUS VCF sample order: expected {}, found {}".format(expected, observed))
                    line = "\t".join(fields[:9] + [member["sample"] for member in trio]) + "\n"
            else:
                if header is None or line.startswith("#") or len(fields) != len(header):
                    raise ValueError("Malformed VCF record at line {}".format(number))
                if any(not value or any(c.isspace() for c in value) for value in fields):
                    raise ValueError("Empty or whitespace-containing VCF field at line {}".format(number))
                if not fields[1].isdigit() or int(fields[1]) < 1 or not re.fullmatch("[ACGTNacgtn]+", fields[3]):
                    raise ValueError("Invalid VCF position or REF at line {}".format(number))
                if fields[5] != ".":
                    try:
                        quality = float(fields[5])
                    except ValueError:
                        raise ValueError("Invalid VCF QUAL at line {}".format(number)) from None
                    if not math.isfinite(quality) or quality < 0:
                        raise ValueError("Invalid VCF QUAL at line {}".format(number))
                if any(not alt for alt in fields[4].split(",")):
                    raise ValueError("Empty VCF ALT at line {}".format(number))
                classes = {allele_class(fields[3], alt) for alt in fields[4].split(",")}
                counts[next(iter(classes)) if len(classes) == 1 else "MIXED"] += 1
                for item in fields[7].split(";"):
                    if item.startswith("SVTYPE="):
                        svtypes[item.split("=", 1)[1]] += 1
                        break
                total += 1
            if output_handle is not None:
                output_handle.write(line)
    if header is None:
        raise ValueError("Final VCF is missing its column header")
    return {"total_records": total, "class_counts": dict(counts), "svtype_counts": dict(sorted(svtypes.items()))}


def publish_vcf(source, output, trio):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=output.parent,
                                         prefix="." + output.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            counts = validate_vcf(source, handle, trio)
        os.replace(temporary, output)
        return counts
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def run(args):
    settings = json.loads(args.settings_json)
    if not isinstance(settings, dict):
        raise ValueError("--settings-json must be a JSON object")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.family) or args.family in (".", ".."):
        raise ValueError("Family ID must contain only letters, digits, underscore, dot or hyphen")
    trio = trio_from_ped(args.ped, args.family, args.samples, args.bams)
    base = Path(args.work_dir).expanduser().absolute() / args.family
    base.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="attempt-", dir=base))
    manifest = {"family": args.family, "trio": trio, "reference": str(Path(args.reference).resolve()),
                "work_dir": str(work), "output": str(Path(args.output).absolute()),
                "started_utc": utc_now(), "status": "started"}
    manifest_path = work / "run.json"
    print("RUFUS work directory: " + str(work), flush=True)
    try:
        sources = stage_inputs(trio, args.reference, work)
        command = caller_command(settings, args.threads, work, sources)
        manifest["command"] = command
        manifest["caller_log"] = str(work / "caller.log")
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        environment = os.environ.copy()
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            environment[key] = "1"
        with (work / "caller.log").open("w") as handle:
            subprocess.run(command, cwd=work, stdout=handle, stderr=subprocess.STDOUT,
                           env=environment, check=True)
        final = work / (FINAL_NAME + ".gz")
        if not final.is_file():
            final = work / FINAL_NAME
        if not final.is_file():
            raise ValueError("RUFUS exited successfully but did not create " + FINAL_NAME + "[.gz]")
        manifest["counts"] = publish_vcf(final, args.output, trio)
        manifest["raw_final_vcf"] = str(final)
        manifest["status"] = "complete"
        print("Published {} records to {}".format(manifest["counts"]["total_records"], args.output), flush=True)
        return 0
    except Exception as error:
        manifest["status"] = "failed"
        manifest["error"] = "{}: {}".format(type(error).__name__, error)
        raise
    finally:
        manifest["finished_utc"] = utc_now()
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def summary(args):
    if len(args.families) != len(args.vcfs) or len(set(args.families)) != len(args.families):
        raise ValueError("Summary requires one VCF per unique family")
    rows = []
    for family, path in zip(args.families, args.vcfs):
        counts = validate_vcf(path)
        rows.append([family, str(counts["total_records"])] +
                    [str(counts["class_counts"][key]) for key in CLASSES] +
                    [json.dumps(counts["svtype_counts"], sort_keys=True, separators=(",", ":"))])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("\t".join(["Family", "TotalRecords"] + list(CLASSES) + ["SVTYPE_counts"]) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    caller = subparsers.add_parser("run")
    caller.add_argument("--ped", required=True)
    caller.add_argument("--family", required=True)
    caller.add_argument("--samples", nargs="+", required=True)
    caller.add_argument("--bams", nargs="+", required=True)
    caller.add_argument("--reference", required=True)
    caller.add_argument("--threads", type=int, required=True)
    caller.add_argument("--settings-json", default="{}")
    caller.add_argument("--work-dir", required=True)
    caller.add_argument("--output", required=True)
    caller.set_defaults(func=run)
    summarizer = subparsers.add_parser("summary")
    summarizer.add_argument("--families", nargs="*", required=True)
    summarizer.add_argument("--vcfs", nargs="*", required=True)
    summarizer.add_argument("--output", required=True)
    summarizer.set_defaults(func=summary)
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, OSError, EOFError, subprocess.CalledProcessError) as error:
        parser.exit(1, "RUFUS integration failed: {}\n".format(error))


if __name__ == "__main__":
    sys.exit(main())
