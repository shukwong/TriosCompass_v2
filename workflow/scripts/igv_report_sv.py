#!/usr/bin/env python3
"""Build an offline trio SV report with pinned igv-reports and a local IGV.js bundle.

VCF coordinates/fields are retained verbatim. PED + BAM SM headers, never merged
VCF column order, establish child/father/mother identities. All source records
are accounted for in manifest.json, including unsupported or malformed records.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import html
import json
import os
from pathlib import Path
import re
import socket
import sys

SUPPORTED = {"DEL", "DUP", "INV", "INS", "BND"}
IGV_REPORTS_VERSION = "1.17.0"
IGV_JS_VERSION = "3.5.2"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def script_json(value):
    # Data belongs inside a script element, not HTML markup.
    return json.dumps(value).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")


def positive_int(value, field):
    if isinstance(value, bool):
        raise ValueError(field + " must be an integer, not an INFO flag")
    number = int(value)
    if number < 1:
        raise ValueError(field + " must be a positive integer")
    return number


def parse_record(line, number, sample_ids, min_sv_length):
    record = dict(record=number, source_record=line, status="skipped", reason="", notes=[])
    fields = line.split("\t")
    if len(fields) < 8:
        record["reason"] = "VCF record has fewer than eight fields"
        return record
    chrom, pos, source_id, ref, alt, qual, filt, raw_info = fields[:8]
    info = {}
    for item in raw_info.split(";"):
        if item and item != ".":
            key, sep, value = item.partition("=")
            info[key] = value if sep else True
    record.update(chrom=chrom, id=source_id, ref=ref, alt=alt, qual=qual, filter=filt,
                  info=info, source_end=info.get("END"), source_svlen=info.get("SVLEN"), samples={})
    if len(fields) > 8:
        keys = fields[8].split(":")
        for sample, data in zip(sample_ids, fields[9:]):
            record["samples"][sample] = dict(zip(keys, data.split(":")))
    try:
        pos = positive_int(pos, "POS")
        record["pos"] = pos
        if not chrom or ":" in chrom:
            raise ValueError("Missing or unsupported chromosome name")
        if "," in alt:
            raise ValueError("Multiallelic ALT is unsupported; normalize to biallelic records first")
        raw_type = str(info.get("SVTYPE", ""))
        svtype = raw_type.split(":")[0].upper()
        symbolic = re.fullmatch(r"<([^>]+)>", alt)
        if not svtype and symbolic:
            svtype = symbolic.group(1).split(":")[0].upper()
        if "[" in alt or "]" in alt:
            if svtype and svtype != "BND":
                raise ValueError("Breakend ALT conflicts with SVTYPE")
            svtype = "BND"
        sequence_alleles = bool(re.fullmatch(r"[ACGTNacgtn]+", ref) and re.fullmatch(r"[ACGTNacgtn]+", alt))
        if not svtype and sequence_alleles:
            delta = len(alt) - len(ref)
            if abs(delta) >= min_sv_length:
                svtype = "INS" if delta > 0 else "DEL"
                record["notes"].append("SVTYPE inferred from sequence allele lengths")
        record["svtype"] = svtype
        if svtype not in SUPPORTED:
            raise ValueError("Unsupported or small variant: " + (svtype or "no structural SVTYPE"))
        source_length = info.get("SVLEN")
        svlen = None
        if source_length is not None:
            try:
                svlen = int(source_length)
            except (TypeError, ValueError):
                record["notes"].append("Noninteger SVLEN retained verbatim")
        end = positive_int(info["END"], "END") if info.get("END") is not None else None
        panels = [dict(chrom=chrom, pos1=pos)]
        if svtype == "BND":
            match = re.fullmatch(r"([ACGTNacgtn.]*)([\[\]])([^:\[\]]+):(\d+)\2([ACGTNacgtn.]*)", alt)
            if match and bool(match[1]) != bool(match[5]):
                chr2, pos2 = match[3], positive_int(match[4], "breakend position")
                if info.get("CHR2") not in (None, chr2):
                    raise ValueError("BND ALT and CHR2 disagree")
                record["orientation"] = dict(bracket=match[2], remote_first=not bool(match[1]))
            elif alt == "<BND>" and info.get("CHR2") and end:
                chr2, pos2 = info["CHR2"], end
                record["notes"].append("BND partner derived from CHR2/END; orientation unavailable")
            else:
                raise ValueError("Malformed or single-ended BND ALT")
            panels.append(dict(chrom=chr2, pos1=pos2))
        elif svtype in {"DEL", "DUP", "INV"}:
            if info.get("CHR2", chrom) != chrom:
                raise ValueError("Interchromosomal event requires BND representation")
            if end is None:
                if svtype == "DEL" and sequence_alleles and len(ref) > len(alt):
                    end = pos + len(ref) - 1
                    record["notes"].append("Second anchor inferred from source REF span; source END absent (allele may be unnormalized)")
                elif svlen not in (None, 0):
                    end = pos + abs(svlen)
                    record["notes"].append("Second anchor inferred from SVLEN; source END absent")
            if end is None or end <= pos:
                raise ValueError("Missing or invalid second breakpoint (END)")
            panels.append(dict(chrom=chrom, pos1=end))
            if svlen is not None and abs(svlen) != end - pos:
                record["notes"].append("SVLEN differs from END-POS; source values retained")
        elif end not in (None, pos):
            record["notes"].append("INS displayed at POS; source END retained")
        if svtype == "BND":
            length = None
        elif svlen is not None:
            length = abs(svlen)
        elif svtype in {"INS", "DEL"} and sequence_alleles:
            length = abs(len(alt) - len(ref))
        elif end is not None and svtype != "INS":
            length = end - pos
        else:
            length = None
        # Explicit structural types are retained; threshold only classifies
        # unannotated sequence alleles in mixed small-variant/SV Rufus VCFs.
        record.update(end=end, panels=panels, length=length, status="pending")
    except (ValueError, TypeError) as error:
        record["reason"] = str(error)
    return record


def read_variants(filename, min_sv_length=50):
    records, samples = [], []
    opener = gzip.open if str(filename).endswith(".gz") else open
    with opener(filename, "rt") as stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                if len(samples) != len(set(samples)):
                    raise ValueError("VCF sample names must be unique")
            if line and not line.startswith("#"):
                records.append(parse_record(line, len(records) + 1, samples, min_sv_length))
    return records, samples


def trio_order(ped, family, samples, bams):
    if len(samples) != len(bams) or len(set(samples)) != len(samples):
        raise ValueError("--samples and --bams must pair one-to-one with unique sample IDs")
    mapping, members, children = dict(zip(samples, bams)), set(), []
    with open(ped) as stream:
        for line in stream:
            row = line.split()
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 6:
                raise ValueError("PED row has fewer than six fields")
            if row[0] == family:
                if row[1] in members:
                    raise ValueError("Duplicate sample in family PED: " + row[1])
                members.add(row[1])
                if row[2] != "0" and row[3] != "0":
                    children.append(row[1:4])
    if len(children) != 1:
        raise ValueError("Expected exactly one complete trio child in family " + family)
    ids = children[0]
    if len(set(ids)) != 3 or not all(sample in mapping and sample in members for sample in ids):
        raise ValueError("PED child and parents must be distinct and present in PED and BAM mapping")
    if len({str(Path(mapping[sample]).absolute()) for sample in ids}) != 3:
        raise ValueError("The trio must use three distinct BAM files")
    return [dict(sample=sample, role=role, label=role.title() + " " + sample,
                 bam=str(Path(mapping[sample]).absolute()))
            for role, sample in zip(("child", "father", "mother"), ids)]


def validate_inputs(trio, records, reference):
    import pysam
    with pysam.FastaFile(reference) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
    for record in records:
        if record["status"] == "pending":
            for panel in record["panels"]:
                if panel["chrom"] not in lengths or panel["pos1"] > lengths[panel["chrom"]]:
                    record.update(status="skipped", reason="Breakpoint outside reference FASTA dictionary")
                    break
    for member in trio:
        with pysam.AlignmentFile(member["bam"], reference_filename=reference) as bam:
            if not bam.has_index():
                raise ValueError("Alignment index missing: " + member["bam"])
            observed = sorted({rg.get("SM", "") for rg in bam.header.to_dict().get("RG", [])})
            if observed != [member["sample"]]:
                raise ValueError("BAM read-group sample mismatch for {}: {}".format(member["sample"], observed))
            bam_lengths = dict(zip(bam.references, bam.lengths))
            for record in records:
                if record["status"] == "pending":
                    for panel in record["panels"]:
                        if bam_lengths.get(panel["chrom"]) != lengths[panel["chrom"]]:
                            raise ValueError("BAM/reference dictionary mismatch at " + panel["chrom"] + " for " + member["sample"])
            member.update(read_group_samples=observed, size_bytes=Path(member["bam"]).stat().st_size)
    return lengths


def support_fields(record, sample):
    fields = {} if record.get("ignore_vcf_samples", False) else record.get("samples", {}).get(sample, {})
    genotype = fields.get("GT", "NA")
    result = {"GT": "NA" if genotype in (".", "./.", ".|.") else genotype}
    for key in ("PR", "SR"):
        value = fields.get(key, ".")
        parts = value.split(",")
        valid_pair = len(parts) == 2 and all(part == "." or re.fullmatch(r"[0-9]+", part) for part in parts)
        result[key] = " / ".join("NA" if part == "." else part for part in parts) if valid_pair else "NA"
    return result


def extraction_windows(panels, lengths, flank):
    """Return merged 1-based inclusive windows for igv-reports/samtools slicing."""
    windows = []
    for panel in panels:
        windows.append(dict(chr=panel["chrom"], start=max(1, panel["pos1"] - flank),
                            end=min(lengths[panel["chrom"]], panel["pos1"] + flank)))
    if len(windows) == 2 and windows[0]["chr"] == windows[1]["chr"]:
        first, second = windows
        if first["start"] <= second["end"] + 1 and second["start"] <= first["end"] + 1:
            windows = [dict(chr=first["chr"], start=min(first["start"], second["start"]),
                            end=max(first["end"], second["end"]))]
    return windows


def initial_locus(panel, bound, width):
    width = min(width, bound["end"] - bound["start"])
    start = max(bound["start"], min(panel["pos1"] - 1 - width // 2, bound["end"] - width))
    return "{}:{}-{}".format(panel["chrom"], start + 1, start + width)


def slice_alignment(member, regions, args):
    """Keep mate-reference IDs and emit each alignment once across both windows.

    igv-reports 1.17's BamReader prunes non-panel @SQ entries and uses separate
    region iterators. Complete headers preserve off-panel mate identities;
    samtools -M prevents double counting a read spanning both selected windows.
    """
    import pysam
    options = ["-b", "-h", "-M", "-F", str(args.exclude_flags)]
    if member["bam"].endswith(".cram"):
        options += ["-T", args.reference]
    options.append(member["bam"])
    options += ["{}:{}-{}".format(region["chr"], region["start"], region["end"]) for region in regions]
    return pysam.view(*options, catch_stdout=True)


def build_sessions(records, trio, args, lengths):
    import pysam
    from igv_reports.datauri import get_data_uri
    sessions, windows_metadata = {}, {}
    with pysam.FastaFile(args.reference) as fasta:
        for record in records:
            if record["status"] != "pending":
                continue
            if args.max_variants and len(sessions) >= args.max_variants:
                record.update(status="skipped", reason="Configured max_variants limit reached")
                continue
            sid = str(len(sessions))
            regions = extraction_windows(record["panels"], lengths, args.flank)
            fasta_chunks = []
            for region in regions:
                start, end = max(1, region["start"] - 500), min(lengths[region["chr"]], region["end"] + 500)
                sequence = fasta.fetch(region["chr"], start - 1, end)
                fasta_chunks.append(">{}:{}-{} @len={}\n{}".format(region["chr"], start, end, lengths[region["chr"]], sequence))
            bounds = []
            for panel in record["panels"]:
                bound = next(region for region in regions if region["chr"] == panel["chrom"] and region["start"] <= panel["pos1"] <= region["end"])
                bounds.append(dict(chr=bound["chr"], start=bound["start"] - 1, end=bound["end"]))
            locus = " ".join(initial_locus(panel, bound, args.initial_window) for panel, bound in zip(record["panels"], bounds))
            bed = "\n".join("{}\t{}\t{}\tBreakpoint {}".format(panel["chrom"], panel["pos1"] - 1, panel["pos1"], index + 1) for index, panel in enumerate(record["panels"])) + "\n"
            tracks = [dict(name="Breakpoint anchors", type="annotation", format="bed", url=get_data_uri(bed),
                           displayMode="EXPANDED", height=45, color="#8c3b93")]
            for member in trio:
                data = slice_alignment(member, regions, args)
                tracks.append(dict(name=member["label"], type="alignment", format="bam", url=get_data_uri(data),
                                   height=args.track_height, coverageTrackHeight=55, min=0, max=args.coverage_max,
                                   autoscale=False, viewAsPairs=True, showSoftClips=True, showInsertionText=True,
                                   colorBy="unexpectedPair", samplingDepth=args.sampling_depth, samplingWindowSize=100,
                                   filter={"mq": 0}, order=len(tracks)))
            session = dict(locus=locus, reference={"fastaURL": get_data_uri("\n".join(fasta_chunks))}, tracks=tracks)
            sessions[sid] = get_data_uri(json.dumps(session))
            windows_metadata[sid] = dict(initialLocus=locus, bounds=bounds)
            record.update(status="rendered", session_id=sid, initial_locus=locus,
                          extraction_windows_1based_inclusive=regions,
                          trio_support={member["sample"]: support_fields(record, member["sample"]) for member in trio})
    return sessions, windows_metadata


def build_table(records, trio):
    headers = ["unique_id", "SV type", "Breakpoint A", "Breakpoint B", "Length (bp)", "Source END", "Source SVLEN"]
    headers += [member["label"] for member in trio] + ["Filter", "Source ID"]
    rows = []
    for record in records:
        if record["status"] != "rendered":
            continue
        anchors = ["{}:{}".format(panel["chrom"], panel["pos1"]) for panel in record["panels"]]
        values = [record["svtype"], anchors[0], anchors[1] if len(anchors) > 1 else "NA",
                  record["length"], record["source_end"], record["source_svlen"]]
        for member in trio:
            support = record["trio_support"][member["sample"]]
            values.append("GT: {} | PR: {} | SR: {}".format(support["GT"], support["PR"], support["SR"]))
        values += [record["filter"], record["id"]]
        # igv-reports' template inserts cells using innerHTML.
        rows.append([record["session_id"]] + [html.escape(str(value)) if value is not None else "NA" for value in values])
    return dict(headers=[html.escape(header) for header in headers], rows=rows)


def render_document(records, trio, sessions, windows, args, manifest):
    import igv_reports
    title = "Trio SV review — " + args.family
    intro = "<h1>" + html.escape(title) + "</h1>"
    intro += ("<p>Tracks: child → father → mother, identified from the PED and verified BAM read groups. "
              "SV candidates and caller support are not independently validated de novo classifications. "
              "PR and SR show source REF / ALT counts only when the VCF sample ID matches a verified family member. "
              "Caller-summary VCF columns are ignored when configured, even if a label matches a family member. "
              "NA means unavailable, not zero or a reference genotype. Reciprocal breakend records are separate source rows; "
              "their counts are never summed.</p><p>Insertions show one anchor; DEL, DUP, INV and BND show both breakpoints, "
              "each in reference-forward orientation. Long-event interiors are not embedded. Source END/SVLEN are retained "
              "separately from any inferred anchor or display length. Every input record, including skips, is listed below.</p>")
    intro += ("<details><summary>Read display settings</summary><p>Embedded reads extend up to ±{} bp around each anchor, "
              "clipped to the contig boundaries. Extraction excludes SAM flags {} without subsampling. Display sampling: "
              "{} reads per 100-base window; MAPQ threshold 0. Shared coverage scale: 0–{}. "
              "Navigation beyond these local windows cannot retrieve additional reads.</p></details>").format(
                  args.flank, args.exclude_flags, args.sampling_depth, args.coverage_max)
    metadata = {record["session_id"]: record for record in records if record["status"] == "rendered"}
    footer = "<details><summary>All source records and skipped-record reasons</summary><pre>" + html.escape(json.dumps(records, indent=2)) + "</pre></details>"
    footer += "<details><summary>Run provenance</summary><pre>" + html.escape(json.dumps(manifest, indent=2)) + "</pre></details>"
    style = "<style>body{margin:20px;font:13px system-ui,sans-serif}h1{font-size:23px}p{max-width:1300px;line-height:1.5}td,th{padding:5px;white-space:nowrap}th{background:#266a94;color:white}#tableSelectorDiv{overflow:auto}pre{white-space:pre-wrap;overflow-wrap:anywhere;max-height:260px;overflow:auto;font:12px monospace}</style>"
    if not sessions:
        return '<!doctype html><html lang="en"><head><meta charset="utf-8"><title>' + html.escape(title) + "</title>" + style + "</head><body>" + intro + "<p>No supported SV candidates to display.</p>" + footer + "</body></html>"
    if not args.igv_js:
        raise ValueError("--igv-js is required for a nonempty report; supply a local IGV.js 3.5.2 bundle")
    library = Path(args.igv_js).read_text()
    if "createBrowser" not in library or len(library) < 100000 or not re.search(r'return\s*["\']3\.5\.2["\']', library) or re.search(r"</script", library, re.I):
        raise ValueError("--igv-js must be an unmodified IGV.js 3.5.2 distribution bundle")
    template = (Path(igv_reports.__file__).parent / "templates" / "variant_template.html").read_text()
    cdn = '<script src="https://cdn.jsdelivr.net/npm/igv@3.5.2/dist/igv.min.js"></script>'
    if cdn not in template:
        raise ValueError("Unexpected igv-reports template; install pinned igv-reports 1.17.0")
    template = template.replace(cdn, "@EMBED_IGV_JS@")
    template = template.replace('"@TABLE_JSON@"', script_json(build_table(records, trio))).replace('"@SESSION_DICTIONARY@"', script_json(sessions))
    template = template.replace("<title>IGV Variant Inspector </title>", "<title>" + html.escape(title) + "</title>")
    template = template.replace("showChromosomeWidget: false,", "showChromosomeWidget: false, loadDefaultGenomes: false,")
    template = template.replace("</head>", style + "</head>")
    template = template.replace("<body>", "<body>" + intro + '<details><summary>Selected source record</summary><pre id="event-details"></pre></details>', 1)
    controls = Path(args.controls_template).read_text()
    controls = controls.replace("__READ_WINDOW_METADATA__", script_json(windows)).replace("__EVENT_METADATA__", script_json(metadata))
    controls = controls.replace("__READ_HEIGHT__", str(args.read_height)).replace("__TRACK_HEIGHT__", str(args.track_height))
    template = template.replace("                initTable()", "                initTable()\n                initReadControls()", 1)
    keyhook = "document.addEventListener('keydown', function (event) {"
    template = template.replace(keyhook, keyhook + "\n        if (event.target.closest('input,select,textarea,button')) return;", 1)
    template = template.replace("</body>", footer + controls + "</body>")
    return template.replace("@EMBED_IGV_JS@", '<script type="text/javascript">' + library + "</script>")


def write_tables(output, records):
    fields = ["record", "id", "chrom", "pos", "svtype", "source_end", "source_svlen", "length", "filter", "status", "reason"]
    for filename, rows in (("variants.tsv", records), ("skipped_records.tsv", [r for r in records if r["status"] == "skipped"])):
        with (output / filename).open("w") as stream:
            writer = csv.DictWriter(stream, fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("vcf", "ped", "family", "reference", "output-dir"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--bams", nargs="+", required=True)
    parser.add_argument("--ignore-vcf-samples", action="store_true", help="Treat VCF columns as caller summaries; never show them as trio genotype/support")
    parser.add_argument("--controls-template", default=str(Path(__file__).parent / "igv_read_controls.html"), help="Path to read-control template (explicit for Snakemake source caching)")
    parser.add_argument("--igv-js", help="Local IGV.js 3.5.2 bundle, embedded without network access")
    for name, default in (("flank", 2000), ("initial-window", 400), ("track-height", 320), ("read-height", 14),
                          ("coverage-max", 200), ("sampling-depth", 500), ("min-sv-length", 50), ("exclude-flags", 1536), ("max-variants", 0)):
        parser.add_argument("--" + name, type=int, default=default)
    args = parser.parse_args(argv)
    if min(args.flank, args.initial_window, args.coverage_max, args.sampling_depth, args.min_sv_length) <= 0 or min(args.exclude_flags, args.max_variants) < 0:
        parser.error("Sizes and sampling limits must be positive; exclude-flags and max-variants must be nonnegative")
    if not 6 <= args.read_height <= 32 or not 180 <= args.track_height <= 800:
        parser.error("read-height must be 6–32; track-height must be 180–800")
    from importlib.metadata import version
    if version("igv-reports") != IGV_REPORTS_VERSION:
        raise ValueError("Use pinned igv-reports " + IGV_REPORTS_VERSION)
    records, samples = read_variants(args.vcf, args.min_sv_length)
    for record in records:
        record["ignore_vcf_samples"] = args.ignore_vcf_samples
    trio = trio_order(args.ped, args.family, args.samples, args.bams)
    lengths = validate_inputs(trio, records, args.reference)
    if any(record["status"] == "pending" for record in records) and not args.igv_js:
        raise ValueError("--igv-js is required for a nonempty report; supply a local IGV.js 3.5.2 bundle")
    output = Path(args.output_dir).absolute()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory must be empty: " + str(output))
    output.mkdir(parents=True, exist_ok=True)
    sessions, windows = build_sessions(records, trio, args, lengths)
    manifest = dict(family=args.family, created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    host=socket.gethostname(), cwd=os.getcwd(), job_id=os.environ.get("SLURM_JOB_ID"),
                    command=[sys.executable] + sys.argv, python=sys.version,
                    versions={name: version(name) for name in ("igv-reports", "pysam")},
                    igv_js_version=IGV_JS_VERSION if sessions else None,
                    igv_js_sha256=sha256(args.igv_js) if args.igv_js else None,
                    source_vcf=str(Path(args.vcf).absolute()), source_vcf_sha256=sha256(args.vcf),
                    ped=str(Path(args.ped).absolute()), ped_sha256=sha256(args.ped), reference=str(Path(args.reference).absolute()),
                    source_vcf_samples=samples, trio=trio, settings=vars(args),
                    coordinate_convention="VCF POS/END and extraction windows: 1-based inclusive; browser navigation bounds: 0-based half-open; BED anchors: [POS-1,POS)",
                    alignment_extraction="samtools multi-region iterator (-M); complete BAM header retained to preserve off-panel mate references",
                    rendered=len(sessions), skipped=sum(r["status"] == "skipped" for r in records), records=records,
                    read_window_metadata=windows)
    document = render_document(records, trio, sessions, windows, args, {k: v for k, v in manifest.items() if k != "records"})
    (output / "index.html").write_text(document)
    # Both entry points are complete standalone HTMLs, so Snakemake's report
    # embedding and a downloaded file work without relative-resource links.
    (output / "variants.html").write_text(document)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_tables(output, records)
    print("Rendered {} SV records; skipped {}. Report: {}".format(manifest["rendered"], manifest["skipped"], output / "index.html"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
