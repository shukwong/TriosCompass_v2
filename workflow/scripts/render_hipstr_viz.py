#!/usr/bin/env python3
"""Render a HipSTR VCF record using its STR coordinates, preserving its VCF identity."""

import argparse
import gzip
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import subprocess
import tempfile


def resolve_locus(vcf, chrom, vcf_pos):
    """Find exactly one VCF record; POS can precede the repeat's INFO/START."""
    opener = gzip.open if str(vcf).endswith(".gz") else open
    matches = []
    with opener(vcf, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0] != chrom or fields[1] != str(vcf_pos):
                continue
            if len(fields) < 8:
                raise ValueError("Malformed VCF record at {}:{}".format(chrom, vcf_pos))
            matches.append(fields)
    if len(matches) != 1:
        raise ValueError("Expected one VCF record at {}:{}, found {}".format(
            chrom, vcf_pos, len(matches)))

    coordinates = {}
    for item in matches[0][7].split(";"):
        key, separator, value = item.partition("=")
        if key not in ("START", "END"):
            continue
        if key in coordinates or not separator or not value.isdigit():
            raise ValueError("Invalid or ambiguous INFO/{} at {}:{}".format(
                key, chrom, vcf_pos))
        coordinates[key] = int(value)
    if "START" not in coordinates or "END" not in coordinates:
        raise ValueError("HipSTR record {}:{} requires INFO/START and INFO/END".format(
            chrom, vcf_pos))
    start, end = coordinates["START"], coordinates["END"]
    if start < 1 or end < start:
        raise ValueError("Invalid STR interval {}:{}-{}".format(chrom, start, end))
    return start, end


class AlignmentHTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.table_classes = set()
        self.reference_cells = 0
        self.read_cells = 0
        self.captions = []
        self.caption = None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.table_classes.update(dict(attrs).get("class", "").split())
        if tag == "td":
            classes = set(dict(attrs).get("class", "").split())
            self.reference_cells += bool(classes & {
                "Areftd", "Creftd", "Greftd", "Treftd", "Nreftd", "vreftd"})
            self.read_cells += bool(classes & {
                "Atd", "Ctd", "Gtd", "Ttd", "Ntd", "-td", "vtd", "inserttd", "snptd"})
        if tag == "caption":
            self.caption = []

    def handle_data(self, data):
        if self.caption is not None:
            self.caption.append(data)

    def handle_endtag(self, tag):
        if tag == "caption" and self.caption is not None:
            self.captions.append("".join(self.caption).split())
            self.caption = None


def validate_alignment_html(text, chrom, start, end):
    """Accept HipSTR HTML fragments, but never a stylesheet-only or wrong-locus view."""
    parsed = AlignmentHTML()
    parsed.feed(text)
    if "reftable" not in parsed.table_classes or not parsed.reference_cells or not parsed.read_cells:
        raise ValueError("VizAln produced no reference/read alignment content")
    expected = [chrom, str(start), str(end)]
    if expected not in parsed.captions:
        raise ValueError("VizAln output does not contain STR caption {}:{}-{}".format(
            chrom, start, end))


def render(vcf, alignments, output, chrom, vcf_pos):
    if vcf_pos < 1:
        raise ValueError("VCF position must be positive")
    start, end = resolve_locus(vcf, chrom, vcf_pos)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".vizaln-", dir=str(output.parent)) as work:
        staged = Path(work) / "alignment.html"
        subprocess.run(["VizAln_rev", str(alignments), str(staged), chrom, str(start)],
                       check=True)
        validate_alignment_html(staged.read_text(), chrom, start, end)
        os.replace(str(staged), str(output))
    return {"chrom": chrom, "vcf_pos": vcf_pos, "str_start": start,
            "str_end": end, "output": str(output)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--alignments", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--vcf-pos", type=int, required=True)
    args = parser.parse_args()
    try:
        result = render(args.vcf, args.alignments, args.output, args.chrom, args.vcf_pos)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, "HipSTR visualization failed: {}\n".format(error))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
