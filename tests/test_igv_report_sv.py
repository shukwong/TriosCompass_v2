"""Identity, coordinate and offline-payload regressions for trio SV viewing.

Run with unittest. The real BAM/HTML integration test uses the pinned IGV
reports environment and IGV_REPORTS_TEST_JS=/path/to/igv.3.5.2.min.js.
"""
import base64
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/igv_report_sv.py"
spec = importlib.util.spec_from_file_location("igv_report_sv", SCRIPT)
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)


def vcf_line(alt="<DEL>", info="SVTYPE=DEL;END=400;SVLEN=-200", pos=200, ident="v1", fmt=None):
    row = ["chr1", str(pos), ident, "A", alt, ".", "PASS", info]
    if fmt:
        row += fmt
    return "\t".join(row)


def decode(uri):
    data = base64.b64decode(uri.split(",", 1)[1])
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


class VariantTests(unittest.TestCase):
    def test_source_coordinates_are_not_rewritten(self):
        r = sv.parse_record(vcf_line(info="SVTYPE=DEL;END=400;SVLEN=-199"), 1, [], 50)
        self.assertEqual((r["source_end"], r["source_svlen"]), ("400", "-199"))
        self.assertEqual(r["panels"], [{"chrom": "chr1", "pos1": 200}, {"chrom": "chr1", "pos1": 400}])
        self.assertTrue(any("differs" in note for note in r["notes"]))

    def test_all_supported_types_and_breakend_orientations(self):
        for kind in ("DEL", "DUP", "INV", "INS"):
            r = sv.parse_record(vcf_line(alt="<" + kind + ">", info="SVTYPE=" + kind + ";END=400"), 1, [], 50)
            self.assertEqual(r["status"], "pending")
            self.assertEqual(len(r["panels"]), 1 if kind == "INS" else 2)
        for alt in ("A[chr2:700[", "A]chr2:700]", "[chr2:700[A", "]chr2:700]A"):
            r = sv.parse_record(vcf_line(alt=alt, info="SVTYPE=BND"), 1, [], 50)
            self.assertEqual(r["status"], "pending", r)
            self.assertEqual(r["panels"][1], {"chrom": "chr2", "pos1": 700})

    def test_malformed_and_unsupported_records_are_accounted_for(self):
        lines = ["broken", vcf_line(pos=0), vcf_line(info="SVTYPE=DEL;END=bad"),
                 vcf_line(alt="<DUP>,<DEL>"), vcf_line(alt="A[chr2:5]", info="SVTYPE=BND"),
                 vcf_line(alt="A.", info="SVTYPE=BND"), vcf_line(info="SVTYPE=CNV;END=400")]
        for line in lines:
            r = sv.parse_record(line, 1, [], 50)
            self.assertEqual(r["status"], "skipped", r)
            self.assertTrue(r["reason"])
            self.assertEqual(r["source_record"], line)

    def test_mixed_rufus_sequence_alleles(self):
        small = sv.parse_record(vcf_line(alt="AT", info="."), 1, [], 50)
        long = sv.parse_record(vcf_line(alt="A" + "T" * 60, info="."), 2, [], 50)
        self.assertEqual(small["status"], "skipped")
        self.assertEqual((long["svtype"], long["length"], len(long["panels"])), ("INS", 60, 1))
        deletion = sv.parse_record(vcf_line(alt="A", info=".").replace("\tA\tA\t", "\t" + "A" * 61 + "\tA\t"), 3, [], 50)
        self.assertEqual(deletion["panels"][1]["pos1"], 260)
        self.assertIsNone(deletion["source_end"])

    def test_complex_sequence_deletion_uses_full_reference_span(self):
        line = "\t".join(["chr1", "200", "complex", "A" * 101, "A" * 51, ".", "PASS", "."])
        r = sv.parse_record(line, 1, [], 50)
        self.assertEqual(r["panels"][1]["pos1"], 300)
        self.assertEqual(r["length"], 50)
        self.assertIsNone(r["source_end"])
        self.assertTrue(any("unnormalized" in note for note in r["notes"]))

    def test_sample_support_never_uses_caller_column_order(self):
        r = sv.parse_record(vcf_line(fmt=["GT:PR:SR", "0/1:10,8:.,4"]), 1, ["caller_label"], 50)
        self.assertEqual(sv.support_fields(r, "child"), {"GT": "NA", "PR": "NA", "SR": "NA"})
        r["samples"]["child"] = r["samples"]["caller_label"]
        self.assertEqual(sv.support_fields(r, "child"), {"GT": "0/1", "PR": "10 / 8", "SR": "NA / 4"})

    def test_joint_caller_column_name_collision_is_still_na(self):
        r = sv.parse_record(vcf_line(fmt=["GT:PR:SR", "0/1:10,8:9,4"]), 1, ["child"], 50)
        r["ignore_vcf_samples"] = True
        self.assertEqual(sv.support_fields(r, "child"), {"GT": "NA", "PR": "NA", "SR": "NA"})
        self.assertEqual(r["samples"]["child"]["PR"], "10,8")

    def test_window_clamping_and_overlap(self):
        panels = [{"chrom": "chr1", "pos1": 1}, {"chrom": "chr1", "pos1": 90}]
        self.assertEqual(sv.extraction_windows(panels, {"chr1": 100}, 60), [{"chr": "chr1", "start": 1, "end": 100}])
        self.assertEqual(sv.initial_locus(panels[0], {"start": 0, "end": 100}, 400), "chr1:1-100")

    def test_script_json_escapes_markup(self):
        self.assertNotIn("<", sv.script_json({"id": "</script><img onerror=alert(1)>"}))

    def test_pedigree_drives_role_order(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "family.ped"
            path.write_text("f mom 0 0 2 1\nf kid dad mom 1 2\nf dad 0 0 1 1\n")
            trio = sv.trio_order(path, "f", ["dad", "mom", "kid"], ["dad.bam", "mom.bam", "kid.bam"])
            self.assertEqual([r["sample"] for r in trio], ["kid", "dad", "mom"])
            with self.assertRaisesRegex(ValueError, "distinct BAM"):
                sv.trio_order(path, "f", ["dad", "mom", "kid"], ["same.bam"] * 3)

    def test_compressed_empty_vcf(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "empty.vcf.gz"
            with gzip.open(path, "wt") as stream:
                stream.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            self.assertEqual(sv.read_variants(path), ([], []))


def make_fixture(root):
    """Entirely synthetic trio; reusable by the optional offline browser check."""
    import pysam
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    reference = root / "reference.fa"
    reference.write_text(">chr1\n" + "ACGT" * 3000 + "\n>chr2\n" + "TGCA" * 3000 + "\n")
    pysam.faidx(str(reference))
    ped = root / "family.ped"
    ped.write_text("f mom 0 0 2 1\nf kid dad mom 1 2\nf dad 0 0 1 1\n")
    samples = ["mom", "dad", "kid"]
    bams = []
    for sample in samples:
        path = root / (sample + ".bam")
        bams.append(str(path))
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chr1", "LN": 12000}, {"SN": "chr2", "LN": 12000}],
                  "RG": [{"ID": "rg1", "SM": sample}]}
        with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
            for tid in (0, 1):
                for pos in range(0, 11900, 50):
                    read = pysam.AlignedSegment()
                    read.query_name = sample + "_" + str(tid) + "_" + str(pos)
                    read.query_sequence = "ACGT" * 25
                    read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                    read.flag = 0
                    read.reference_id = tid
                    read.reference_start = pos
                    read.mapping_quality = 60
                    read.cigar = [(0, 100)]
                    read.set_tag("RG", "rg1")
                    bam.write(read)
                    if tid == 0 and pos == 1800:
                        read.query_name += "_mate_off_panel"
                        read.flag = 65
                        read.next_reference_id = 1
                        read.next_reference_start = 7000
                        bam.write(read)
                    elif tid == 0 and pos == 1900:
                        read.query_name += "_spans_both"
                        read.cigar = [(0, 50), (2, 7000), (0, 50)]
                        bam.write(read)
        pysam.index(str(path))
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcaller_label\n"
    lines = [vcf_line(pos=2000, info="SVTYPE=DEL;END=7000;SVLEN=-4999", ident="del"),
             vcf_line(pos=3000, alt="<DUP>", info="SVTYPE=DUP;END=9000", ident="dup"),
             vcf_line(pos=4000, alt="<INV>", info="SVTYPE=INV;END=8000", ident="inv"),
             vcf_line(pos=5000, alt="<INS>", info="SVTYPE=INS;SVLEN=100", ident="ins"),
             vcf_line(pos=6000, alt="A]chr2:7000]", info="SVTYPE=BND;MATEID=mate", ident="bnd"),
             vcf_line(pos=1, info="SVTYPE=DEL;END=11999", ident="edge"),
             vcf_line(pos=200, alt="AT", info=".", ident="small"), "malformed"]
    vcf = root / "calls.vcf"
    vcf.write_text(header + "\n".join(line + "\tGT:PR:SR\t0/1:10,8:9,2" if line != "malformed" else line for line in lines) + "\n")
    return dict(reference=str(reference), ped=str(ped), samples=samples, bams=bams, vcf=str(vcf))


def fixture_args(fixture, output):
    args = ["--vcf", fixture["vcf"], "--ped", fixture["ped"], "--family", "f", "--reference", fixture["reference"],
            "--samples"] + fixture["samples"] + ["--bams"] + fixture["bams"] + ["--output-dir", str(output), "--flank", "500"]
    return args


@unittest.skipUnless(importlib.util.find_spec("pysam") and importlib.util.find_spec("igv_reports"), "requires pinned IGV reports environment")
class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.fixture = make_fixture(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def test_empty_vcf_creates_valid_standalone_report(self):
        Path(self.fixture["vcf"]).write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        output = self.root / "empty"
        sv.main(fixture_args(self.fixture, output))
        self.assertIn("No supported SV candidates", (output / "index.html").read_text())
        self.assertEqual(json.loads((output / "manifest.json").read_text())["rendered"], 0)

    def test_swapped_bam_identity_fails(self):
        trio = sv.trio_order(self.fixture["ped"], "f", self.fixture["samples"], list(reversed(self.fixture["bams"])))
        with self.assertRaisesRegex(ValueError, "sample mismatch"):
            sv.validate_inputs(trio, [], self.fixture["reference"])

    def test_slice_preserves_off_panel_mate_and_deduplicates_spanning_read(self):
        import pysam
        member = {"bam": self.fixture["bams"][2]}
        regions = [{"chr": "chr1", "start": 1500, "end": 2500}, {"chr": "chr1", "start": 6500, "end": 7500}]
        data = sv.slice_alignment(member, regions, SimpleNamespace(exclude_flags=1536, reference=self.fixture["reference"]))
        sliced = self.root / "sliced.bam"
        sliced.write_bytes(data)
        with pysam.AlignmentFile(str(sliced), "rb") as bam:
            reads = list(bam)
            self.assertEqual(bam.references, ("chr1", "chr2"))
            mate = next(read for read in reads if read.query_name.endswith("_mate_off_panel"))
            self.assertEqual((mate.next_reference_id, mate.next_reference_start), (1, 7000))
            self.assertEqual(sum(read.query_name.endswith("_spans_both") for read in reads), 1)

    @unittest.skipUnless(os.environ.get("IGV_REPORTS_TEST_JS"), "set IGV_REPORTS_TEST_JS to local IGV.js 3.5.2")
    def test_real_standalone_payloads_cover_all_sv_types(self):
        output = self.root / "report"
        # Snakemake can cache the script and HTML into separate directories.
        cached_helper = self.root / "source_cache" / "igv_report_sv.py"
        cached_helper.parent.mkdir()
        shutil.copyfile(SCRIPT, cached_helper)
        completed = subprocess.run([sys.executable, str(cached_helper)] + fixture_args(self.fixture, output) +
                                   ["--igv-js", os.environ["IGV_REPORTS_TEST_JS"], "--controls-template", str(SCRIPT.with_name("igv_read_controls.html"))],
                                   capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertEqual((manifest["rendered"], manifest["skipped"]), (6, 2))
        self.assertEqual([m["sample"] for m in manifest["trio"]], ["kid", "dad", "mom"])
        document = (output / "index.html").read_text()
        self.assertFalse(re.search(r'<(?:script|link)[^>]+(?:src|href)="https?://', document))
        self.assertNotIn("__READ_WINDOW_METADATA__", document)
        sessions = json.JSONDecoder().raw_decode(document.split("const sessionDictionary = ", 1)[1])[0]
        for sid, uri in sessions.items():
            session = json.loads(decode(uri))
            tracks = [track for track in session["tracks"] if track["type"] == "alignment"]
            self.assertEqual([t["name"] for t in tracks], ["Child kid", "Father dad", "Mother mom"])
            self.assertTrue(all(decode(track["url"]).startswith(b"BAM\x01") for track in tracks))
            sequence = decode(session["reference"]["fastaURL"]).decode()
            for piece in sequence.split(">")[1:]:
                h, bases = piece.split("\n", 1)
                match = re.search(r":(\d+)-(\d+) @len=(\d+)", h)
                start, end, length = map(int, match.groups())
                self.assertGreaterEqual(start, 1)
                self.assertLessEqual(end, length)
                self.assertEqual(len(bases.strip()), end - start + 1)
        self.assertEqual(manifest["records"][0]["source_svlen"], "-4999")
        self.assertTrue(all(v == {"GT": "NA", "PR": "NA", "SR": "NA"} for v in manifest["records"][0]["trio_support"].values()))


if __name__ == "__main__":
    unittest.main()
