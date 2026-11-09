"""Regress HipSTR's distinct VCF POS and repeat START coordinate contracts."""

import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/render_hipstr_viz.py"
SPEC = importlib.util.spec_from_file_location("render_hipstr_viz", SCRIPT)
viz = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(viz)


class HipstrVisualizationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.vcf = Path(self.temporary.name) / "calls.vcf.gz"

    def write_records(self, records):
        with gzip.open(self.vcf, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, pos, info in records:
                handle.write("{}\t{}\t.\tA\tAT\t.\tPASS\t{}\n".format(chrom, pos, info))

    def test_shifted_repeat_uses_explicit_start_not_vcf_pos(self):
        self.write_records([("chr1", 100, "START=104;END=150")])
        self.assertEqual(viz.resolve_locus(self.vcf, "chr1", 100), (104, 150))

    def test_unshifted_record_preserves_start(self):
        self.write_records([("chr1", 100, "START=100;END=150")])
        self.assertEqual(viz.resolve_locus(self.vcf, "chr1", 100), (100, 150))

    def test_missing_start_is_not_inferred_from_pos(self):
        self.write_records([("chr1", 100, "END=150")])
        with self.assertRaisesRegex(ValueError, "requires INFO/START"):
            viz.resolve_locus(self.vcf, "chr1", 100)

    def test_ambiguous_record_is_not_selected_arbitrarily(self):
        self.write_records([("chr1", 100, "START=104;END=150"),
                            ("chr1", 100, "START=105;END=151")])
        with self.assertRaisesRegex(ValueError, "found 2"):
            viz.resolve_locus(self.vcf, "chr1", 100)

    def test_nearby_record_is_not_used_for_missing_position(self):
        self.write_records([("chr1", 101, "START=104;END=150")])
        with self.assertRaisesRegex(ValueError, "found 0"):
            viz.resolve_locus(self.vcf, "chr1", 100)

    def test_invalid_or_duplicate_coordinates_are_rejected(self):
        for info in ("START=104;END=103", "START=0;END=150", "START=104,105;END=150",
                     "START=104;START=105;END=150", "START=104", "START=bad;END=150"):
            with self.subTest(info=info):
                self.write_records([("chr1", 100, info)])
                with self.assertRaises(ValueError):
                    viz.resolve_locus(self.vcf, "chr1", 100)

    def test_substantive_fragment_is_accepted_without_html_wrapper(self):
        fragment = ('<style>.readtable {color: black}</style>'
                    '<table class="reftable"><caption>chr1\t104\t150</caption>'
                    '<tr><td class="Areftd">A</td></tr>'
                    '<tr><td class="Atd">A</td></tr></table>')
        viz.validate_alignment_html(fragment, "chr1", 104, 150)
        with self.assertRaisesRegex(ValueError, "caption"):
            viz.validate_alignment_html(fragment, "chr1", 100, 150)

    def test_stylesheet_only_output_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "alignment content"):
            viz.validate_alignment_html('<style>.readtable {} .reftable {}</style>',
                                        "chr1", 104, 150)


if __name__ == "__main__":
    unittest.main()
