#!/usr/bin/env python3
"""Validate substantive HipSTR fragments and reject empty or mismatched reports."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/validate_full_workflow.py'
spec = importlib.util.spec_from_file_location('full_workflow_validator', SCRIPT)
validator = importlib.util.module_from_spec(spec); spec.loader.exec_module(validator)
SAMPLES = ['child', 'father', 'mother']
LOCUS = ['chr1', 104, 110]
STYLE = '<style>.reftable { color: white; } .Atd { color: purple; }</style>'
REFERENCE = '<tr><td class="Areftd">A</td></tr>'
def sample(name):
    return '<tr><td class="samplename"><font color="red">'+name+': 0|1</font></td></tr><tr><td class="Atd">A</td></tr>'
def fragment():
    return STYLE+'<div><table class="reftable"><caption>chr1 104 110</caption>'+REFERENCE+''.join(sample(n) for n in SAMPLES)+'</table><br></div>'

class VizAlnTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.path = Path(self.temporary.name)/'chr1_100.html'
    def tearDown(self): self.temporary.cleanup()
    def inspect(self, text, **kwargs):
        self.path.write_text(text)
        return validator.vizaln_html_document(self.path, SAMPLES, LOCUS, **kwargs)
    def rejects(self, text, pattern):
        with self.assertRaisesRegex(ValueError, pattern): self.inspect(text)
    def test_real_fragment_contract(self):
        d=self.inspect(fragment());self.assertTrue(d['fragment']);self.assertEqual(d['alignment_rows_by_sample'],dict.fromkeys(SAMPLES,1))
    def test_full_document_wrapper(self):
        self.assertFalse(self.inspect('<!doctype html><html><body>'+fragment()+'</body></html>')['fragment'])
    def test_stylesheet_only_is_rejected(self): self.rejects(STYLE,'alignment table')
    def test_plain_error_is_rejected(self): self.rejects('ERROR: tabix failed','alignment table')
    def test_truncated_tail_is_rejected(self): self.rejects(fragment()[:-6],'Unclosed')
    def test_mismatched_closure_is_rejected(self): self.rejects(fragment().replace('</table>','</div>'),'Unbalanced')
    def test_missing_reference_is_rejected(self): self.rejects(fragment().replace(REFERENCE,''),'reference sequence')
    def test_wrong_start_is_rejected(self): self.rejects(fragment().replace('104 110','100 110'),'caption')
    def test_wrong_chromosome_is_rejected(self): self.rejects(fragment().replace('chr1 104','chr2 104'),'caption')
    def test_wrong_end_is_rejected(self): self.rejects(fragment().replace('104 110','104 111'),'caption')
    def test_missing_sample_is_rejected(self): self.rejects(fragment().replace(sample('mother'),''),'sample sections')
    def test_duplicate_sample_is_rejected(self): self.rejects(fragment().replace('mother:','father:'),'sample sections')
    def test_wrong_sample_is_rejected(self): self.rejects(fragment().replace('mother:','other:'),'sample sections')
    def test_no_read_rows_for_one_sample_is_rejected(self):
        self.rejects(fragment().replace(sample('mother'),'<tr><td class="samplename">mother: 0|1</td></tr>'),'alignment rows')
    def test_unexpected_markup_is_rejected(self): self.rejects(fragment()+'<script>error()</script>','Unexpected VizAln tag')
    def test_bounded_html_size(self):
        with self.assertRaisesRegex(ValueError,'exceeds HTML validation bound'):self.inspect(fragment()+' '*1048576,max_html_mb=1)
    def test_generic_document_check_stays_strict(self):
        self.path.write_text(fragment())
        with self.assertRaisesRegex(ValueError,'Missing HTML header'):validator.html_document(self.path)
        self.path.write_text('<html><body>report</body></html>');self.assertTrue(validator.html_document(self.path)['html'])
    def test_start_mapping_retains_pos_filename(self):
        p=self.path.with_suffix('.vcf');p.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tT\t.\tPASS\tSTART=104;END=110\n')
        self.assertEqual(validator.hipstr_vizaln_loci(p),{'chr1_100':LOCUS})
    def test_missing_start_fails(self):
        p=self.path.with_suffix('.vcf');p.write_text('chr1\t100\t.\tA\tT\t.\tPASS\tEND=110\n')
        with self.assertRaisesRegex(ValueError,'INFO/START'):validator.hipstr_vizaln_loci(p)
    def test_ambiguous_start_fails(self):
        p=self.path.with_suffix('.vcf');p.write_text('chr1\t100\t.\tA\tT\t.\tPASS\tSTART=104;START=105;END=110\n')
        with self.assertRaisesRegex(ValueError,'unique integer'):validator.hipstr_vizaln_loci(p)

if __name__ == '__main__':
    unittest.main()
