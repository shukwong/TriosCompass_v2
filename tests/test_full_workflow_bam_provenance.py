"""End-to-end auditor regressions using small, indexed synthetic trio BAMs.

These tests exercise the real CLI, samtools retrieval and IGV payload checks.
They require the existing CI dependencies PyYAML, pysam and samtools.
"""
import base64
import contextlib
import gzip
import importlib.util
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest import mock

try:
    import pysam
    import yaml
except ImportError:
    pysam = yaml = None

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/validate_full_workflow.py'
spec = importlib.util.spec_from_file_location('bam_provenance_validator', SCRIPT)
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)
SAMTOOLS = shutil.which('samtools')
TRIO = ['Z_child', 'A_father', 'B_mother']
VCF_HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'


def data_uri(content):
    return 'data:application/octet-stream;base64,' + base64.b64encode(content).decode()


@unittest.skipUnless(pysam and yaml and SAMTOOLS, 'requires pysam, PyYAML and samtools')
class BamProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='trios-audit-provenance-')
        self.root = Path(self.temporary.name)
        self.out = self.root / 'custom-output'
        self.raw = {sample: self.root / 'input' / (sample + '.bam') for sample in TRIO}
        self.fixed = {sample: self.out / 'fixed-rg' / (sample + '.bam') for sample in TRIO}
        self.reference = self.root / 'reference.fa'
        self.reference.write_text('>chr1\n' + 'A' * 100 + '\n')
        Path(str(self.reference) + '.fai').write_text('chr1\t100\t6\t100\t101\n')
        self.report_dir = self.out / 'igv_reports' / 'F'
        self.manifest_path = self.report_dir / 'manifest.json'
        self.rufus_manifest = self.out / 'rufus' / 'work' / 'F' / 'attempt-test' / 'run.json'

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, path, content):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)

    def write_bam(self, path, sample, *, length=100, count=1):
        path.parent.mkdir(parents=True, exist_ok=True)
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': 'chr1', 'LN': length}],
                  'RG': [{'ID': sample, 'SM': sample}]}
        with pysam.AlignmentFile(str(path), 'wb', header=header) as handle:
            for index in range(count):
                read = pysam.AlignedSegment(handle.header)
                read.query_name = '{}-read{}'.format(path.stem, index)
                read.query_sequence = 'AAAAA'
                read.flag = 0
                read.reference_id = 0
                read.reference_start = index
                read.mapping_quality = 60
                read.cigarstring = '5M'
                read.query_qualities = pysam.qualitystring_to_array('IIIII')
                read.set_tag('RG', sample)
                handle.write(read)
        pysam.index(str(path))

    def write_vcf(self, relative, samples=()):
        path = self.out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        header = VCF_HEADER + ('\tFORMAT\t' + '\t'.join(samples) if samples else '') + '\n'
        if path.suffix == '.gz':
            with gzip.open(path, 'wt') as handle:
                handle.write(header)
            Path(str(path) + '.tbi').write_bytes(b'TBI\x01synthetic-empty-index')
        else:
            path.write_text(header)
        return path

    def make_run(self, reset=True):
        self.config = {'output_dir': 'custom-output', 'ped_dir': 'input/ped',
                       'pepfile': 'input/pep.yaml', 'ref': {'sequence': 'reference.fa'},
                       'bam_input': {'reset_RG': reset},
                       'igv_reports': {'enable': True, 'source': 'joint'},
                       'rufus': {'enable': True},
                       'bam_qc': {'flagstat': {'enable': True}}}
        self.config_path = self.root / 'config/config.yaml'
        self.write(self.config_path, json.dumps(self.config))
        self.write(self.root / 'input/pep.yaml', json.dumps({'sample_table': 'samples.csv'}))
        self.write(self.root / 'input/samples.csv', 'SAMPLE_ID,BAM\n' + ''.join(
            '{},input/{}.bam\n'.format(sample, sample) for sample in reversed(TRIO)))
        self.write(self.root / 'input/ped/F.ped',
                   'F Z_child A_father B_mother 1 2\nF A_father 0 0 1 1\nF B_mother 0 0 2 1\n')
        for sample in TRIO:
            self.write_bam(self.raw[sample], 'old_' + sample if reset else sample)
            if reset:
                self.write_bam(self.fixed[sample], sample)
            for caller, suffix in [('gatkhc_pb', '.gatk.g.vcf.gz'),
                                   ('deepvariant_pb', '.deepvariant.g.vcf.gz')]:
                self.write_vcf(caller + '/' + sample + suffix, [sample])
            self.write(self.out / 'flagstat' / (sample + '.flagstat'), '1 + 0 in total (QC-passed reads + QC-failed reads)\n')
            self.write(self.out / 'flagstat' / (sample + '.idxstat'), 'chr1\t100\t1\t0\n*\t0\t0\t0\n')
        self.write_vcf('dnm_vcf/F.dnm.vcf.gz', [TRIO[0]])
        self.write(self.out / 'dnm_vcf_summary/DNM_summary.txt', 'TrioID\tDNM_Cnt\nF\t0\n')
        for relative in ['glnexus/F.dv_combined.vcf.gz', 'gatk_cgp/F.cgp_norm.vcf.gz']:
            self.write_vcf(relative, TRIO)
        rufus_vcf = self.write_vcf('rufus/F.rufus.vcf', TRIO)
        self.write(self.out / 'rufus/rufus_summary.tsv', 'Family\tTotalRecords\nF\t0\n')
        self.effective = self.fixed if reset else self.raw
        self.write(self.rufus_manifest, json.dumps({
            'family': 'F', 'output': str(rufus_vcf), 'reference': str(self.reference),
            'trio': [{'sample': sample, 'bam': str(self.effective[sample])} for sample in TRIO],
            'status': 'complete', 'finished_utc': '2000-01-01T00:00:00Z',
            'counts': {'total_records': 0}}))
        self.write_report()

    def write_report(self, *, paths=None, payloads=None):
        paths = paths or self.effective
        payloads = payloads or self.effective
        source_vcf = self.write_vcf('joint_dnSV/F.dnSV.vcf')
        source_record = 'chr1\t1\tv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5'
        with source_vcf.open('a') as handle:
            handle.write(source_record + '\n')
        record = {'source_record': source_record, 'status': 'rendered', 'svtype': 'DEL',
                  'session_id': '0', 'initial_locus': 'chr1:1-5', 'ignore_vcf_samples': True,
                  'extraction_windows_1based_inclusive': [{'chr': 'chr1', 'start': 1, 'end': 5}],
                  'trio_support': {sample: {'GT': 'NA', 'PR': 'NA', 'SR': 'NA'} for sample in TRIO}}
        manifest = {'family': 'F', 'trio': [
            {'sample': sample, 'bam': str(paths[sample]), 'label': sample} for sample in TRIO],
            'reference': str(self.reference), 'source_vcf': str(source_vcf),
            'source_vcf_sha256': validator.sha256(source_vcf), 'records': [record],
            'rendered': 1, 'skipped': 0, 'settings': {'exclude_flags': 0}}
        session = {'locus': 'chr1:1-5',
                   'reference': {'fastaURL': data_uri(b'>chr1:1-5 @len=100\nAAAAA\n')},
                   'tracks': [{'type': 'alignment', 'name': sample,
                               'url': data_uri(payloads[sample].read_bytes())} for sample in TRIO]}
        document = '<html><body>read-zoom-in read-row-height<script>const sessionDictionary = ' + json.dumps(
            {'0': data_uri(json.dumps(session).encode())}) + ';</script></body></html>'
        self.write(self.manifest_path, json.dumps(manifest))
        for name in ['index.html', 'variants.html']:
            self.write(self.report_dir / name, document)
        self.write(self.report_dir / 'variants.tsv', 'record\nv1\n')
        self.write(self.report_dir / 'skipped_records.tsv', 'record\n')

    def run_audit(self, *extra):
        command = [str(SCRIPT), '--run-dir', str(self.root), '--family', 'F',
                   '--allow-disabled', '--skip-extra-igv', '--samtools', SAMTOOLS, *extra]
        with mock.patch.object(sys, 'argv', command), contextlib.redirect_stdout(io.StringIO()):
            code = validator.main()
        report = self.root / 'validation/validation.json'
        return code, json.loads(report.read_text()) if report.exists() else None

    def assert_failure(self, result, branch, message):
        code, report = result
        self.assertEqual(code, 1)
        matches = [row for row in report['checks'] if row['branch'] == branch and
                   row['status'] == 'FAIL' and message in row.get('error', '')]
        self.assertTrue(matches, report['checks'])

    def test_reset_enabled_uses_corrected_bams_and_checks_real_payloads(self):
        self.make_run()
        code, report = self.run_audit()
        self.assertEqual(code, 0, report['checks'])
        headers = [row for row in report['checks'] if row['check'] == 'bam_header_structure']
        self.assertEqual(len(headers), 3)
        self.assertTrue(all(not row['details']['sample_identity_checked'] for row in headers))
        payload = next(row for row in report['checks'] if row['branch'] == 'igv_joint' and
                       row['check'] == 'source_provenance_and_bounded_payload_validation')
        self.assertEqual(payload['details']['bam_payloads_checked'], 3)
        snapshot = next(row['details'] for row in report['checks'] if row['check'] == 'source_metadata_snapshot')
        source_paths = {row['path'] for row in snapshot['files']}
        self.assertTrue({str(path.resolve()) for path in self.raw.values()}.issubset(source_paths))
        self.assertFalse({str(path.resolve()) for path in self.fixed.values()} & source_paths)

    def test_reset_disabled_uses_original_inputs(self):
        self.make_run(reset=False)
        code, report = self.run_audit()
        self.assertEqual(code, 0, report['checks'])
        self.assertFalse(any(row['branch'] == 'workflow_bams' for row in report['checks']))
        self.assertEqual(sum(row['check'] == 'bam_trio_identity' for row in report['checks']), 3)

    def test_reset_disabled_still_requires_original_sample_identity(self):
        self.make_run(reset=False)
        self.write_bam(self.raw[TRIO[0]], 'wrong')
        self.assert_failure(self.run_audit(), 'source_identity', 'BAM RG/SM')

    def test_missing_corrected_bam_fails(self):
        self.make_run()
        self.fixed[TRIO[0]].unlink()
        self.assert_failure(self.run_audit(), 'workflow_bams', 'Missing required file')

    def test_wrong_corrected_sample_fails(self):
        self.make_run()
        self.write_bam(self.fixed[TRIO[0]], 'wrong')
        self.assert_failure(self.run_audit(), 'workflow_bams', 'BAM RG/SM')

    def test_missing_corrected_index_fails(self):
        self.make_run()
        Path(str(self.fixed[TRIO[0]]) + '.bai').unlink()
        self.assert_failure(self.run_audit(), 'workflow_bams', 'Missing required file')

    def test_corrected_dictionary_must_match_original(self):
        self.make_run()
        self.write_bam(self.fixed[TRIO[0]], TRIO[0], length=101)
        self.assert_failure(self.run_audit(), 'workflow_bams', 'changed the sequence dictionary')

    def test_corrected_counts_must_match_original(self):
        self.make_run()
        self.write_bam(self.fixed[TRIO[0]], TRIO[0], count=2)
        self.assert_failure(self.run_audit(), 'workflow_bams', 'changed alignment counts')

    def test_igv_original_path_is_rejected_after_reset(self):
        self.make_run()
        self.write_report(paths=self.raw)
        self.assert_failure(self.run_audit(), 'igv_joint', 'IGV uses a different BAM')

    def test_igv_corrected_path_is_rejected_without_reset(self):
        self.make_run(reset=False)
        self.write_report(paths=self.fixed)
        self.assert_failure(self.run_audit(), 'igv_joint', 'IGV uses a different BAM')

    def test_igv_payload_dictionary_must_match_effective_source(self):
        self.make_run()
        other = self.root / 'other.bam'
        self.write_bam(other, TRIO[0], length=101)
        self.write_report(payloads={**self.fixed, TRIO[0]: other})
        self.assert_failure(self.run_audit(), 'igv_joint', 'sequence dictionary differs from source')

    def test_igv_payload_records_must_match_effective_source(self):
        self.make_run()
        other = self.root / 'other.bam'
        self.write_bam(other, TRIO[0], count=2)
        self.write_report(payloads={**self.fixed, TRIO[0]: other})
        self.assert_failure(self.run_audit(), 'igv_joint', 'records differ from effective workflow source')

    def test_rufus_original_path_is_rejected_after_reset(self):
        self.make_run()
        manifest = json.loads(self.rufus_manifest.read_text())
        manifest['trio'][0]['bam'] = str(self.raw[TRIO[0]])
        self.rufus_manifest.write_text(json.dumps(manifest))
        self.assert_failure(self.run_audit(), 'RUFUS', 'RUFUS uses a different BAM')

    def test_rufus_wrong_reference_is_rejected(self):
        self.make_run()
        manifest = json.loads(self.rufus_manifest.read_text())
        manifest['reference'] = str(self.root / 'wrong.fa')
        self.rufus_manifest.write_text(json.dumps(manifest))
        self.assert_failure(self.run_audit(), 'RUFUS', 'RUFUS reference mismatch')

    def test_source_snapshot_can_precede_read_group_replacement(self):
        self.make_run()
        for bam in self.fixed.values():
            bam.unlink()
            Path(str(bam) + '.bai').unlink()
        code, report = self.run_audit('--capture-sources', 'before.json')
        self.assertEqual(code, 0)
        self.assertIsNone(report)
        for sample, bam in self.fixed.items():
            self.write_bam(bam, sample)
        code, report = self.run_audit('--source-snapshot', 'before.json')
        self.assertEqual(code, 0, report['checks'])

    def test_raw_source_change_is_still_rejected_after_reset(self):
        self.make_run()
        code, _ = self.run_audit('--capture-sources', 'before.json')
        self.assertEqual(code, 0)
        self.write_bam(self.raw[TRIO[0]], 'changed_legacy_name')
        self.assert_failure(self.run_audit('--source-snapshot', 'before.json'),
                            'provenance', 'Source size,mtime or small-file hash changed')

    def test_raw_bam_integrity_is_still_checked_after_reset(self):
        self.make_run()
        self.raw[TRIO[0]].write_bytes(b'corrupt BAM')
        self.assert_failure(self.run_audit(), 'source_identity', 'Command failed')

    def test_raw_bam_index_is_still_required_after_reset(self):
        self.make_run()
        Path(str(self.raw[TRIO[0]]) + '.bai').unlink()
        self.assert_failure(self.run_audit(), 'source_identity', 'Missing required file')


if __name__ == '__main__':
    unittest.main()
