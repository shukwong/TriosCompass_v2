#!/usr/bin/env python3
"""Post-run TriosCompass artifact/provenance audit; never a GIAB truth benchmark.

Requires Python >=3.9, PyYAML and samtools on PATH; no pysam/browser required.
Default: --run-dir RUN reads RUN/config/config.yaml and writes RUN/validation/.
Missing/partial outputs are recorded and return exit1. Raw RUFUS defects are a
separate integrity warning even when all required final artifacts validate.
Large intermediate VCFs and BAMs are checked with explicitly bounded samples.
"""
import argparse, base64, collections, csv, gzip, hashlib, html, io, json, math
import os, re, shutil, socket, struct, subprocess, sys, tempfile
from datetime import datetime, timezone
from pathlib import Path
from html.parser import HTMLParser

MISSING = {'', '.', '?', 'NA', 'N/A', 'null', 'None'}
VCF_FIXED = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def merge_dict(left, right):
    for key, value in right.items():
        if isinstance(value, dict) and isinstance(left.get(key), dict):
            merge_dict(left[key], value)
        else:
            left[key] = value
    return left


def require(condition, message):
    if not condition:
        raise ValueError(message)


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()


def noncomment_table(path, delimiter='\t', allow_empty=False):
    with Path(path).open(newline='') as handle:
        lines = [line for line in handle if line.strip() and not line.startswith('#')]
    require(lines or allow_empty, 'Missing table header')
    if not lines:
        return []
    rows = list(csv.DictReader(lines, delimiter=delimiter))
    require(all(None not in row and all(v is not None for v in row.values()) for row in rows), 'Ragged table')
    return rows


def finite(value, name, minimum=None):
    number = float(value)
    require(math.isfinite(number), name + ' is nonfinite')
    require(minimum is None or number >= minimum, name + ' is below minimum')
    return number


def parse_vcf(path, sample_set=None, sample_order=None, limit=0, keep_keys=False):
    """Stream structural validation; limit0 scans every record including gzip CRC."""
    header = None
    count = 0
    types = collections.Counter()
    keys = []
    malformed = []
    completed = True
    with open_text(path) as handle:
        first = handle.readline()
        require(re.fullmatch(r'##fileformat=VCFv4\.\d+\r?\n?', first), 'Invalid/missing VCF fileformat')
        for line_number, line in enumerate(handle, 2):
            if line.startswith('##'):
                require(header is None, 'Metadata after column header')
                continue
            if line.startswith('#CHROM\t'):
                require(header is None, 'Duplicate column header')
                header = line.rstrip('\r\n').split('\t')
                require(header[:8] == VCF_FIXED and len(header) != 9, 'Invalid VCF columns')
                require(len(header) == 8 or header[8] == 'FORMAT', 'Missing FORMAT column')
                require(len(header[9:]) == len(set(header[9:])), 'Duplicate sample IDs')
                continue
            require(header is not None, 'Record before VCF column header')
            fields = line.rstrip('\r\n').split('\t')
            reason = None
            if len(fields) != len(header):
                reason = 'column count {} != {}'.format(len(fields), len(header))
            elif any(not field or any(c.isspace() for c in field) for field in fields):
                reason = 'empty or whitespace-containing field'
            elif not fields[1].isdigit() or int(fields[1]) < 1 or not re.fullmatch('[ACGTNacgtn]+', fields[3]):
                reason = 'invalid POS/REF'
            elif fields[5] != '.':
                try:
                    finite(fields[5], 'QUAL', 0)
                except ValueError:
                    reason = 'invalid QUAL'
            if reason:
                if len(malformed) < 20:
                    malformed.append({'line': line_number, 'reason': reason, 'chrom': fields[0], 'pos': fields[1] if len(fields) > 1 else None})
            elif len(fields) >= 8:
                for item in fields[7].split(';'):
                    if item.startswith('SVTYPE='):
                        types[item.split('=', 1)[1]] += 1
                if keep_keys:
                    keys.append('{}:{}:{}:{}'.format(fields[0], fields[1], fields[3], fields[4]))
            count += 1
            if limit and count >= limit:
                completed = False
                break
    require(header is not None, 'Missing VCF column header')
    if sample_set is not None:
        require(set(header[9:]) == set(sample_set), 'VCF samples {} != {}'.format(header[9:], sample_set))
    if sample_order is not None:
        require(header[9:] == list(sample_order), 'VCF order {} != required {}'.format(header[9:], sample_order))
    require(not malformed, 'Malformed VCF rows: ' + json.dumps(malformed))
    result = {'records_examined': count, 'full_scan': completed, 'samples': header[9:], 'svtype_counts': dict(types)}
    if keep_keys:
        result['variant_keys'] = keys
    return result


def inspect_index(path):
    path = Path(path)
    require(path.is_file() and path.stat().st_size > 0, 'Missing/empty index')
    with path.open('rb') as handle:
        prefix = handle.read(4)
    if prefix[:2] == b'\x1f\x8b':
        with gzip.open(path, 'rb') as handle:
            prefix = handle.read(4)
    require(prefix in (b'TBI\x01', b'CSI\x01', b'BAI\x01'), 'Unrecognized index signature')
    return {'magic': prefix.hex(), 'bytes': path.stat().st_size}


def read_ped(path, family):
    members = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        row = line.split()
        require(len(row) >= 6, 'Malformed PED row')
        if row[0] == family:
            require(row[1] not in members, 'Duplicate PED sample')
            members[row[1]] = row
    children = [row for row in members.values() if row[2] != '0' and row[3] != '0']
    require(len(children) == 1, 'Expected one complete trio child')
    child = children[0]
    order = [child[1], child[2], child[3]]
    require(len(set(order)) == 3 and all(s in members for s in order), 'Missing/distinct trio members')
    return order


class Audit:
    def __init__(self, args):
        self.args = args
        self.run = Path(args.run_dir).absolute()
        self.report_dir = Path(args.report_dir).absolute() if args.report_dir else self.run / 'validation'
        self.checks = []
        self.files = {}
        self.samtools = args.samtools or shutil.which('samtools')
        self.config = {}
        self.bam_headers = {}

    def path(self, value):
        value = Path(value).expanduser()
        return value if value.is_absolute() else self.run / value

    def check(self, branch, name, function, path=None, required=True, kind='validation'):
        row = {'branch': branch, 'check': name, 'required': required, 'kind': kind}
        if path is not None:
            row['path'] = str(path)
        try:
            row['details'] = function()
            row['status'] = 'PASS'
        except Exception as error:
            row['status'] = 'FAIL' if required else 'WARN'
            row['error'] = '{}: {}'.format(type(error).__name__, error)
        self.checks.append(row)
        return row.get('details') if row['status'] == 'PASS' else None

    def file(self, branch, path, empty=False, required=True):
        path = Path(path)
        def inspect():
            require(path.is_file(), 'Missing required file')
            stat = path.stat()
            require(empty or stat.st_size > 0, 'Empty required file')
            info = {'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
            if stat.st_size <= self.args.hash_limit_mb * 1024 * 1024:
                info['sha256'] = sha256(path)
            self.files[str(path)] = info
            return info
        return self.check(branch, 'artifact_exists', inspect, path, required, 'artifact') is not None

    def command(self, argv, binary=False, timeout=180):
        result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
        require(result.returncode == 0, 'Command failed {}: {}'.format(argv, result.stderr.decode(errors='replace')[-3000:]))
        return result.stdout if binary else result.stdout.decode()

    def sam(self, *arguments, binary=False):
        require(bool(self.samtools), 'samtools not available on PATH; pass --samtools')
        return self.command([self.samtools, *map(str, arguments)], binary=binary)

    def vcf(self, branch, path, samples=None, order=None, full=True, keys=False, indexed=False):
        if not self.file(branch, path):
            return None
        parsed = self.check(branch, 'vcf_structure_and_samples', lambda: parse_vcf(path, samples, order, 0 if full else self.args.intermediate_records, keys), path)
        if indexed:
            candidates = [Path(str(path) + suffix) for suffix in ('.tbi', '.csi')]
            index = next((p for p in candidates if p.exists()), candidates[0])
            if self.file(branch, index):
                self.check(branch, 'index_signature', lambda: inspect_index(index), index)
        return parsed

    def bam_header(self, path, sample):
        header = self.sam('view', '-H', path)
        samples = {field[3:] for line in header.splitlines() if line.startswith('@RG\t') for field in line.split('\t')[1:] if field.startswith('SM:')}
        require(samples == {sample}, 'BAM RG/SM {} != {}'.format(sorted(samples), sample))
        reference = [line for line in header.splitlines() if line.startswith('@SQ\t')]
        require(reference, 'BAM has no sequence dictionary')
        self.bam_headers[str(path)] = {'sample': sample, 'reference': reference}
        return {'samples': sorted(samples), 'contigs': len(reference), 'header_sha256': hashlib.sha256(header.encode()).hexdigest()}

    def bam(self, branch, path, sample, filtering=None):
        if not self.file(branch, path):
            return None
        self.check(branch, 'bam_quickcheck', lambda: {'output': self.sam('quickcheck', '-v', path)}, path)
        self.check(branch, 'bam_trio_identity', lambda: self.bam_header(path, sample), path)
        index = next((p for p in (Path(str(path)+'.bai'), Path(path).with_suffix('.bai'), Path(str(path)+'.csi')) if p.is_file()), Path(str(path)+'.bai'))
        if self.file(branch, index):
            self.check(branch, 'bam_index_signature', lambda: inspect_index(index), index)
        def indexed_counts():
            require(index.is_file(), 'No BAM index; refusing a whole-file idxstats fallback')
            rows = [line.split('\t') for line in self.sam('idxstats', path).splitlines()]
            require(rows and all(len(r) == 4 for r in rows), 'Malformed idxstats')
            total = sum(int(r[2])+int(r[3]) for r in rows)
            require(all(int(r[2]) >= 0 and int(r[3]) >= 0 for r in rows), 'Negative BAM counts')
            mapped = next((r for r in rows if r[0] != '*' and int(r[2]) > 0), None)
            examined = 0
            if mapped:
                region = '{}:1-{}'.format(mapped[0], min(int(mapped[1]), 1000000))
                require(bool(self.samtools), 'samtools missing')
                process = subprocess.Popen([self.samtools, 'view', str(path), region], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
                limited = False
                try:
                    for line in process.stdout:
                        f = line.rstrip().split('\t')
                        require(len(f) >= 11, 'Malformed SAM record in indexed retrieval')
                        flag, mapq, cigar = int(f[1]), int(f[4]), f[5]
                        if filtering:
                            require(not flag & 0x400, 'Duplicate read in redup export')
                        if filtering == 'hmq':
                            require(flag & 2 and mapq >= 60 and not flag & 0x500 and 'M' in cigar and not re.search('[SHPN]', cigar), 'Read violates HMQ export filter')
                        examined += 1
                        if examined >= self.args.bam_records:
                            limited = True
                            break
                finally:
                    process.stdout.close()
                    if process.poll() is None and (limited or sys.exc_info()[0] is not None):
                        process.terminate()
                    code = process.wait(timeout=30)
                require(limited or code == 0, 'Indexed BAM retrieval failed')
            return {'indexed_total_records': total, 'sampled_alignment_records': examined, 'sampling_scope': 'first mapped contig, first1Mb, at most{} records'.format(self.args.bam_records)}
        return self.check(branch, 'bam_index_counts_and_bounded_filter_check', indexed_counts, path)

    def source_snapshot(self, paths):
        snapshot = {'schema_version': 1, 'files': []}
        for path in paths:
            path = Path(path).resolve()
            require(path.is_file(), 'Missing source: ' + str(path))
            st = path.stat()
            item = {'path': str(path), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}
            if st.st_size <= self.args.hash_limit_mb * 1024 * 1024:
                item['sha256'] = sha256(path)
            snapshot['files'].append(item)
        return snapshot

    def report(self, source, directory, trio, bams, reference):
        branch = 'igv_' + source
        manifest_file = directory / 'manifest.json'
        available = all([self.file(branch, directory/name) for name in ('index.html', 'variants.html', 'manifest.json', 'variants.tsv', 'skipped_records.tsv')])
        if not available:
            return
        def validate():
            manifest = json.loads(manifest_file.read_text())
            require(manifest['family'] == self.args.family, 'IGV family mismatch')
            require([m['sample'] for m in manifest['trio']] == trio, 'IGV trio order mismatch')
            for m in manifest['trio']:
                require(Path(m['bam']).resolve() == bams[m['sample']].resolve(), 'IGV uses a different BAM')
            require(Path(manifest['reference']).resolve() == reference.resolve(), 'IGV reference mismatch')
            expected_vcf = self.output / {'joint':'joint_dnSV/{f}.dnSV.vcf','manta':'dnSV_manta/{f}.dnSV.manta.vcf','rufus':'rufus/{f}.rufus.vcf'}[source].format(f=self.args.family)
            require(Path(manifest['source_vcf']).resolve() == expected_vcf.resolve(), 'IGV source is not this full run')
            require(sha256(expected_vcf) == manifest['source_vcf_sha256'], 'IGV source hash changed')
            require(Path(str(reference)+'.fai').is_file(), 'Reference index missing; refusing to create shared-source index')
            require(sha256(directory/'index.html') == sha256(directory/'variants.html'), 'Standalone report entries differ')
            records = manifest['records']
            with open_text(expected_vcf) as source_handle:
                source_lines = source_handle.read().splitlines()
            original_rows = [line for line in source_lines if line and not line.startswith('#')]
            require(original_rows == [r['source_record'] for r in records], 'Report omitted/reordered/changed source VCF records')
            source_samples = next(line for line in source_lines if line.startswith('#CHROM')).split('\t')[9:]
            for record, raw in zip(records, original_rows):
                fields = raw.split('\t')
                source_fields = {sample: dict(zip(fields[8].split(':'), val.split(':'))) for sample, val in zip(source_samples, fields[9:])}
                if record['status'] == 'rendered':
                    for sample in trio:
                        expected = source_support({} if source == 'joint' else source_fields.get(sample, {}))
                        require(record['trio_support'][sample] == expected, 'Report GT/PR/SR differs from source VCF')
            rendered = [r for r in records if r['status'] == 'rendered']
            require(len(rendered) == manifest['rendered'] and len(records)-len(rendered) == manifest['skipped'], 'IGV manifest counts disagree')
            require(len(noncomment_table(directory/'variants.tsv')) == len(records), 'IGV manifest/TSV row count differs')
            if source == 'joint':
                require(all(r.get('ignore_vcf_samples') for r in records), 'Native joint caller columns incorrectly treated as trio genotypes')
                require(all(all(all(v == 'NA' for v in fields.values()) for fields in r.get('trio_support',{}).values()) for r in rendered), 'Joint report inferred biological trio support')
            require((directory/'index.html').stat().st_size <= self.args.max_html_mb*1024*1024, 'HTML exceeds bounded validator size; rerun with explicitly larger --max-html-mb')
            document = (directory/'index.html').read_text()
            require('<html' in document.lower(), 'Missing HTML document')
            require(not re.search(r'<script\b[^>]*\bsrc=["\x27]https?://', document, re.I), 'Report has a remote script dependency')
            if not rendered:
                return {'rendered':0,'skipped':manifest['skipped'],'payload_sessions_checked':0,'browser_tested':False}
            require('read-zoom-in' in document and 'read-row-height' in document, 'Missing read controls')
            sessions = json.JSONDecoder().raw_decode(document.split('const sessionDictionary = ',1)[1])[0]
            require(len(sessions) == len(rendered), 'Embedded session count differs')
            selected=[]
            for row in rendered:
                if row['svtype'] not in {r['svtype'] for r in selected}:
                    selected.append(row)
            for row in rendered:
                if len(selected) >= self.args.report_sessions:
                    break
                if row not in selected:
                    selected.append(row)
            selected=selected[:self.args.report_sessions]
            payloads=0;references=0
            def decode(uri, raw=False):
                value=base64.b64decode(uri.split(',',1)[1])
                return value if raw or not value.startswith(b'\x1f\x8b') else gzip.decompress(value)
            with tempfile.TemporaryDirectory(prefix='payload-',dir=self.report_dir) as temp:
                for record in selected:
                    session=json.loads(decode(sessions[record['session_id']]))
                    require(session['locus'] == record['initial_locus'], 'IGV locus differs')
                    for chunk in decode(session['reference']['fastaURL']).decode().split('>')[1:]:
                        lines=chunk.splitlines();match=re.fullmatch(r'([^:]+):(\d+)-(\d+) @len=(\d+)',lines[0])
                        require(bool(match), 'Invalid embedded FASTA coordinate')
                        fetched=self.sam('faidx',reference,'{}:{}-{}'.format(match[1],match[2],match[3]))
                        require(''.join(lines[1:]) == ''.join(fetched.splitlines()[1:]), 'Embedded reference differs')
                        references+=1
                    tracks=[t for t in session['tracks'] if t.get('type')=='alignment']
                    require([t['name'] for t in tracks] == [m['label'] for m in manifest['trio']], 'Embedded track order differs')
                    regions=['{}:{}-{}'.format(r['chr'],r['start'],r['end']) for r in record['extraction_windows_1based_inclusive']]
                    for index,(track,member) in enumerate(zip(tracks,manifest['trio'])):
                        path=Path(temp)/('payload{}.bam'.format(index));path.write_bytes(decode(track['url'],True))
                        self.bam_header(path,member['sample'])
                        require(self.bam_headers[str(path)]['reference'] == self.bam_headers[str(bams[member['sample']])]['reference'], 'Embedded BAM sequence dictionary differs from source')
                        actual=collections.Counter(self.sam('view',path).splitlines())
                        expected=collections.Counter(self.sam('view','-M','-F',str(manifest['settings']['exclude_flags']),member['bam'],*regions).splitlines())
                        require(actual == expected, 'Embedded BAM records differ from original indexed source')
                        payloads+=1
            return {'rendered':len(rendered),'skipped':manifest['skipped'],'payload_sessions_checked':len(selected),'all_sessions_payload_checked':len(selected)==len(rendered),'bam_payloads_checked':payloads,'reference_segments_checked':references,'browser_tested':False}
        self.check(branch, 'source_provenance_and_bounded_payload_validation', validate, manifest_file)

    def finish(self):
        self.report_dir.mkdir(parents=True,exist_ok=True)
        failures=[row for row in self.checks if row['status']=='FAIL']
        warnings=[row for row in self.checks if row['status']=='WARN']
        missing=[row for row in failures if row['kind']=='artifact']
        # Detailed variant keys aid internal crosschecks but need not inflate the report.
        for row in self.checks:
            if isinstance(row.get('details'),dict):
                row['details'].pop('variant_keys',None)
        result={'schema_version':1,'created_utc':datetime.now(timezone.utc).isoformat(),'host':socket.gethostname(),'allocation_id':os.environ.get('SLURM_JOB_ID'),'command':[sys.executable]+sys.argv,'run_dir':str(self.run),'output_dir':str(getattr(self,'output',self.run/'output')),'config_files':[str(p) for p in getattr(self,'config_paths',[])],'validation_passed':not failures,'workflow_artifacts_complete':not missing and any(r['kind']=='artifact' and r['required'] for r in self.checks) and not any(r['branch']=='configuration' for r in failures),'scheduler_completion':'Not inferred from artifact presence; reconcile controller/scheduler receipt separately.','accuracy_benchmark':{'status':'NOT_ASSESSED','reason':'No truth-set comparison, precision/recall, sensitivity, or de novo biological validation is performed.'},'raw_rufus_integrity':'WARN' if any(r['branch']=='rufus_raw' and r['status']=='WARN' for r in self.checks) else 'PASS' if any(r['branch']=='rufus_raw' and r['status']=='PASS' for r in self.checks) else 'NOT_ASSESSED','limits':{'intermediate_vcf_records_per_file':self.args.intermediate_records,'bam_records_per_file':self.args.bam_records,'report_sessions_per_report':self.args.report_sessions,'browser_tested':False,'large_file_content_hashing':False,'bam_filter_validation':'Bounded indexed sample, not a full-genome read-by-read comparison.'},'failure_count':len(failures),'warning_count':len(warnings),'checks':self.checks,'file_provenance':self.files}
        (self.report_dir/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
        with (self.report_dir/'checks.tsv').open('w',newline='') as f:
            writer=csv.writer(f,delimiter='\t');writer.writerow(['branch','check','status','required','path','details'])
            for row in self.checks:
                writer.writerow([row['branch'],row['check'],row['status'],row['required'],row.get('path',''),row.get('error',json.dumps(row.get('details',{}),sort_keys=True))])
        lines=['Validation: '+('PASS' if not failures else 'FAIL'),'Required artifact presence: '+('complete' if result['workflow_artifacts_complete'] else 'incomplete'),'Failures: {}; warnings: {}'.format(len(failures),len(warnings)),'RUFUS raw integrity: '+result['raw_rufus_integrity'],'Accuracy benchmark: not assessed; no truth comparison.','Browser testing: not performed.']
        for row in failures+warnings:
            lines.append('{} {} / {}: {}'.format(row['status'],row['branch'],row['check'],row.get('error','')))
        (self.report_dir/'summary.txt').write_text('\n'.join(lines)+'\n')
        print('\n'.join(lines[:6]));print('Evidence: '+str(self.report_dir/'validation.json'))
        return 1 if failures else 0


def run(audit):
    args=audit.args
    audit.report_dir.mkdir(parents=True,exist_ok=True)
    audit.config_paths=[audit.path(p) for p in (args.config or ['config/config.yaml'])]
    def load_config():
        import yaml
        config={}
        for path in audit.config_paths:
            require(path.is_file(),'Missing config '+str(path))
            value=yaml.safe_load(path.read_text());require(isinstance(value,dict),'Config is not a mapping')
            merge_dict(config,value)
        return config
    config=audit.check('configuration','load_configuration',load_config)
    if config is None:
        return audit.finish()
    audit.config=config
    audit.output=audit.path(args.output_dir or config.get('output_dir','output'))
    out=audit.output;family=args.family
    ped=audit.path(args.ped) if args.ped else audit.path(config.get('ped_dir','input/ped'))/(family+'.ped')
    trio=audit.check('source_identity','pedigree_trio',lambda:read_ped(ped,family),ped)
    if trio is None:
        return audit.finish()
    child,father,mother=trio;family_order=[father,mother,child]
    def input_bams():
        import yaml
        pep=audit.path(config['pepfile']);metadata=yaml.safe_load(pep.read_text())
        table=Path(metadata['sample_table']);table=table if table.is_absolute() else pep.parent/table
        with table.open() as f: rows=list(csv.DictReader(f))
        mapping={row['SAMPLE_ID']:audit.path(row['BAM']) for row in rows}
        require(len(mapping)==len(rows),'Duplicate input manifest samples')
        require(set(trio).issubset(mapping),'PED samples missing from PEP')
        require(len({mapping[s].resolve() for s in trio})==3,'Trio BAM paths are not distinct')
        return {'bams':{s:str(mapping[s]) for s in trio},'pep':str(pep),'sample_table':str(table)}
    inputs=audit.check('source_identity','input_manifest',input_bams)
    if inputs is None:
        return audit.finish()
    bams={s:Path(p) for s,p in inputs['bams'].items()}
    reference=audit.path(config['ref']['sequence'])
    source_paths=[*audit.config_paths,ped,Path(inputs['pep']),Path(inputs['sample_table']),reference,Path(str(reference)+'.fai')]
    source_paths.extend(bams.values())
    for bam in bams.values():
        index=next((p for p in (Path(str(bam)+'.bai'),bam.with_suffix('.bai'),Path(str(bam)+'.csi')) if p.is_file()),Path(str(bam)+'.bai'))
        source_paths.append(index)
    current=audit.check('provenance','source_metadata_snapshot',lambda:audit.source_snapshot(source_paths))
    if args.capture_sources:
        require(current is not None,'Cannot capture incomplete source snapshot')
        destination=audit.path(args.capture_sources)
        require(not destination.exists(),'Refusing to overwrite prior source snapshot')
        destination.parent.mkdir(parents=True,exist_ok=True);destination.write_text(json.dumps(current,indent=2)+'\n')
        print('Captured source snapshot: '+str(destination));return 0
    if args.source_snapshot:
        def compare_snapshot():
            previous=json.loads(audit.path(args.source_snapshot).read_text())
            before={r['path']:r for r in previous['files']};after={r['path']:r for r in current['files']}
            require(before==after,'Source size,mtime or small-file hash changed; inspect snapshot differences')
            return {'files':len(after),'method':'size+mtime_ns and SHA256 for files within hash limit; not complete BAM content hashing'}
        audit.check('provenance','sources_unchanged_since_capture',compare_snapshot)
    else:
        audit.check('provenance','source_immutability_baseline',lambda:require(False,'No pre-run source snapshot supplied; current identities/statistics recorded without claiming byte immutability'),required=False)
    enable=lambda *keys: bool(config.get(keys[0],{}).get(keys[1],{}).get('enable',False)) if len(keys)==2 else bool(config.get(keys[0],{}).get('enable',False))
    switches={name:enable(*path) for name,path in {'phasing':('phasing',),'JIGV':('JIGV',),'dnSTR':('dnSTR',),'HipSTR':('dnSTR','hipstr'),'dnSV':('dnSV',),'GraphTyper+smoove':('dnSV','GraphTyper2+smoove'),'Manta':('dnSV','Manta'),'rufus':('rufus',),'IGV':('igv_reports',),'MultiQC':('multiqc',),'flagstat':('bam_qc','flagstat'),'WGSmetrics':('bam_qc','collectwgsmetrics'),'multiplemetrics':('bam_qc','collectmultiplemetrics'),'Telomerecat':('bam_qc','telomerecat'),'BAMexports':('bam_qc','bam2bam')}.items()}
    audit.check('configuration','entire_bam_workflow_enabled',lambda:(require(not config.get('fastq_input',{}).get('enable',False),'Expected existing BAM input'),require(all(switches.values()) or args.allow_disabled,'Disabled required branches: '+','.join(k for k,v in switches.items() if not v)),switches)[-1])
    audit.check('runtime','samtools_version',lambda:{'version':audit.sam('--version').splitlines()[0]})
    for sample,bam in bams.items():
        audit.bam('source_identity',bam,sample)
    dnm=audit.vcf('small_variants',out/'dnm_vcf'/f'{family}.dnm.vcf.gz',order=[child],keys=True,indexed=True)
    summary=out/'dnm_vcf_summary/DNM_summary.txt'
    if audit.file('small_variants',summary):
        def dnm_summary():
            rows=noncomment_table(summary);matched=[r for r in rows if r.get('TrioID')==family]
            require(len(matched)==1,'Missing/duplicate DNM family summary')
            count=int(matched[0]['DNM_Cnt']);require(dnm is not None and count==dnm['records_examined'],'DNM summary count differs')
            return {'records':count}
        audit.check('small_variants','summary_matches_vcf',dnm_summary,summary)
    for sample in trio:
        for caller,suffix in [('gatkhc_pb','.gatk.g.vcf.gz'),('deepvariant_pb','.deepvariant.g.vcf.gz')]:
            audit.vcf('small_variants',out/caller/(sample+suffix),order=[sample],full=False,indexed=True)
    for path in [out/'glnexus'/f'{family}.dv_combined.vcf.gz',out/'gatk_cgp'/f'{family}.cgp_norm.vcf.gz']:
        audit.vcf('small_variants',path,samples=trio,full=False,indexed=True)
    for sample in trio:
        qc=config.get('bam_qc',{})
        if enable('bam_qc','flagstat'):
            flag=out/'flagstat'/(sample+'.flagstat');idx=out/'flagstat'/(sample+'.idxstat')
            if audit.file('QC',flag) and audit.file('QC',idx):
                def counts(flag=flag,idx=idx):
                    first=flag.read_text().splitlines()[0];m=re.match(r'(\d+) \+ (\d+) in total',first);require(m is not None,'Unrecognized flagstat header')
                    total=int(m[1])+int(m[2]);rows=[x.split('\t') for x in idx.read_text().splitlines()]
                    require(all(len(x)==4 for x in rows),'Malformed idxstats');n=sum(int(x[2])+int(x[3]) for x in rows)
                    require(total==n and n>0,'flagstat/idxstats counts disagree or empty full BAM')
                    return {'reads':n}
                audit.check('QC','flagstat_idxstats_consistency',counts,flag)
        if enable('bam_qc','collectwgsmetrics'):
            path=out/'collectwgsmetrics'/(sample+'.collect_wgs_metrics.txt')
            if audit.file('QC',path):
                def wgs(path=path):
                    rows=noncomment_table_block(path);require(rows,'Missing WGS metric row');row=rows[0]
                    values={k:finite(row[k],k,0) for k in ['GENOME_TERRITORY','MEAN_COVERAGE','SD_COVERAGE','MEDIAN_COVERAGE','PCT_1X','PCT_10X']}
                    require(values['GENOME_TERRITORY']>0 and values['MEAN_COVERAGE']>0,'Empty WGS coverage')
                    return {'core_metrics':values,'optional_missing_fields':[k for k,v in row.items() if v in MISSING]}
                audit.check('QC','finite_core_wgs_metrics',wgs,path)
        if enable('bam_qc','collectmultiplemetrics'):
            path=out/'collectmultiplemetrics'/sample/'sequencingArtifact.pre_adapter_summary_metrics.txt'
            if audit.file('QC',path):
                def artifacts(path=path,sample=sample):
                    rows=noncomment_table_block(path);require(rows,'No artifact metrics')
                    require({r['SAMPLE_ALIAS'] for r in rows}=={sample},'Artifact metric sample mismatch')
                    qfields=[k for k in rows[0] if 'QSCORE' in k]
                    require(qfields,'Missing artifact quality fields')
                    for row in rows:
                        for key in qfields:finite(row[key],key,0)
                    return {'rows':len(rows),'finite_quality_fields':qfields}
                audit.check('QC','finite_artifact_metrics_and_identity',artifacts,path)
        if enable('bam_qc','telomerecat'):
            path=out/'telomerecat'/(sample+'.csv')
            if audit.file('QC',path):
                def telomere(path=path):
                    rows=noncomment_table(path,',');require(rows,'No Telomerecat data rows')
                    keys=[k for k in rows[0] if 'length' in k.lower()]
                    require(keys,'No length column in Telomerecat CSV')
                    values=[{k:finite(row[k],k,0) for k in keys} for row in rows]
                    return {'rows':len(rows),'length_fields':values}
                audit.check('QC','finite_telomere_length',telomere,path)
        if enable('bam_qc','bam2bam'):
            counts=[]
            for suffix,filtering in [('redup.bam','redup'),('redup.HMQ.sorted.bam','hmq')]:
                value=audit.bam('BAMexports',out/'bam2bam'/(sample+'.triocompass.addRG.'+suffix),sample,filtering)
                counts.append(value)
            if all(counts):
                audit.check('BAMexports','HMQ_count_not_above_redup',lambda c=counts:(require(c[1]['indexed_total_records']<=c[0]['indexed_total_records'],'HMQ count exceeds redup'),{'redup':c[0]['indexed_total_records'],'hmq':c[1]['indexed_total_records']})[1])
    if enable('phasing'):
        path=out/'phase_DNMs'/f'{family}.parental_origin.tab'
        if audit.file('phasing',path):
            def phase():
                rows=noncomment_table(path);keys=[r['VariantID'] for r in rows]
                require(dnm is not None and collections.Counter(keys)==collections.Counter(dnm['variant_keys']),'Phasing variant identities differ from final DNM VCF')
                return {'rows':len(rows),'origins':dict(collections.Counter(r['ParentalOrigin'] for r in rows))}
            audit.check('phasing','phase_table_covers_final_DNMs',phase,path)
    if enable('JIGV'):
        path=out/'call_JIGV'/f'{family}.JIGV.html'
        if audit.file('JIGV',path):audit.check('JIGV','html_document',lambda:html_document(path),path)
    if enable('dnSTR') and config['dnSTR'].get('hipstr',{}).get('enable',False):
        chunks=int(config['dnSTR']['split_n'])
        for i in range(chunks):
            prefix=out/'hipstr'/str(i).zfill(5)
            audit.vcf('STR',Path(str(prefix)+'.vcf.gz'),samples=trio,full=False,indexed=True)
            for suffix in ('.aln.viz.gz','.aln.viz.gz.tbi'):
                audit.file('STR',Path(str(prefix)+suffix))
            for suffix in ('.all_mutations.tab','.locus_summary.tab'):
                audit.file('STR',out/'monstr'/('hipstr_'+str(i).zfill(5)+suffix),empty=True)
        filtered=out/'monstr_filter/hipstr.filtered.tab';summ=out/'dnSTR_summary/hipstr.dnSTR_summary.txt'
        if audit.file('STR',filtered) and audit.file('STR',summ):
            def str_summary():
                rows=noncomment_table(filtered);require(all(r['child'] in trio for r in rows),'Unknown STR child')
                summary=noncomment_table(summ);matched=[r for r in summary if r.get('child')==child and r.get('trio_id')==family]
                require(len(matched)==1 and int(matched[0]['countn'])==sum(r['child']==child for r in rows),'STR summary disagrees with filtered calls')
                return {'filtered_records':len(rows),'trio_records':int(matched[0]['countn'])}
            audit.check('STR','summary_matches_filtered_calls',str_summary,summ)
        recalled=audit.vcf('STR',out/'vizaln'/family/(family+'.hipstr.vcf.gz'),samples=trio,full=True,keys=True,indexed=True)
        audit.file('STR',out/'vizaln'/family/'hipstr/DONE',empty=True)
        markers=list((out/'vizaln'/family/'hipstr/variants').glob('*.dnm'))
        vizaln_loci = audit.check('STR', 'VizAln_locus_coordinates',
            lambda: hipstr_vizaln_loci(out/'vizaln'/family/(family+'.hipstr.vcf.gz'))) if recalled is not None else {}
        if recalled is not None:
            def marker_check():
                expected={'_'.join(key.split(':')[:2]) for key in recalled['variant_keys']}
                require({m.stem for m in markers}==expected,'VizAln markers differ from recalled STR loci')
                return {'loci':len(expected)}
            audit.check('STR','VizAln_marker_coverage',marker_check)
        for marker in markers:
            if audit.file('STR',marker.with_suffix('.html')):
                audit.check('STR','VizAln_html_document',lambda p=marker.with_suffix('.html'):vizaln_html_document(p, trio, (vizaln_loci or {}).get(p.stem), args.max_html_mb),marker.with_suffix('.html'))
    if enable('dnSV') and enable('dnSV','GraphTyper2+smoove'):
        for relative in [f'smoove/{family}/{family}-smoove.genotyped.vcf.gz',f'graphtyper/{family}.gt2.vcf.gz',f'graphtyper_filter/{family}.gt2_filter.vcf.gz']:
            audit.vcf('SV_sample_order',out/relative,order=family_order,full=False,indexed=True)
        for caller in ('smoove','gt2'):
            audit.vcf('SV',out/'dnSV'/f'{family}.{caller}.dnSV.vcf',order=family_order)
        joint=audit.vcf('SV',out/'joint_dnSV'/f'{family}.dnSV.vcf')
        check_count_summary(audit,'SV',out/'dnSV_summary/dnSV_summary.txt',family,joint)
    if enable('dnSV') and enable('dnSV','Manta'):
        audit.vcf('SV',out/'manta_joint'/family/'results/variants/diploidSV.vcf.gz',samples=trio,full=False,indexed=True)
        manta=audit.vcf('SV',out/'dnSV_manta'/f'{family}.dnSV.manta.vcf',samples=trio)
        check_count_summary(audit,'SV',out/'dnSV_manta_summary/dnSV_summary.txt',family,manta,True)
    if enable('rufus'):
        published=out/'rufus'/f'{family}.rufus.vcf'
        rufus=audit.vcf('RUFUS',published,order=trio)
        summary=out/'rufus/rufus_summary.tsv'
        if audit.file('RUFUS',summary):
            def rufus_summary():
                rows=noncomment_table(summary);matched=[r for r in rows if r.get('Family')==family]
                require(len(matched)==1 and rufus is not None and int(matched[0]['TotalRecords'])==rufus['records_examined'],'RUFUS summary count differs')
                return {'records':int(matched[0]['TotalRecords'])}
            audit.check('RUFUS','summary_matches_final_vcf',rufus_summary,summary)
        roots=[audit.path(p) for p in args.rufus_work] or [audit.path(config['rufus'].get('work_dir',str(out/'rufus/work'))),out/'rufus']
        manifests=[]
        for work in roots:
            if work.is_dir():
                for p in work.glob('**/run.json'):
                    try:
                        m=json.loads(p.read_text())
                        if m.get('family')==family and Path(m.get('output','')).resolve()==published.resolve():manifests.append((p,m))
                    except (ValueError,OSError):pass
        def manifest_check():
            complete=[(p,m) for p,m in manifests if m.get('status')=='complete']
            require(complete,'No completed RUFUS attempt manifest for published VCF')
            p,m=max(complete,key=lambda pair:pair[1].get('finished_utc',''))
            require([r['sample'] for r in m['trio']]==trio,'RUFUS manifest trio order mismatch')
            require(rufus is not None and m['counts']['total_records']==rufus['records_examined'],'RUFUS manifest count differs')
            return {'manifest':str(p),'caller_log':m.get('caller_log'),'command':m.get('command'),'work_dir':m.get('work_dir')}
        chosen=audit.check('RUFUS','completed_attempt_provenance',manifest_check)
        if chosen:
            parent=Path(chosen['manifest']).parent
            candidates=sorted(set(parent.glob('*hashcount.fastq.bam.vcf'))|set(parent.glob('*hashcount.fastq.bam.vcf.gz')))
            if not candidates:
                audit.check('rufus_raw','raw_vcf_available',lambda:require(False,'Raw VCF absent from retained attempt; final validation does not establish raw integrity'),required=False)
            for raw in candidates:
                audit.file('rufus_raw',raw,required=False)
                audit.check('rufus_raw','raw_vcf_structure',lambda p=raw:parse_vcf(p),raw,required=False)
    if enable('multiqc'):
        path=audit.path(config['multiqc']['output_dir'])/'multiqc.html'
        if audit.file('MultiQC',path):audit.check('MultiQC','html_document',lambda:html_document(path),path)
    if enable('igv_reports'):
        primary=config['igv_reports'].get('source','joint')
        reports={primary:audit.path(config['igv_reports'].get('output_dir',str(out/'igv_reports')))/family}
        if not args.skip_extra_igv:
            if enable('dnSV','Manta'):reports.setdefault('manta',out/'igv_manta'/family)
            if enable('rufus'):reports.setdefault('rufus',out/'igv_rufus'/family)
        for item in args.extra_igv:
            source,sep,path=item.partition('=');require(sep and source in ('joint','manta','rufus'),'--extra-igv requires joint|manta|rufus=DIR')
            reports[source]=audit.path(path)
        for source,directory in reports.items():audit.report(source,directory,trio,bams,reference)
    return audit.finish()


def source_support(fields):
    genotype=fields.get('GT','NA')
    result={'GT':'NA' if genotype in ('.','./.','.|.') else genotype}
    for key in ('PR','SR'):
        parts=fields.get(key,'.').split(',')
        valid=len(parts)==2 and all(part=='.' or part.isdigit() for part in parts)
        result[key]=' / '.join('NA' if part=='.' else part for part in parts) if valid else 'NA'
    return result


def noncomment_table_block(path):
    """Only first metrics table; Picard may append a separate histogram."""
    header=None;lines=[]
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('#'):continue
            if not line.strip():
                if header is not None:break
                continue
            if header is None:header=line
            lines.append(line)
    rows=list(csv.DictReader(lines,delimiter='\t'))
    require(all(None not in row and all(v is not None for v in row.values()) for row in rows),'Ragged metric table')
    return rows


def hipstr_vizaln_loci(path):
    """Map POS-named artifacts to HipSTR's actual repetitive START/END interval."""
    loci = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            require(len(fields) >= 8, 'Malformed recalled HipSTR VCF record')
            info = collections.defaultdict(list)
            for item in fields[7].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info[key].append(value)
            require(all(len(info[key]) == 1 and info[key][0].isdigit() for key in ('START', 'END')),
                    'HipSTR visualization requires unique integer INFO/START and INFO/END')
            start, end = int(info['START'][0]), int(info['END'][0])
            require(0 < start <= end, 'Invalid HipSTR visualization interval')
            key = fields[0] + '_' + fields[1]
            locus = [fields[0], start, end]
            require(key not in loci or loci[key] == locus, 'Ambiguous HipSTR interval for ' + key)
            loci[key] = locus
    return loci


class VizAlnFragmentParser(HTMLParser):
    """Validate the emitted alignment-table dialect, including explicit closures.

    HipSTR generate_aln_html.py emits a style block and a div/table fragment,
    without html/doctype. HTMLParser alone tolerates truncation, so maintain an
    explicit stack and verify reference, locus, sample and read-row content.
    """
    allowed = {'html', 'head', 'body', 'style', 'div', 'table', 'caption', 'tr', 'td', 'font', 'br', 'meta'}
    void = {'br', 'meta'}
    read_classes = {'Atd', 'Ctd', 'Gtd', 'Ttd', 'Ntd', '-td', 'vtd', 'spacertd', 'inserttd', 'snptd'}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.tags = collections.Counter()
        self.styles = []
        self.captions = []
        self.samples = []
        self.read_rows = collections.Counter()
        self.reference_cells = 0
        self.current_sample = None
        self.cell = None
        self.row_read_cells = 0
        self.row_reference_cells = 0
        self.row_sample_label = False

    def handle_starttag(self, tag, attrs):
        require(tag in self.allowed, 'Unexpected VizAln tag: ' + tag)
        attributes = dict(attrs)
        require(not any(key.lower().startswith('on') for key in attributes), 'Unexpected executable HTML attribute')
        if tag == 'table':
            require('reftable' in attributes.get('class', '').split(), 'Missing VizAln reference table class')
        if tag in ('caption', 'tr'):
            require(self.stack and self.stack[-1] == 'table', 'VizAln ' + tag + ' outside table')
        if tag == 'td':
            require(self.stack and self.stack[-1] == 'tr', 'VizAln cell outside row')
            self.cell = {'classes': set(attributes.get('class', '').split()), 'text': []}
        if tag == 'tr':
            self.row_read_cells = self.row_reference_cells = 0
            self.row_sample_label = False
        if tag == 'style':
            self.styles.append('')
        if tag == 'caption':
            self.captions.append('')
        self.tags[tag] += 1
        if tag not in self.void:
            self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        if tag not in self.void:
            self.handle_endtag(tag)

    def handle_data(self, data):
        if self.stack and self.stack[-1] == 'style':
            self.styles[-1] += data
        if 'caption' in self.stack:
            self.captions[-1] += data
        if self.cell is not None:
            self.cell['text'].append(data)

    def handle_endtag(self, tag):
        require(self.stack and self.stack[-1] == tag, 'Unbalanced/truncated VizAln tag: ' + tag)
        if tag == 'td':
            classes, text = self.cell['classes'], ''.join(self.cell['text']).strip()
            if 'samplename' in classes:
                label, separator, genotype = text.rpartition(':')
                require(separator and label.strip() and genotype.strip(), 'Malformed VizAln sample label')
                self.current_sample = label.strip()
                self.samples.append(self.current_sample)
                self.row_sample_label = True
            elif any(re.fullmatch(r'[ACGTNv]reftd', c) for c in classes):
                require(re.fullmatch(r'[ACGTN*]', text) is not None, 'Invalid VizAln reference cell')
                self.row_reference_cells += 1
                self.reference_cells += 1
            elif classes & self.read_classes:
                require(re.fullmatch(r'[ACGTN*x-]', text) is not None, 'Invalid VizAln alignment cell')
                self.row_read_cells += 1
            self.cell = None
        if tag == 'tr' and self.row_read_cells:
            require(self.current_sample is not None and not self.row_sample_label and not self.row_reference_cells,
                    'VizAln alignment row lacks a distinct sample section')
            self.read_rows[self.current_sample] += 1
        self.stack.pop()


def vizaln_html_document(path, samples, expected_locus, max_html_mb=256):
    path = Path(path)
    require(path.stat().st_size <= max_html_mb * 1024 * 1024, 'VizAln report exceeds HTML validation bound')
    require(expected_locus is not None, 'No recalled HipSTR interval for VizAln report')
    document = path.read_text(encoding='utf-8')
    require('\x00' not in document, 'NUL byte in VizAln report')
    parser = VizAlnFragmentParser()
    parser.feed(document)
    parser.close()
    require(not parser.stack, 'Unclosed/truncated VizAln markup: ' + ','.join(parser.stack))
    require(parser.tags['table'] == 1 and len(parser.captions) == 1, 'Missing/surplus VizAln alignment table or locus caption')
    require(parser.styles and '.reftable' in ''.join(parser.styles), 'Missing VizAln alignment styling')
    parts = parser.captions[0].split()
    require(len(parts) == 3 and all(p.isdigit() for p in parts[1:]), 'Malformed VizAln locus caption')
    locus = [parts[0], int(parts[1]), int(parts[2])]
    require(locus == list(expected_locus), 'VizAln caption {} != recalled HipSTR interval {}'.format(locus, expected_locus))
    require(parser.reference_cells > 0, 'Missing VizAln reference sequence cells')
    require(parser.samples == list(samples), 'VizAln sample sections {} != required trio {}'.format(parser.samples, samples))
    require(all(parser.read_rows[s] > 0 for s in samples), 'Missing VizAln alignment rows for one or more trio members')
    return {'html': True, 'fragment': not bool(parser.tags['html']), 'bytes': path.stat().st_size,
            'locus': locus, 'samples': parser.samples, 'reference_cells': parser.reference_cells,
            'alignment_rows_by_sample': dict(parser.read_rows), 'rendered_in_browser': False}


def html_document(path):
    with Path(path).open(errors='replace') as handle:prefix=handle.read(65536).lower()
    require('<html' in prefix or '<!doctype html' in prefix,'Missing HTML header')
    return {'html':True,'bytes':Path(path).stat().st_size,'rendered_in_browser':False}


def check_count_summary(audit,branch,path,family,vcf,manta=False):
    if not audit.file(branch,path):return
    def inspect():
        rows=noncomment_table(path);matched=[r for r in rows if r.get('TrioID')==family]
        require(len(matched)==1 and vcf is not None,'Missing family summary or VCF')
        column=next(k for k in matched[0] if k!='TrioID')
        count=int(matched[0][column]);expected=vcf['records_examined']
        bnd=vcf['svtype_counts'].get('BND',0)
        if manta:expected-=bnd//2
        require(count==expected,'SV summary count differs from its implemented record/event convention')
        return {'summary_count':count,'source_records':vcf['records_examined'],'BND_source_rows':bnd,'BND_mate_pairing_independently_validated':False}
    audit.check(branch,'summary_matches_call_records',inspect,path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',required=True)
    parser.add_argument('--config',action='append',help='Merge complete YAML configs in order; default RUN/config/config.yaml')
    parser.add_argument('--family',default='GIAB')
    parser.add_argument('--ped');parser.add_argument('--output-dir');parser.add_argument('--report-dir')
    parser.add_argument('--samtools');parser.add_argument('--rufus-work',action='append',default=[])
    parser.add_argument('--extra-igv',action='append',default=[],metavar='SOURCE=DIR')
    parser.add_argument('--skip-extra-igv',action='store_true',help='Validate only the primary configured report before extra report stages')
    parser.add_argument('--allow-disabled',action='store_true',help='Audit configured scope without claiming all workflow branches enabled')
    parser.add_argument('--capture-sources',metavar='JSON',help='Capture pre-run source stats/small-file hashes only; refuses overwrite')
    parser.add_argument('--source-snapshot',metavar='JSON',help='Compare post-run sources to --capture-sources JSON')
    parser.add_argument('--intermediate-records',type=int,default=50)
    parser.add_argument('--bam-records',type=int,default=1000)
    parser.add_argument('--report-sessions',type=int,default=12)
    parser.add_argument('--max-html-mb',type=int,default=256)
    parser.add_argument('--hash-limit-mb',type=int,default=16)
    args=parser.parse_args()
    require(min(args.intermediate_records,args.bam_records,args.report_sessions,args.max_html_mb,args.hash_limit_mb)>0,'Validation bounds must be positive')
    audit=Audit(args)
    try:return run(audit)
    except Exception as error:
        audit.check('validator','unexpected_validation_error',lambda:(_ for _ in ()).throw(error))
        return audit.finish()


if __name__=='__main__':
    sys.exit(main())
