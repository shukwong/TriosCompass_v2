# Validate workflow outputs

`scripts/validate_full_workflow.py` audits TriosCompass output structure, sample
identity, report contents and source provenance. It supports one complete trio
from an existing-BAM workflow. It does not run callers or modify scientific
outputs.

Use Python 3.9 or newer, PyYAML, and `samtools` on `PATH`. Pass `--samtools` to
select another executable or a wrapper that invokes samtools in your container.
The focused validator tests use only the Python standard library.

Run from this repository, pointing `--run-dir` at the directory containing the
workflow's configuration and outputs:

```bash
python scripts/validate_full_workflow.py \
  --run-dir /path/to/run \
  --family FAM001
```

The default family is `GIAB`; specify `--family` for other pedigrees. The auditor
reads `RUN/config/config.yaml`, the configured PEP/sample table, and the family's
PED file. It resolves relative configuration paths against `--run-dir`; a PEP
`sample_table` path is relative to the PEP file. Use `--config` repeatedly to
merge configuration files in order. `--ped` and `--output-dir` override their
configured locations.

By default, the audit requires all supported BAM workflow branches to be
enabled. For an intentionally smaller configured workflow, add
`--allow-disabled`; this audits the enabled scope without claiming that every
branch ran. Use `--skip-extra-igv` before the additional Manta/RUFUS report stages,
or `--extra-igv manta=/path/to/reports/FAM001` to specify an extra report directory.

## Reports and exit status

Each invocation writes or replaces these files in `RUN/validation/`, or in the
directory selected with `--report-dir`:

- `validation.json`: check results, artifact metadata and explicit scope limits.
- `checks.tsv`: one row per check.
- `summary.txt`: required failures and separate warnings.

Exit status is `0` when required checks pass and `1` otherwise. Missing or
malformed required outputs fail validation. Warnings remain visible without
turning into required failures. In particular, a malformed retained raw RUFUS
VCF is reported separately from the published final VCF's validation.

The audit checks VCF/sample structure, indexes, bounded BAM retrieval, QC tables,
count summaries, phasing coverage and configured visualization reports. HipSTR
VizAln emits an HTML fragment without an `<html>` envelope. The dedicated check
requires a balanced alignment table, reference cells, ordered trio sample
sections and read rows. Its caption must match HipSTR's `INFO/START` and
`INFO/END`, which can differ from the VCF `POS` used in the report filename.
Stylesheet-only or truncated fragments fail; a generic HTML header is not enough.

## Source provenance and bounds

Capture source identities before execution, then compare them during the final
audit:

```bash
python scripts/validate_full_workflow.py \
  --run-dir /path/to/run --family FAM001 \
  --capture-sources provenance/source_snapshot.json

python scripts/validate_full_workflow.py \
  --run-dir /path/to/run --family FAM001 \
  --source-snapshot provenance/source_snapshot.json
```

Capture refuses to overwrite an existing snapshot. Snapshots record file sizes,
modification times and SHA256 hashes for files within `--hash-limit-mb` (default
16 MB). They do not establish full BAM/reference byte identity. Without a prior
snapshot, the audit records current identities and a provenance warning.

Large intermediate VCF and BAM checks use bounded samples. The relevant options
are `--intermediate-records`, `--bam-records`, `--report-sessions`,
`--max-html-mb` and `--hash-limit-mb`; `--help` lists their defaults. Final VCF
checks scan all records. Reports are inspected structurally, not rendered in a
browser.

Artifact validation does not establish scheduler completion, whole-genome call
set equivalence, biological validity, or truth-set accuracy. Reconcile workflow
and scheduler status separately; precision/recall requires a separate benchmark.

Run the portable fragment regressions with:

```bash
python -m unittest discover -s tests -p 'test_full_workflow_validation.py' -v
```
