TriosCompass expects to use the parent folder of TriosCompass_v2 (the repo clone folder) as the working space, so as to separate the Snakemake workflow from the working space.  

Three configure files are required:
1. *Profile yaml file*: workflow/profiles/<PROFILE_NAME>/config.yaml
2. *Config yaml file* 
3. *Sample yaml file* 

The *config yaml file* can be specified by Snakemake command-line argument "\-\-configfile" or in the "*profile yaml file*, for example:

```yml
configfile: TriosCompass_v2/config/fullbam_config.yaml
snakefile: TriosCompass_v2/workflow/Snakefile
```

In turn, the *sample yaml file* is specified in the *config yaml file*, to [define sample input via PEPs](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#configuring-scientific-experiments-via-peps)

```yml
pepfile: "config/fullbam_pep.yaml"
pepschema: "../schemas/bam_schema.yaml"
```

+ Example of the PEP configure file for FASTQ input files
  + config/fastq_pep.yaml
    ```yml
    pep_version: 2.0.0
    sample_table: sample_fastq.csv
    
    # In manifest file, Sample_ID + Flowcell should be unique
    sample_modifiers:
      append:
        sample_name: "sn"
      derive:
        attributes: [sample_name]
        sources:
          sn: "{SAMPLE_ID}_{FLOWCELL}"
    ```

  + config/sample_fastq.csv 
    ```csv
    SAMPLE_ID,FLOWCELL,LANE,INDEX,R1,R2
    HG002,BH2JWTDSX5,1,CGGTTGTT-GTGGTATG,data/fq/HG002_NA24385_son_80X_R1.fq.gz,data/fq/HG002_NA24385_son_80X_R2.fq.gz
    HG003,BH2JWTDSX5,1,GCGTCATT-CAGACGTT,data/fq/HG003_NA24149_father_80X_R1.fq.gz,data/fq/HG003_NA24149_father_80X_R2.fq.gz
    HG004,BH2JWTDSX5,1,CTGTTGAC-ACCTCAGT,data/fq/HG004_NA24143_mother_80X_R1.fq.gz,data/fq/HG004_NA24143_mother_80X_R2.fq.gz
    ```
  + workflow/schemas/fastq_schema.yaml 
(schemas to validate config/fastq_pep.yaml) 
    ```yml
    workflow/schemas/fastq_schema.yaml 
    description: A example schema for a pipeline.
    imports:
      - http://schema.databio.org/pep/2.0.0.yaml
      # - TriosCompass_v2/workflow/schemas/2.0.0.yaml
      
    properties:
      samples:
        type: array
        items:
          type: object
          properties:
            SAMPLE_ID:
              type: string
              description: "sample id"
            FLOWCELL:
              type: string
              description: "Flowcell"
            INDEX:
              type: string
              description: "Library index"
            LANE:
              type: string
              description: "Lane number in flowcell"
              enum: ["1", "2"]
            R1:
              type: string
              description: "path to the R1 fastq file"
            R2:
              type: string
              description: "path to the R2 fastq file"
          required:
            - FLOWCELL
            - SAMPLE_ID
            - INDEX
            - R1
            - R2
    ```
+ Example of the PEP configure file for BAM input files
  + config/bam_pep.yaml
    ```yml
    pep_version: 2.0.0
    sample_table: sample_bam.csv
    
    
    sample_modifiers:
        append:
        sample_name: "sn"
        derive:
        attributes: [sample_name]
        sources:
            sn: "{SAMPLE_ID}"
    ```

  + config/sample_bam.csv 
    ```csv
    SAMPLE_ID,BAM
    HG002,sorted_bam/HG002_NA24385_son_80X.bam
    HG003,sorted_bam/HG003_NA24149_father_80X.bam
    HG004,sorted_bam/HG004_NA24143_mother_80X.bam
    ```
  + workflow/schemas/fastq_schema.yaml 
(schemas to validate config/bam_pep.yaml) 
    ```yml
    description: A example schema for a pipeline.
    imports:
      - http://schema.databio.org/pep/2.0.0.yaml
      # - TriosCompass_v2/workflow/schemas/2.0.0.yaml
      
    properties:
      samples:
        type: array
        items:
          type: object
          properties:
            SAMPLE_ID:
              type: string
              description: "sample id"
            BAM:
              type: string
              description: "path to the bam file"
          required:
            - SAMPLE_ID
            - BAM
    ```

## Portable v1.1.0 migration

`config/config.yaml` supplies defaults. Pass a run-specific YAML with
`--configfile`; Snakemake merges it over those defaults. Replace all example
BAM/FASTQ, reference, BED and pedigree paths before running. The new example
PEPs/PED are synthetic: `example_bam_pep.yaml`, `rerf_fastq_pep.yaml`, and
`example_ped/FAMILY.ped`. No study manifest is included in this migration.

The FASTQ schema `workflow/schemas/rerf_fastq_schema.yaml` is self-contained
and allows lanes 1–8. `pepfile` is relative to the repository root (absolute
paths also work); sample-table paths are relative to the PEP file. Set
`fastq_input.enable: true`, the FASTQ PEP, and
`pepschema: ../schemas/rerf_fastq_schema.yaml` for FASTQ input. BAM paths in
sample tables, reference paths and output paths are relative to the run working
directory unless absolute.

The shared SV exclusion BED is now `dnSV.exclude_bed`. A nonempty legacy
`dnSV.GraphTyper2+smoove.exclude_bed` remains supported when the top-level value
is absent or null. Global `dnSV.enable` gates the existing Manta/GraphTyper2/smoove
branches. Rufus and IGV have independent enable flags. `all` is the default target
and works with MultiQC disabled.

### RUFUS

```yaml
rufus:
  enable: true
  executable: /path/to/RUFUS/runRufus.sh
  # For a container, executable names the program INSIDE that image.
  container: null
  runtime: apptainer
  runtime_args: []
  extra_args: []
  kmer_size: 25
  threads: 24
  mem_mb: 256000
  runtime_minutes: 5760
```

RUFUS receives the workflow reference and PED-derived child/father/mother BAMs.
The helper checks BAM sample names, stages distinct filenames plus indexes,
and retains each attempt under `{output_dir}/rufus/work/{family}/attempt-*`.
The reference needs `.fai` and BWA `.amb/.ann/.bwt/.pac/.sa` indexes; the workflow
provides indexing dependencies. Use the correctness-fixed RUFUS installation
validated for your inputs; this PR integrates the caller, not its source patch.
Record the chosen build/version in the run configuration or provenance.

Successful calls are published to `rufus/{family}.rufus.vcf`; the raw caller
VCF, log and `run.json` remain with the attempt. The exported VCF uses PED sample
IDs in verified child/father/mother order. `rufus/rufus_summary.tsv` reports total
records, allele classes and source `SVTYPE` counts. Total records are **not** an
SV-only count or a count of validated de novo events. Failed attempts do not
publish a new final VCF. Reruns retain previous work; large k-mer databases need
an explicit retention/cleanup policy for production cohorts.

The default CPU/memory/time values are starting requests for WGS, not measured
resource guarantees. Override them for your data. Native RUFUS must have all its
external dependencies installed. `container` may name a local image or a runtime
URI; input/reference parents are bound read-only and the attempt directory is
bound read-write. Additional installation binds go in `runtime_args`.

### Telomerecat and filtered BAM exports

```yaml
bam_qc:
  telomerecat:
    enable: true
    threads: 2
    mem_mb: 16000
    runtime_minutes: 720
  bam2bam:
    enable: true
    threads: 4
    mem_mb: 16000
    runtime_minutes: 720
  collectmultiplemetrics:
    enable: true
    mode: seq_artifact  # use all to retain the previous all-metrics behavior
containers:
  telomerecat: docker://quay.io/wtsicgp/telomerecat:4.0.0
```

Telomerecat runs `bam2length` separately for each BAM, passing the allocated
process count. CSVs are written to `telomerecat/{sample}.csv`; per-attempt working
files live under `telomerecat/work/{sample}`. This reproduces the reference's
per-sample analysis, not pooled cohort correction. CSVs are workflow outputs;
they are not automatically incorporated into MultiQC. The pinned CLI and image
are documented by [Telomerecat](https://github.com/cancerit/telomerecat/tree/4.0.0).

`bam2bam` exports two indexed BAMs under `bam2bam/`: one removes duplicate-flagged
reads; the other also requires a proper pair, MAPQ >=60, excludes secondary
alignments, and accepts CIGARs containing M but no S/H/P/N. Supplementary records
are retained when they meet those filters, matching the reference. Inputs must
be coordinate sorted. Despite the legacy `.addRG` filename, this step does not
change read groups. These are additional products and do not replace the BAMs
used by variant callers or Telomerecat. The filters remove evidence useful for
SV review, so IGV uses the original workflow BAMs.

### Interactive SV reports with igv-reports

```yaml
igv_reports:
  enable: true
  source: joint  # joint, manta or rufus
  igv_js: /path/to/igv.3.5.2.min.js
  flank: 2000
  initial_window: 400
  coverage_max: 200
  sampling_depth: 500
  min_sv_length: 50
  max_variants: 0  # all eligible records; set a limit explicitly if needed
  mem_mb: 4000
  runtime_minutes: 60
```

Enable the selected caller when its VCF still needs to be generated. Existing
VCFs can be supplied through `igv_reports.vcf_template`, for example
`/results/{fam}.vcf.gz`, without enabling that caller. Use `source: joint` for
SURVIVOR/caller-summary files; the viewer suppresses genotype/support attribution
for those columns. `vcf_samples_are_trio: true` can explicitly enable exact-ID
mapping for a custom true multisample VCF. Missing GT/PR/SR remain NA, never an
inferred reference genotype or zero support. Parent evidence absent from the
input VCF cannot be recovered from that VCF.

Download [IGV.js 3.5.2](https://cdn.jsdelivr.net/npm/igv@3.5.2/dist/igv.min.js)
once to persistent software storage and supply its local path. It is embedded
at generation time; the resulting HTML can be opened offline. The report records
the asset hash. The pinned `workflow/envs/igv-reports.yaml` environment is selected
with `--use-conda`. Runtime images similarly require `--use-singularity` or an
already provisioned native tool environment.

Open `{output_dir}/igv_reports/{family}/index.html`. The report embeds the reference
and child/father/mother read slices and provides zoom, pan, reset, read-size and
track-height controls. DEL/DUP/INV/BND use paired anchors; INS uses its insertion
anchor. Both panels retain reference-forward orientation. Browsing is limited to
embedded windows; long-event interiors are not included. Source END/SVLEN and
raw records are retained, and unsupported/malformed records are listed in
`skipped_records.tsv`. Unannotated sequence alleles shorter than `min_sv_length`
are skipped; explicit supported SVTYPE annotations remain eligible regardless of
length. Empty inputs produce a valid empty report. Manifest/source hashes and
displayed support are exported alongside the HTML.

This replaces the experimental Samplot workflow. IGV reports are inspection
artifacts, not independent validation of germline de novo status.

### Runtime changes and Biowulf

The main environment now matches the validated Snakemake 7.32.4/Python 3.11
runtime, with `pulp<2.8` and `setuptools==69.1.1` for Snakemake 7 compatibility
(including its `pkg_resources` version check), and Peppy 0.40.4.
TRTools retains 5.0.2 with Python 3.10,
SciPy 1.11.4 and `libdeflate<1.25`. FastQC defaults to the reference's
`v1.3.2/bio/fastqc` wrapper and four threads; override
`fastq_input.fastqc.wrapper` or `.threads` when needed. HipSTR outputs and indexes
are staged beside their configured destination and published only after the
commands succeed. GATK normalization also works with hard filtering disabled.

`containers` supports `samtools`, `picard`, `telomerecat`, `hipstr`, `manta`,
`manta_joint` and `multiqc`; replace defaults with locally tested immutable
images/digests as appropriate. The original HipSTR image still uses `latest`
because the supplied sandbox did not identify an equivalent public release.
For a different joint-Manta image set `dnSV.Manta.config_command` to its actual
`configManta.py` path. `parabricks.singularity_cmd` remains the configurable
runtime command; the default image is 4.3.1-1 and
`parabricks.fq2bam_extra` defaults to `--num-gpus 1 --low-memory`.

Use the new Snakemake 7 profile from the repository root:

```bash
snakemake --snakefile workflow/Snakefile --configfile /path/to/run.yaml --dry-run all
snakemake --profile workflow/profiles/biowulf --configfile /path/to/run.yaml all
```

The second command submits jobs. The existing `slurm` profile remains available
for its original site settings. The Biowulf profile uses `gpu` rather than
`gpuq`, submits at bounded rates and preserves the new rules' configurable
resource requests. GPU declarations are translated into a Slurm GRES request,
while explicit GPU/GRES choices remain authoritative. Configure container binds
for all input/output files; cache environments and images in persistent storage.

Unit tests use synthetic data. Full DAG checks and mocked launchers do not
establish whole-genome caller concordance, telomere accuracy, or GPU performance.


## Container-only execution
The supplied Slurm jobscript forwards `CUDA_VISIBLE_DEVICES` into Apptainer/Singularity, including an explicitly empty assignment, so clean container environments retain the scheduler GPU allocation.


Use `workflow/profiles/biowulf-container` to execute scientific tools through
Apptainer/Singularity with `--cleanenv`. This profile disables host Conda and
environment modules and forces `container_only=true`. Workflow parsing rejects
missing images, native executable rules, and mixed Conda execution before any
analysis is submitted. The Snakemake controller and scheduler helpers still run
in their orchestration environment.

Provide local SIF paths in the run configuration. Container files should contain
the tools and their dependencies; binding a host Conda prefix into a generic
image does not make a self-contained runtime. The tested GIAB images package the
verified environments at `/opt/env`, with `PATH`, `CONDA_PREFIX`, and (for GATK)
`JAVA_HOME` set inside the image. GNU Parallel and `bc` are included where used.

```yaml
runtime:
  container_only: true
  executable: apptainer  # singularity is also supported
  args: [--cleanenv]
containers:
  utility: /absolute/path/utility.sif
  bcftools: /absolute/path/bcftools.sif
  tabix: /absolute/path/tabix.sif
  gatk4: /absolute/path/gatk4.sif
  glnexus: /absolute/path/glnexus.sif
  slivar: /absolute/path/slivar.sif
  csvtk: /absolute/path/csvtk.sif
  graphtyper: /absolute/path/graphtyper.sif
  vcflib: /absolute/path/vcflib.sif
  survivor: /absolute/path/survivor.sif
  whatshap: /absolute/path/whatshap.sif
  trtools: /absolute/path/trtools.sif
  igv-reports: /absolute/path/igv-reports.sif
  rufus: /absolute/path/rufus.sif
  bwa: /absolute/path/bwa.sif
parabricks:
  singularity_cmd: apptainer exec --cleanenv --nv /absolute/path/parabricks.sif
rufus:
  enable: true
  executable: /opt/rufus/runRufus.sh
  container: null
```

Configure only the optional branches you enable, plus the common calling tools.
Existing container tools also accept local images through `samtools`, `picard`,
`telomerecat`, `hipstr`, `manta`, `manta_joint`, `smoove`, `svimmer`, `monstr`,
`jigv`, and `multiqc`. FASTQ processing additionally requires `fastp` and, when
enabled, `fastqc`. BWA is required for FASTQ alignment or RUFUS reference indexing.
The utility image must include Python/Python3, Perl, GNU Parallel, `bc`, and standard
GNU shell utilities. A GraphTyper image must include GNU Parallel.
The HipSTR image must include Python 3, `VizAln_rev` on `PATH`, and its
adjacent `generate_aln_html.py` helper. STR report generation uses the VCF
`INFO/START` and `INFO/END` coordinates while retaining VCF `POS` filenames.

Parabricks runs through its explicit GPU container command. Its `bgzip` and `tabix`
steps use a separate configured HTS image, avoiding a container within a container.
`containers.rufus` contains the complete rule, including its Python adapter; do
not also set `rufus.container`, which is the older caller-only mode. Set the
RUFUS executable to its path inside the image.

```bash
snakemake --snakefile workflow/Snakefile \
  --profile workflow/profiles/biowulf-container \
  --configfile /absolute/path/run.yaml all
```

Supply filesystem binds through `APPTAINER_BINDPATH`/`SINGULARITY_BINDPATH` in
the launch environment and jobscript. On Biowulf include both `/data` and its
`/vf` targets and allocated scratch. Forward intentional thread/RUFUS tuning
variables with `APPTAINERENV_`/`SINGULARITYENV_` prefixes; `--cleanenv` removes
ordinary host variables, including conflicting Python packages.

The IGV environment pins Python 3.11, pysam 0.23.3, and igv-reports 1.17.0.
The 1.17.0 package was unavailable in the configured Conda channels, so the
recipe installs that same pinned release through pip. The container packages
that verified environment; workflow execution requires no pip/Conda installation.
