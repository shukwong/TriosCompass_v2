# Optional sample QC. Filtered BAMs are exports, never implicit caller inputs.
_bam_qc = config.get("bam_qc", {})
_samtools_image = config.get("containers", {}).get("samtools", "docker://euformatics/samtools:1.19.2")

if _bam_qc.get("flagstat", {}).get("enable", False):
    rule flagstat:
        input: get_bam_by_subj
        output:
            flagstat=output_dir + "/flagstat/{subj}.flagstat",
            idxstat=output_dir + "/flagstat/{subj}.idxstat"
        threads: config.get("threads", {}).get("flagstat", 2)
        params: workers=lambda w, threads: max(0, threads - 1)
        singularity: _samtools_image
        shell:
            "samtools flagstat -@ {params.workers} {input:q} > {output.flagstat:q}; "
            "samtools idxstats {input:q} > {output.idxstat:q}"
    optional_output.extend([expand(output_dir + "/flagstat/{subj}.flagstat", subj=subjs), expand(output_dir + "/flagstat/{subj}.idxstat", subj=subjs)])
    qc_output.append(expand(output_dir + "/flagstat/{subj}.idxstat", subj=subjs))

if _bam_qc.get("collectwgsmetrics", {}).get("enable", False):
    rule collectwgsmetrics:
        input: bam=get_bam_by_subj, ref=genome
        output: output_dir + "/collectwgsmetrics/{subj}.collect_wgs_metrics.txt"
        threads: config.get("threads", {}).get("collectwgsmetrics", 1)
        resources:
            mem_mb=_bam_qc["collectwgsmetrics"].get("mem_mb", 64000),
            runtime=_bam_qc["collectwgsmetrics"].get("runtime_minutes", 1440)
        params: heap_mb=lambda w, resources: max(1, int(resources.mem_mb * 0.8))
        benchmark: output_dir + "/benchmark/collectwgsmetrics/{subj}.tsv"
        singularity: config.get("containers", {}).get("picard", "docker://quay.io/biocontainers/picard:2.27.3--hdfd78af_0")
        shell:
            """
            unset _JAVA_OPTIONS JDK_JAVA_OPTIONS
            export JAVA_TOOL_OPTIONS="-Xmx{params.heap_mb}m -XX:+UseSerialGC"
            picard CollectWgsMetrics I={input.bam:q} O={output:q} R={input.ref:q}
            """
    optional_output.append(expand(output_dir + "/collectwgsmetrics/{subj}.collect_wgs_metrics.txt", subj=subjs))
    qc_output.append(expand(output_dir + "/collectwgsmetrics/{subj}.collect_wgs_metrics.txt", subj=subjs))

if _bam_qc.get("collectmultiplemetrics", {}).get("enable", False):
    _metrics_mode = _bam_qc["collectmultiplemetrics"].get("mode", "seq_artifact")
    if _metrics_mode not in ("seq_artifact", "all"):
        raise ValueError("bam_qc.collectmultiplemetrics.mode must be seq_artifact or all")
    rule collectmultiplemetrics:
        input: bam=get_bam_by_subj, ref=genome
        output: output_dir + "/collectmultiplemetrics/{subj}/sequencingArtifact.pre_adapter_summary_metrics.txt"
        threads: config.get("threads", {}).get("collectmultiplemetrics", 8)
        resources: gpu=1
        benchmark: output_dir + "/benchmark/collectmultiplemetrics/{subj}.tsv"
        params:
            prefix=output_dir + "/collectmultiplemetrics/{subj}",
            singularity_cmd=config["parabricks"]["singularity_cmd"],
            mode="--gen-all-metrics" if _metrics_mode == "all" else "--gen-seq-artifact"
        shell:
            "{params.singularity_cmd} pbrun collectmultiplemetrics --ref {input.ref:q} "
            "--bam {input.bam:q} --out-qc-metrics-dir {params.prefix:q} {params.mode}"
    optional_output.append(expand(output_dir + "/collectmultiplemetrics/{subj}/sequencingArtifact.pre_adapter_summary_metrics.txt", subj=subjs))
    qc_output.append(expand(output_dir + "/collectmultiplemetrics/{subj}/sequencingArtifact.pre_adapter_summary_metrics.txt", subj=subjs))

if _bam_qc.get("telomerecat", {}).get("enable", False):
    _telomerecat = _bam_qc["telomerecat"]
    rule telomerecat:
        input: bam=get_bam_by_subj
        output: csv=output_dir + "/telomerecat/{subj}.csv"
        threads: _telomerecat.get("threads", config.get("threads", {}).get("telomerecat", 2))
        resources:
            mem_mb=_telomerecat.get("mem_mb", 16000),
            runtime=_telomerecat.get("runtime_minutes", 720)
        params:
            bam=lambda w, input: os.path.abspath(input.bam),
            csv=lambda w, output: os.path.abspath(output.csv),
            work=lambda w: os.path.abspath(output_dir + "/telomerecat/work/" + w.subj),
            executable=_telomerecat.get("executable", "telomerecat")
        singularity: config.get("containers", {}).get("telomerecat", "docker://quay.io/wtsicgp/telomerecat:4.0.0")
        log: output_dir + "/logs/telomerecat/{subj}.log"
        benchmark: output_dir + "/benchmark/telomerecat/{subj}.tsv"
        shell:
            """
            export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
            mkdir -p {params.work:q}
            work=$(mktemp -d {params.work:q}/attempt.XXXXXX)
            log_file=$(realpath {log:q})
            cd "$work"
            echo "Work directory: $work" > "$log_file"
            {params.executable:q} bam2length -p {threads} -v2 --output "$work/length.csv" {params.bam:q} >> "$log_file" 2>&1
            test -s "$work/length.csv"
            mv "$work/length.csv" {params.csv:q}
            """
    optional_output.append(expand(output_dir + "/telomerecat/{subj}.csv", subj=subjs))

if _bam_qc.get("bam2bam", {}).get("enable", False):
    _bam2bam = _bam_qc["bam2bam"]
    rule bam2bam:
        input: bam=get_bam_by_subj
        output:
            redup_bam=output_dir + "/bam2bam/{subj}.triocompass.addRG.redup.bam",
            redup_idx=output_dir + "/bam2bam/{subj}.triocompass.addRG.redup.bam.bai",
            hmq_bam=output_dir + "/bam2bam/{subj}.triocompass.addRG.redup.HMQ.sorted.bam",
            hmq_idx=output_dir + "/bam2bam/{subj}.triocompass.addRG.redup.HMQ.sorted.bam.bai"
        threads: _bam2bam.get("threads", config.get("threads", {}).get("bam2bam", 4))
        resources:
            mem_mb=_bam2bam.get("mem_mb", 16000),
            runtime=_bam2bam.get("runtime_minutes", 720)
        params:
            workers=lambda w, threads: max(0, threads - 1),
            sort_mem=lambda w, threads, resources: max(1, int(resources.mem_mb * 0.7 / max(1, threads))),
            work=lambda w: output_dir + "/bam2bam/work/" + w.subj
        singularity: _samtools_image
        log: output_dir + "/logs/bam2bam/{subj}.log"
        benchmark: output_dir + "/benchmark/bam2bam/{subj}.tsv"
        shell:
            """
            exec 2> {log:q}
            mkdir -p {params.work:q}
            work=$(mktemp -d {params.work:q}/attempt.XXXXXX)
            echo "Work directory: $work" >&2
            samtools view -@ {params.workers} -b -F 0x400 -o "$work/redup.bam" {input.bam:q}
            samtools index -@ {params.workers} "$work/redup.bam"
            samtools view -@ {params.workers} -b -f 0x2 -q 60 -F 0x500 -e 'cigar !~ "[SHPN]" && cigar =~ "M"' -o "$work/selected.bam" "$work/redup.bam"
            samtools sort -m {params.sort_mem}M -@ {params.workers} -T "$work/sort" -o "$work/hmq.bam" "$work/selected.bam"
            samtools index -@ {params.workers} "$work/hmq.bam"
            mv "$work/redup.bam" {output.redup_bam:q}
            mv "$work/redup.bam.bai" {output.redup_idx:q}
            mv "$work/hmq.bam" {output.hmq_bam:q}
            mv "$work/hmq.bam.bai" {output.hmq_idx:q}
            rm "$work/selected.bam"
            rmdir "$work"
            """
    for _suffix in (".redup.bam", ".redup.bam.bai", ".redup.HMQ.sorted.bam", ".redup.HMQ.sorted.bam.bai"):
        optional_output.append(expand(output_dir + "/bam2bam/{subj}.triocompass.addRG" + _suffix, subj=subjs))
