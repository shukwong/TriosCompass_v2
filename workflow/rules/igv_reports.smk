# Opt-in standalone SV inspection; pedigree and BAM headers establish identities.
_igv_config = config.get("igv_reports", {})
_igv_vcf = _igv_config.get("vcf_template")
_igv_source = _igv_config.get("source", "joint")
# Native joint columns are caller summaries, even if one looks like a sample ID.
_igv_trio_samples = _igv_config.get("vcf_samples_are_trio", _igv_source != "joint")
if not _igv_vcf and _igv_source == "joint":
    _igv_trio_samples = False
if not _igv_vcf:
    _igv_sources = {
        "joint": output_dir + "/joint_dnSV/{fam}.dnSV.vcf",
        "manta": output_dir + "/dnSV_manta/{fam}.dnSV.manta.vcf",
        "rufus": output_dir + "/rufus/{fam}.rufus.vcf",
    }
    if _igv_source not in _igv_sources:
        raise ValueError("igv_reports.source must be joint, manta or rufus, or supply vcf_template")
    _igv_vcf = _igv_sources[_igv_source]
_igv_output = _igv_config.get("output_dir", output_dir + "/igv_reports")
_igv_js = _igv_config.get("igv_js")

rule igv_reports_dnsv:
    input:
        vcf=_igv_vcf,
        ped=config["ped_dir"] + "/{fam}.ped",
        bams=get_bams_by_family,
        ref=genome,
        fai=genome_fai,
        igv_js=[_igv_js] if _igv_js else [],
        helper=workflow.source_path("../scripts/igv_report_sv.py"),
        controls=workflow.source_path("../scripts/igv_read_controls.html")
    output:
        gallery=report(
            directory(_igv_output + "/{fam}"),
            patterns=["{gallery,index}.html"],
            caption="../report/igv_reports.rst",
            category="De novo SVs",
            subcategory="IGV reports",
            labels={"Family": "{fam}", "File type": "Interactive HTML"}
        )
    params:
        samples=lambda w: [person.id for person in families[w.fam]],
        js_args=["--igv-js", _igv_js] if _igv_js else [],
        sample_args=[] if _igv_trio_samples else ["--ignore-vcf-samples"],
        flank=_igv_config.get("flank", 2000),
        initial_window=_igv_config.get("initial_window", 400),
        track_height=_igv_config.get("track_height", 320),
        read_height=_igv_config.get("read_height", 14),
        coverage_max=_igv_config.get("coverage_max", 200),
        sampling_depth=_igv_config.get("sampling_depth", 500),
        min_sv_length=_igv_config.get("min_sv_length", 50),
        exclude_flags=_igv_config.get("exclude_flags", 1536),
        max_variants=_igv_config.get("max_variants", 0)
    threads: 1
    resources:
        mem_mb=_igv_config.get("mem_mb", 4000),
        runtime=_igv_config.get("runtime_minutes", 60)
    log:
        output_dir + "/logs/igv_reports/{fam}.log"
    benchmark:
        output_dir + "/benchmark/igv_reports/{fam}.tsv"
    conda:
        "../envs/igv-reports.yaml"
    container:
        tc_container("igv-reports")
    shell:
        """
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
        python {input.helper:q} \
            --vcf {input.vcf:q} --ped {input.ped:q} --family {wildcards.fam:q} \
            --samples {params.samples:q} --bams {input.bams:q} \
            --reference {input.ref:q} --output-dir {output.gallery:q} --controls-template {input.controls:q} \
            {params.js_args:q} {params.sample_args:q} \
            --flank {params.flank:q} --initial-window {params.initial_window:q} \
            --track-height {params.track_height:q} --read-height {params.read_height:q} \
            --coverage-max {params.coverage_max:q} --sampling-depth {params.sampling_depth:q} \
            --min-sv-length {params.min_sv_length:q} --exclude-flags {params.exclude_flags:q} \
            --max-variants {params.max_variants:q} > {log:q} 2>&1
        """

optional_output.append(expand(_igv_output + "/{fam}", fam=fam_ids))
