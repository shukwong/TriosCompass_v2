import json

_rufus_config = config.get("rufus", {})
_rufus_work = _rufus_config.get("work_dir", output_dir + "/rufus/work")
_rufus_image = tc_container("rufus")
if _rufus_image and _rufus_config.get("container"):
    raise ValueError("Use containers.rufus for the entire rule or rufus.container for the caller, not both")


rule run_rufus:
    input:
        bams=get_bams_by_family,
        ped=config["ped_dir"] + "/{fam}.ped",
        ref=genome,
        fai=genome_fai,
        bwa=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        helper=workflow.source_path("../scripts/run_rufus.py")
    output:
        vcf=report(
            output_dir + "/rufus/{fam}.rufus.vcf",
            caption="../report/dnSV_rufus.rst",
            category="RUFUS variants",
            subcategory="Calls",
            labels={"Family": "{fam}", "File type": "VCF"}
        )
    params:
        samples=lambda w: [person.id for person in families[w.fam]],
        settings=json.dumps(_rufus_config),
        work_dir=_rufus_work
    threads:
        _rufus_config.get("threads", config.get("threads", {}).get("rufus", 24))
    resources:
        mem_mb=_rufus_config.get("mem_mb", 256000),
        runtime=_rufus_config.get("runtime_minutes", 5760)
    log:
        output_dir + "/logs/rufus/{fam}.log"
    benchmark:
        output_dir + "/benchmark/rufus/{fam}.tsv"
    container:
        _rufus_image
    shell:
        """
        python {input.helper:q} run \
            --ped {input.ped:q} --family {wildcards.fam:q} \
            --samples {params.samples:q} --bams {input.bams:q} \
            --reference {input.ref:q} --threads {threads} \
            --settings-json {params.settings:q} \
            --work-dir {params.work_dir:q} --output {output.vcf:q} > {log:q} 2>&1
        """


rule rufus_summary:
    input:
        vcfs=expand(output_dir + "/rufus/{fam}.rufus.vcf", fam=fam_ids),
        helper=workflow.source_path("../scripts/run_rufus.py")
    output:
        report(
            output_dir + "/rufus/rufus_summary.tsv",
            caption="../report/dnSV_rufus.rst",
            category="RUFUS variants",
            subcategory="Record counts"
        )
    params:
        families=fam_ids
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30
    container:
        tc_container("utility")
    shell:
        """
        python {input.helper:q} summary \
            --families {params.families:q} --vcfs {input.vcfs:q} --output {output:q}
        """


optional_output.append(expand(output_dir + "/rufus/{fam}.rufus.vcf", fam=fam_ids))
optional_output.append(output_dir + "/rufus/rufus_summary.tsv")
