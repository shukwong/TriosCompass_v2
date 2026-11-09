import uuid

rule deepvariant_pb:
    input:
        bam=get_bam_by_subj,
        ref=genome
    output:
        gvcf=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf",
        gz=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz",
        tbi=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz.tbi"
    threads: config["threads"]["deepvariant_pb"]
    benchmark:
        output_dir +"/benchmark/deepvariant_pb/{subj}.tsv"
    conda: "../envs/tabix.yaml"
    params:
        singularity_cmd = tc_parabricks(),
        bgzip_cmd = tc_exec("tabix", "bgzip"),
        tabix_cmd = tc_exec("tabix", "tabix")
    shell: """
        {params.singularity_cmd} \
        pbrun deepvariant \
        --gvcf  \
        --ref {input.ref} \
        --in-bam {input.bam} \
        --out-variants {output.gvcf} 
        {params.bgzip_cmd} -c {output.gvcf} > {output.gz}
        {params.tabix_cmd} {output.gz}
    """

rule glnexus_dv: 
    input:
        gvcf=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz", id=[person.id for person in families[w.fam]]),
        tbi=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz.tbi", id=[person.id for person in families[w.fam]]),
        ref=genome
    output:
        vcf = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz",
        tbi = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz.tbi"
    threads: config["threads"]["glnexus_dv"]
    benchmark:
        output_dir +"/benchmark/glnexus_dv/dv_{fam}.tsv"
    conda: "../envs/glnexus.yaml"
    container: tc_container("glnexus")
    params: 
        tempdir=temp(directory(output_dir +"/glnexus/GLnexus.DB_{}")).format(uuid.uuid4()),
        config=config["deepvariant"]["glnexus_config"] 
    shell: """
        glnexus_cli -t {threads} --config {params.config} --dir {params.tempdir} {input.gvcf} |  bcftools norm -f {input.ref} -m - -O z -o {output.vcf} 
        tabix -p vcf {output.vcf}
        rm -fr {params.tempdir}
    """