def _gatk_java_heap_mb(wildcards, resources):
    mem_mb = int(getattr(resources, "mem_mb", 16000))
    if mem_mb < 512:
        raise ValueError("gatk_combine_gvcf requires at least 512 MB of memory")
    # Reserve memory for the JVM and native libraries, including small jobs.
    reserve_mb = min(4096, max(256, mem_mb // 5))
    return mem_mb - reserve_mb


rule gatkhc_pb:
    input:
        bam=get_bam_by_subj,
        ref=genome,
        ref_dict=genome_dict
    output:
        gvcf=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf",
        gz=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf.gz",
        tbi=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf.gz.tbi"
    threads: config["threads"]["gatkhc_pb"]
    resources:
        gpu=1
    conda: "../envs/tabix.yaml"
    params: 
        singularity_cmd = tc_parabricks(),
        bgzip_cmd = tc_exec("tabix", "bgzip"),
        tabix_cmd = tc_exec("tabix", "tabix")
    benchmark:
        output_dir +"/benchmark/gatkhc_pb/{subj}.tsv"
    shell: """
        {params.singularity_cmd} \
        pbrun haplotypecaller \
            --gvcf \
            --ref {input.ref} \
            --in-bam {input.bam} \
            --out-variants {output.gvcf}   
        {params.bgzip_cmd} -c {output.gvcf} > {output.gz}
        {params.tabix_cmd} {output.gz}
    """

rule gatk_combine_gvcf:
    input: 
        vcfs=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.g.vcf.gz", id= [person.id for person in families[w.fam]]),
        tbis=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.g.vcf.gz.tbi", id= [person.id for person in families[w.fam]]),
        ref=genome,
        fai=genome_fai,
    output: output_dir+"/gatk_combine_gvcf/{fam}.combine_gvcf.g.vcf.gz"
    threads: config["threads"]["gatk_combine_gvcf"]
    resources:
        mem_mb=16000
    conda: "../envs/gatk4.yaml"
    container: tc_container("gatk4")
    benchmark:
        output_dir + "/benchmark/gatk_combine_gvcf/{fam}.tsv"
    params: 
        v=lambda w, input: " -V ".join(input.vcfs),
        java_heap_mb=_gatk_java_heap_mb
    shell: """
        unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
        : "${{CONDA_PREFIX:?The GATK Conda environment must be activated}}"
        export PATH="$CONDA_PREFIX/bin:$PATH"
        export JAVA_HOME="$CONDA_PREFIX"
        hash -r
        "$CONDA_PREFIX/bin/gatk" --java-options "-Xmx{params.java_heap_mb}m" CombineGVCFs \
            -V {params.v} \
            -R {input.ref} \
            -O {output} 
    """

rule gatk_genotype_gvcf_pb:
    input: 
        gvcf=output_dir+"/gatk_combine_gvcf/{fam}.combine_gvcf.g.vcf.gz",
        ref=genome
    output:
        vcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf",
        gz=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
    threads: config["threads"]["gatk_genotype_gvcf_pb"]
    resources:
        gpu=1
    conda: "../envs/tabix.yaml"
    benchmark:
        output_dir + "/benchmark/gatk_genotype_gvcf_pb/{fam}.tsv"
    params:
        singularity_cmd = tc_parabricks(),
        bgzip_cmd = tc_exec("tabix", "bgzip"),
        tabix_cmd = tc_exec("tabix", "tabix")
    shell: """
        {params.singularity_cmd} \
        pbrun genotypegvcf \
            --in-gvcf {input.gvcf} \
            --ref {input.ref} \
            --out-vcf {output.vcf} 
        {params.bgzip_cmd} -@ {threads} -c {output.vcf}> {output.gz} && \
        {params.tabix_cmd} {output.gz}
    """

# drop CalculateGenotypePosteriors completely
# gatk CalculateGenotypePosteriors \
        #     -V {input.gvcf} \
        #     -O {output.gvcf} 
        #     -ped {input.ped} --skip-population-priors


if config["gatk_hc"]["gatk_hard_filter"]["enable"]:
    rule gatk_cgp:
        input: 
            gvcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
            # ped=config["ped_dir"]+"/{fam}.ped",
            ref=genome
        output: 
            # gvcf=output_dir+"/gatk_cgp/{fam}.cgp.g.vcf.gz",
            vcf=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz",
            tbi=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz.tbi",
        threads: config["threads"]["gatk_cgp"]
        conda: "../envs/bcftools.yaml"
        container: tc_container("bcftools")
        benchmark:
            output_dir + "/benchmark/gatk_cgp/{fam}.tsv"
        params: filter=config["gatk_hc"]["gatk_hard_filter"]["filter"]
        resources: 
        shell: """
            bcftools norm -f {input.ref} -m -  {input.gvcf} | bcftools view -i'ALT!="*"' | bcftools filter -i '{params.filter}' -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}
        """

else: 
    rule gatk_cgp:
        input: 
            gvcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
            ped=config["ped_dir"]+"/{fam}.ped",
            ref=genome
        output: 
            vcf=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz",
            tbi=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz.tbi",
        threads: config["threads"]["gatk_cgp"]
        conda: "../envs/bcftools.yaml"
        container: tc_container("bcftools")
        benchmark:
            output_dir + "/benchmark/gatk_cgp/{fam}.tsv"
        resources: 
        shell: """
            bcftools norm -f {input.ref} -m -  {input.gvcf} | bcftools view -i'ALT!="*"' -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}
        """
