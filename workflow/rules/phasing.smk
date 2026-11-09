
PHASE_WIN = config["phasing"]["window_size"]

#############################################
# rules for phase DNMs
#############################################
checkpoint scatter_dnms:
    container: tc_container("utility")
    input: output_dir +"/dnm_vcf/{fam}.dnm.vcf.gz"
    output:
        directory(output_dir + "/phase_DNMs/{fam}/variants")
    shell: """
        mkdir -p {output:q}
        zcat {input:q} | awk -v FS='\\t' -v OFS=':' -v outdir={output:q} '
            !/^#/ && NF {{
                out = outdir "/" $1 "_" $2 ".dnm"
                if (out in seen) {{
                    print $1,$2,$4,$5 >>out
                }} else {{
                    print $1,$2,$4,$5 >out
                    seen[out] = 1
                }}
                if (close(out) != 0) exit 1
            }}'
    """

### Extract region around each DNM from the VCF file
rule vcf4phase:
    input:
        dnm = output_dir + "/phase_DNMs/{fam}/variants/{chr}_{pos}.dnm",
        vcf = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz",
        tbi = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz.tbi"
    output:
        vcf=output_dir + "/phase_DNMs/{fam}/variants/DV.{chr}_{pos}.vcf.gz",
        tbi=output_dir + "/phase_DNMs/{fam}/variants/DV.{chr}_{pos}.vcf.gz.tbi"
    conda: "../envs/bcftools.yaml"
    container: tc_container("bcftools")
    benchmark:
        output_dir +"/benchmark/vcf4phase/{fam}.{chr}_{pos}.tsv"
    params: 
        reg= lambda w: '%s:%d-%d' % (w.chr, int(w.pos)-PHASE_WIN, int(w.pos)+PHASE_WIN)
    shell: """
        bcftools view -r {params.reg} {input.vcf} | bcftools annotate -Oz -x ID -I +"%CHROM:%POS:%REF:%ALT" -Oz -o {output.vcf};
        tabix {output.vcf};
    """

rule phase_child:
    input: 
        bams=get_child_bam_by_family,
        ref=genome,
        vcf = output_dir + "/phase_DNMs/{fam}/variants/DV.{chr}_{pos}.vcf.gz"
    output: 
        vcf= output_dir +"/phase_DNMs/{fam}/variants/{fam}.{chr}_{pos}.phase_child.vcf.gz"
    params:
        id = lambda w: CHILD_DICT[w.fam]
    benchmark:
        output_dir +"/benchmark/phase_child/{fam}.{chr}_{pos}.tsv"
    threads: config["threads"]["phase_child"]
    conda: "../envs/whatshap.yaml"
    container: tc_container("whatshap")
    shell: """
        whatshap phase -o {output.vcf} --tag=PS --indels --reference={input.ref} --sample {params.id} {input.vcf} {input.bams}
    """

rule phase_trios:
    input: 
        bams=get_bams_by_family,
        ped=ped_dir + "/{fam}.ped",
        ref=genome,
        vcf = output_dir + "/phase_DNMs/{fam}/variants/DV.{chr}_{pos}.vcf.gz"
    output: 
        vcf= output_dir +"/phase_DNMs/{fam}/variants/{fam}.{chr}_{pos}.phase_trios.vcf.gz"
    benchmark:
        output_dir +"/benchmark/phase_trios/{fam}.{chr}_{pos}.tsv"
    threads: config["threads"]["phase_trios"]
    conda: "../envs/whatshap.yaml"
    container: tc_container("whatshap")
    shell: """
        whatshap phase -o {output.vcf} --tag=PS --indels --ped {input.ped} --reference={input.ref} {input.vcf} {input.bams} 
    """


rule extract_parental_origin:
    container: tc_container("utility")
    input:
        dnm = output_dir + "/phase_DNMs/{fam}/variants/{chr}_{pos}.dnm",
        ped = ped_dir + "/{fam}.ped",
        phase_child = output_dir +"/phase_DNMs/{fam}/variants/{fam}.{chr}_{pos}.phase_child.vcf.gz",
        phase_trios = output_dir +"/phase_DNMs/{fam}/variants/{fam}.{chr}_{pos}.phase_trios.vcf.gz"
    output:
        output_dir + "/phase_DNMs/{fam}/variants/{chr}_{pos}.parental_origin.txt"
    benchmark:
        output_dir +"/benchmark/extract_parental_origin/{fam}.{chr}_{pos}.tsv"
    params: 
        perl_cmd = config["phasing"]["perl_cmd"]
    shell: """
        {params.perl_cmd}  {input.dnm} {input.ped} {input.phase_child} {input.phase_trios} {output}
    """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.scatter_dnms.get(**wildcards).output[0]
    return expand(output_dir + "/phase_DNMs/{fam}/variants/{chr_pos}.parental_origin.txt",
           fam=wildcards.fam,
           chr_pos=glob_wildcards(os.path.join(checkpoint_output, "{chr_pos}.dnm")).chr_pos)


### should be replace "_" with ":" after VID is fixed
rule aggregate_phase:
    container: tc_container("utility")
    input: aggregate_input
    output: 
        report(
            output_dir + "/phase_DNMs/{fam}.parental_origin.tab",
            caption="../report/parental_origin.rst",
            category="De novo mutations",
            subcategory="Phasing",
            labels={
                "Family": "{fam}",
                "Desc": "parental origin of DNMs",
                "File type": "tab/tsv"
            }
        )
    shell: """
        (
            printf 'VariantID\\tParentalOrigin\\tPhase\\tHaplotypeBlockSize\\t InformativeSiteCnt\\tFM_count\\tMF_count\\tFM_status\\tParentOriginChange\\n'
            set -- {input:q}
            if [ "$#" -gt 0 ]; then
                cat -- "$@" | sort -t ":" -Vs -k1,1 -k2,2n
            fi
        ) > {output:q}
    """

optional_output.append(expand(output_dir + "/phase_DNMs/{fam}.parental_origin.tab", fam=fam_ids))