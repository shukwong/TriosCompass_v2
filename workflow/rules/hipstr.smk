## split bed files first
# skip chrX for the time being
rule split_bed_hipstr: 
    input: 
        config["dnSTR"]["hipstr"]["ref_panel"]
    output:
        expand(output_dir+"/splitted_panel/hipstr_{chunk}.bed", chunk=CHUNKS)
    params: 
        split_total = config["dnSTR"]["split_n"],
        prefix=output_dir+"/splitted_panel/hipstr_"
    shell: '''
        # split cannot take piped input, leading to the error of "cannot determine file size"
        # so we use tmpfile here to get around the issue
        tmpfile=$(mktemp /tmp/abc-script.XXXXXX); grep -v -e "^chrX" -e "chrY" {input} > $tmpfile || true ; split --numeric-suffixes=0 -n l/{params.split_total} --suffix-length=5  --additional-suffix=".bed" $tmpfile {params.prefix}
    '''

rule hipstr:
    input: 
        bams = [get_bam(subj) for subj in final_subjs],
        ref=genome,
        reg = output_dir+"/splitted_panel/hipstr_{chunk}.bed"
    output:
        vcf= output_dir + "/hipstr/{chunk}.vcf.gz",
        vcf_index= output_dir + "/hipstr/{chunk}.vcf.gz.tbi",
        viz= output_dir + "/hipstr/{chunk}.aln.viz.gz",
        viz_index= output_dir + "/hipstr/{chunk}.aln.viz.gz.tbi",
    log: output_dir + "/hipstr/{chunk}.log"
    benchmark:
        output_dir + "/benchmark/hipstr/{chunk}.tsv"
    singularity: config.get("containers", {}).get("hipstr", "docker://cgrlab/hipstr:latest")
    params: 
        bams=lambda w, input: ",".join(input.bams)
    shell: """
        output_parent=$(dirname {output.vcf:q})
        log_path={log:q}
        mkdir -p "$output_parent" "$(dirname "$log_path")"
        tmpdir=$(mktemp -d "$output_parent/.tmp.XXXXXX")
        hipstr_cleanup() {{
            rc=$?
            if [ -f "$tmpdir/hipstr.log" ]; then
                cp -- "$tmpdir/hipstr.log" "$log_path" || true
            fi
            rm -rf -- "$tmpdir" || true
            exit "$rc"
        }}
        trap hipstr_cleanup EXIT
        trap 'exit 130' INT
        trap 'exit 143' TERM

        HipSTR --output-gls --bams {params.bams:q} --fasta {input.ref:q} \
            --regions {input.reg:q} --str-vcf "$tmpdir/calls.vcf.gz" \
            --viz-out "$tmpdir/alignments.aln.viz.gz" --def-stutter-model \
            --log "$tmpdir/hipstr.log"
        tabix -p bed "$tmpdir/alignments.aln.viz.gz"
        tabix -p vcf "$tmpdir/calls.vcf.gz"
        mv -- "$tmpdir/calls.vcf.gz" {output.vcf:q}
        mv -- "$tmpdir/calls.vcf.gz.tbi" {output.vcf_index:q}
        mv -- "$tmpdir/alignments.aln.viz.gz" {output.viz:q}
        mv -- "$tmpdir/alignments.aln.viz.gz.tbi" {output.viz_index:q}
    """

########################################################
# generate VizAln plots for every hipstr prediction
########################################################

### extract bed file from the reference panel by chr, pos from MonSTR output
# skip output_dir + "/vizaln/{fam}/{fam}_recall.both.bed"

rule bedfile_for_hipstr_call:
    input: 
        tab=output_dir +"/monstr_filter/hipstr.filtered.tab",
        ref_panel = config["dnSTR"]["hipstr"]["ref_panel"]
    output: 
        tmp = temp(output_dir + "/vizaln/{fam}/{fam}_recall.hipstr.tmp"),
        bed = output_dir + "/vizaln/{fam}/{fam}_recall.hipstr.bed"
    params:
        cmd = "sed 's/^/chr/'",
        id = lambda w: CHILD_DICT[w.fam]
    shell: """
        awk -v FS='\\t' -v OFS='\\t' '{{if(NR>1 && $6=="{params.id}") print $1, $2}}' {input.tab} | {params.cmd} > {output.tmp}
        if [ -s {output.tmp} ]; then 
            grep -f {output.tmp} {input.ref_panel} > {output.bed} || true;
        else
            echo -n > {output.bed}
        fi
    """



rule hipstr_recall:
    input: 
        bams = get_bams_by_family,
        ref=genome,
        reg = output_dir + "/vizaln/{fam}/{fam}_recall.hipstr.bed"
    output:
        vcf = output_dir + "/vizaln/{fam}/{fam}.hipstr.vcf.gz",
        vcf_index = output_dir + "/vizaln/{fam}/{fam}.hipstr.vcf.gz.tbi",
        viz = output_dir + "/vizaln/{fam}/{fam}.hipstr.aln.viz.gz",
        viz_index = output_dir + "/vizaln/{fam}/{fam}.hipstr.aln.viz.gz.tbi"
    log: output_dir + "/vizaln/{fam}/{fam}.hipstr.hipstr.log"
    params:
        bams=lambda w, input: ",".join(input.bams)
    singularity: config.get("containers", {}).get("hipstr", "docker://cgrlab/hipstr:latest")
    benchmark:
        output_dir +"/benchmark/hipstr_recall/{fam}_hipstr.tsv"
    shell: """
        output_parent=$(dirname {output.vcf:q})
        log_path={log:q}
        mkdir -p "$output_parent" "$(dirname "$log_path")"
        tmpdir=$(mktemp -d "$output_parent/.tmp.XXXXXX")
        hipstr_cleanup() {{
            rc=$?
            if [ -f "$tmpdir/hipstr.log" ]; then
                cp -- "$tmpdir/hipstr.log" "$log_path" || true
            fi
            rm -rf -- "$tmpdir" || true
            exit "$rc"
        }}
        trap hipstr_cleanup EXIT
        trap 'exit 130' INT
        trap 'exit 143' TERM

        HipSTR --output-gls --bams {params.bams:q} --fasta {input.ref:q} \
            --regions {input.reg:q} --str-vcf "$tmpdir/calls.vcf.gz" \
            --viz-out "$tmpdir/alignments.aln.viz.gz" --def-stutter-model --min-reads 6 \
            --log "$tmpdir/hipstr.log"
        tabix -p bed "$tmpdir/alignments.aln.viz.gz"
        tabix -p vcf "$tmpdir/calls.vcf.gz"
        mv -- "$tmpdir/calls.vcf.gz" {output.vcf:q}
        mv -- "$tmpdir/calls.vcf.gz.tbi" {output.vcf_index:q}
        mv -- "$tmpdir/alignments.aln.viz.gz" {output.viz:q}
        mv -- "$tmpdir/alignments.aln.viz.gz.tbi" {output.viz_index:q}
    """

checkpoint scatter_chr_pos:
    input: 
        output_dir + "/vizaln/{fam}/{fam}.hipstr.vcf.gz"
    output:
        dir = directory(output_dir + "/vizaln/{fam}/hipstr/variants/"),
        tmp = temp(output_dir + "/vizaln/{fam}/hipstr/variants/tmp")
    shell: """
        mkdir -p {output.dir}
        zgrep -v "^#" {input}  > {output.tmp} || true
        if [ -s {output.tmp} ]; then 
            awk -v FS='\\t' 'system("touch {output.dir}/"$1"_"$2".dnm")' {output.tmp}
        else
            # touch "$(dirname {output.dir})/DONE"
            echo "Do noting now"
        fi
        
    """

rule vizaln:
    input: 
        dummy=output_dir + "/vizaln/{fam}/hipstr/variants/{chr}_{pos}.dnm",
        viz=output_dir + "/vizaln/{fam}/{fam}.hipstr.aln.viz.gz"
    output: 
        report(
            output_dir + "/vizaln/{fam}/hipstr/variants/{chr}_{pos}.html",
            caption="../report/visaln.rst",
            category="Visulization of dnSTR",
            subcategory="{fam}",
            labels={
                "chr": "{chr}",
                "pos": "{pos}",
                "Desc": "VizAln of dnSTR",
                "File type": "HTML"
            }
        )
    benchmark:
        output_dir +"/benchmark/vizaln/{fam}_hipstr.{chr}_{pos}.tsv"
    singularity: config.get("containers", {}).get("hipstr", "docker://cgrlab/hipstr:latest")
    shell: """
        VizAln_rev {input.viz} {output} {wildcards.chr} {wildcards.pos}
    """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.scatter_chr_pos.get(**wildcards).output[0]
    return expand(output_dir + "/vizaln/{fam}/hipstr/variants/{chr_pos}.html",
           fam=wildcards.fam,
           chr_pos=glob_wildcards(os.path.join(checkpoint_output, "{chr_pos}.dnm")).chr_pos)


rule aggregate_visaln:
    input: aggregate_input
    output: 
        output_dir + "/vizaln/{fam}/hipstr/DONE"
    shell: """
        touch {output}
    """


optional_output.append( expand(output_dir + "/vizaln/{fam}/hipstr/DONE", fam=fam_ids))