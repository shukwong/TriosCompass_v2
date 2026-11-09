# FASTQ alignment and RUFUS consume BWA indices. Other BAM-only runs can
# request reference preparation explicitly by configuring containers.bwa.
if (config["fastq_input"]["enable"] or config.get("rufus", {}).get("enable", False)
        or config.get("containers", {}).get("bwa") or not TC_CONTAINER_ONLY):
    if TC_CONTAINER_ONLY or config.get("containers", {}).get("bwa"):
        rule bwa_index:
            input: genome
            output: idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa")
            threads: config["threads"]["bwa_index"]
            container: tc_container("bwa")
            shell: "bwa index {input:q}"
    else:
        rule bwa_index:
            input: genome
            output: idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa")
            threads: config["threads"]["bwa_index"]
            wrapper: "v2.3.2/bio/bwa/index"

rule genome_dict:
    input: genome
    output: genome_dict
    container: tc_container("samtools", "docker://euformatics/samtools:1.19.2")
    shell: "samtools dict {input:q} > {output:q}"

rule genome_faidx:
    input: genome
    output: genome_fai
    container: tc_container("samtools", "docker://euformatics/samtools:1.19.2")
    shell: "samtools faidx {input:q}"
