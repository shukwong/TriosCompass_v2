import peppy
from itertools import compress
import os
import peds

# Process fastq sequences using fastp, fastqc, and fastq_screen

pepfile: os.path.join(workflow.basedir, "..", config["pepfile"])

### Disable pepschema, due to:
# Could not load remote file: http://schema.databio.org/pep/2.0.0.yaml. Original exception: URLError(TimeoutError(110, 'Connection timed out'))
# pepschema: config["pepschema"]

samples = pep.sample_table
subjs = list(set(samples["SAMPLE_ID"].tolist())) # 
sample_bam_dict=dict(zip(samples.SAMPLE_ID, samples.BAM))

# print (sample_bam_dict)

if config["bam_input"]["reset_RG"]:
    if TC_CONTAINER_ONLY or (workflow.use_singularity and not workflow.use_conda):
        rule replace_rg:
            input: lambda w: sample_bam_dict[w.subj]
            output: output_dir + "/fixed-rg/{subj}.bam"
            threads: config["threads"]["replace_rg"]
            resources: mem_mb=16000
            params: heap_mb=lambda w, resources: max(1, int(resources.mem_mb * 0.8))
            benchmark: output_dir + "/benchmark/replace_rg/{subj}.tsv"
            container: tc_container("picard", "docker://quay.io/biocontainers/picard:2.27.3--hdfd78af_0")
            shell:
                """
                unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
                picard -Xms1m -Xmx{params.heap_mb}m -XX:+UseSerialGC \
                    AddOrReplaceReadGroups I={input:q} O={output:q} \
                    RGLB=lib1 RGPL=illumina RGPU={wildcards.subj:q} \
                    RGSM={wildcards.subj:q} CREATE_INDEX=true
                """
    else:
        rule replace_rg:
            input: lambda w: sample_bam_dict[w.subj]
            output:
                output_dir + "/fixed-rg/{subj}.bam"
            params:
                extra="--RGLB lib1 --RGPL illumina --RGPU {subj} --RGSM {subj} --CREATE_INDEX true "
            benchmark:
                output_dir +"/benchmark/replace_rg/{subj}.tsv"
            threads: config["threads"]["replace_rg"]
            wrapper:
                "v1.25.0/bio/picard/addorreplacereadgroups"

    def get_bam(subj):
        return output_dir +"/fixed-rg/{}.bam".format(subj) 


else:
    def get_bam(subj):
        return sample_bam_dict[subj] 

    # print(get_bam("SC501095"))

