def get_fastq(wildcards):
    return samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

rule bwa_mem:
    input:
        reads = get_fastq
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=config["ref"]["index"],
        extra=r"-R '@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tSM:{sample}'",
        sorting="none"
    threads: 4
    resources: cpus=4, mem_mb=20000
    wrapper:
        "0.77.0/bio/bwa/mem"

rule sambamba_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    params:
        "-m 4G"  # optional parameters
    log:
        "logs/sambamba-sort/{sample}.log"
    threads: 4
    resources: cpus=4, mem_mb=20000
    wrapper:
        "0.77.0/bio/sambamba/sort"
