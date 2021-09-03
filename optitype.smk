def get_fastq1(wildcards):
    return samples.loc[wildcards.sample, ["fq1"]].dropna()
def get_fastq2(wildcards):
    return samples.loc[wildcards.sample, ["fq2"]].dropna()

ids = [1,2]    
rule razers3:
    input:
        # list of input reads
        reads= get_fastq1
    output:
        # output format is automatically inferred from file extension. Can be bam/sam or other formats.
        "mapped/{sample}.fished_1.bam"
    log:
        "logs/razers3/{sample}.fished_1.log"
    params:
        # the reference genome
        genome=config["ref"]["hla"],
        # additional parameters
        extra="-i 95 -m 1 -dr 0"
    threads: 4
    resources: cpus=4, mem_mb=20000
    wrapper:
        "0.77.0/bio/razers3"

rule razers3_2:
    input:
        # list of input reads
        reads= get_fastq2
    output:
        # output format is automatically inferred from file extension. Can be bam/sam or other formats.
        "mapped/{sample}.fished_2.bam"
    log:
        "logs/razers3/{sample}.fished_2.log"
    params:
        # the reference genome
        genome=config["ref"]["hla"],
        # additional parameters
        extra="-i 95 -m 1 -dr 0"
    threads: 4
    resources: cpus=4, mem_mb=20000
    wrapper:
        "0.77.0/bio/razers3"

rule samtools_bam2fq_interleaved_1:
    input:
        "mapped/{sample}.fished_1.bam"
    output:
        "preprocessed_optitype/{sample}.fished_1.fastq"
    params:
        " "
    threads: 2
    wrapper:
        "0.77.0/bio/samtools/bam2fq/interleaved"

rule samtools_bam2fq_interleaved_2:
    input:
        "mapped/{sample}.fished_2.bam"
    output:
        "preprocessed_optitype/{sample}.fished_2.fastq"
    params:
        " "
    threads: 2
    wrapper:
        "0.77.0/bio/samtools/bam2fq/interleaved"

rule optitype:
    input:
        # list of input reads
        reads=["preprocessed_optitype/{sample}.fished_1.fastq", "preprocessed_optitype/{sample}.fished_2.fastq"]
    output:
        multiext("optitype_result/{sample}", "_coverage_plot.pdf", "_result.tsv")
    log:
        "logs/optitype/{sample}.log"
    params:
        # Type of sequencing data. Can be 'dna' or 'rna'. Default is 'dna'.
        sequencing_type="dna",
        # additional parameters
        extra="-v"
    threads: 4
    resources: cpus=4
    wrapper:
        "0.77.0/bio/optitype"
