rule mark_duplicates_spark:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        bam="VC_recal/dedup_{sample}.sorted.bam",
        metrics="VC_recal/{sample}.metrics.txt"
    log:
        "logs/gatk_prepro/mark_dup_{sample}.log"
    params:
        extra="--conf 'spark.executor.cores=4'",  # optional
        java_opts="-Xmx4g"
    resources:
        mem_mb=20000
    threads: 4
    shell:
    	"gatk --java-options '{params.java_opts}' MarkDuplicatesSpark "
    	"{params.extra} "
    	"-I {input} "
    	"-O {output.bam} "
    	"-M {output.metrics} &> {log}"

rule samtools_index:
    input:
        "VC_recal/dedup_{sample}.sorted.bam"
    output:
        "VC_recal/dedup_{sample}.sorted.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"

rule gatk_baserecalibratorspark:
    input:
        bam="VC_recal/dedup_{sample}.sorted.bam",
        ref=config["gatk"]["fasta"],
        dict=config["gatk"]["dict"],
        known=config["gatk"]["known"]  # optional known sites
    output:
        recal_table="VC_recal/{sample}.grp"
    log:
        "logs/gatk_prepro/baserecalibrator_{sample}.log"
    params:
        extra="--conf 'spark.executor.cores=4'",  # optional
        java_opts="-Xmx4g"
    resources:
        mem_mb=20000
    threads: 4
    shell:
    	"gatk --java-options '{params.java_opts}' BaseRecalibratorSpark {params.extra} "
    	"-R {input.ref} -I {input.bam} "
    	"--output {output.recal_table} --known-sites {input.known} &> {log}"

rule gatk_applybqsr_spark:
    input:
        bam="VC_recal/dedup_{sample}.sorted.bam",
        ref=config["gatk"]["fasta"],
        dict=config["gatk"]["dict"],
        recal_table="VC_recal/{sample}.grp"
    output:
        bam="preprocessed_VC/{sample}.bam"
    log:
        "logs/gatk_prepro/gatk_applybqsr_spark_{sample}.log"
    params:
        extra="--conf 'spark.executor.cores=4'",  # optional
        java_opts="-Xmx4g"
    resources:
        mem_mb=20000
    threads: 4
    shell:
    	"gatk --java-options '{params.java_opts}' ApplyBQSRSpark {params.extra} "
	"--reference {input.ref} --input {input.bam} "
	"--bqsr-recal-file {input.recal_table} "
	"--output {output.bam} &> {log}"
