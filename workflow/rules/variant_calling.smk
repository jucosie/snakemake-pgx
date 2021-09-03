rule haplotype_caller:
    input:
        # single or list of bam files
        bam="preprocessed_VC/{sample}.bam",
        ref=config["gatk"]["fasta"]
    output:
        gvcf="preprocessed_VC/{sample}.g.vcf.gz"
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="-L {} --output-mode EMIT_ALL_ACTIVE_SITES".format(config["gatk"]["intervals"]),  # optional
    threads: 4
    resources: cpus=4
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller {params.extra} "
        "-R {input.ref} -I {input.bam} "
        "-ERC BP_RESOLUTION -O {output.gvcf} > {log}"
        
rule variant_annotator:
    input:
    	gvcf="preprocessed_VC/{sample}.g.vcf.gz",
    	dbsnp=config["gatk"]["dbsnp"]
    output:
    	anngvcf="VC_interm_dir/{sample}.annotated.g.vcf"
    log:
        "logs/gatk/variantannotator/{sample}.log"
    params:
        extra="-A RMSMappingQuality -A ReadPosRankSumTest -A QualByDepth -A StrandOddsRatio"  # optional
    threads: 4
    resources: cpus=4
    shell:
    	"gatk VariantAnnotator -V {input.gvcf} " 
    	"-O {output.anngvcf} --dbsnp {input.dbsnp} {params.extra} > {log}"


rule genomics_db_import:
    input:
        gvcfs=expand("VC_interm_dir/{sample}.annotated.g.vcf", sample=samples.index)
    output:
        db=directory("db"),
    log:
        "logs/gatk/genomicsdbimport.log"
    params:
        intervals=config["gatk"]["intervals"],
        db_action="--genomicsdb-workspace-path", # optional
        gvcfs_list=list(map("--variant {}".format, expand("VC_interm_dir/{sample}.annotated.g.vcf", sample=samples.index))),
        extra="--merge-input-intervals",  # optional
        java_opts="-Xmx4g",  # optional
    threads: 4
    resources:
        mem_mb=1024, cpus=4
    shell:
        "gatk --java-options '{params.java_opts}' GenomicsDBImport {params.extra} "
    	"{params.gvcfs_list} --intervals {params.intervals} "
    	"{params.db_action} {output.db}"

rule genotype_gvcfs:
    input:
        genomicsdb=rules.genomics_db_import.output.db,  # combined gvcf over multiple samples
    # N.B. gvcf or genomicsdb must be specified
    # in the latter case, this is a GenomicsDB data store
        ref=config["gatk"]["fasta"]
    output:
        vcf="VC_calls/output.vcf",
    log:
        "logs/gatk/genotypegvcfs.log"
    params:
        extra="--include-non-variant-sites --dbsnp {} -L {}".format(config["gatk"]["dbsnp"], config["gatk"]["intervals"]),  # optional
        java_opts="-Xmx4g", # optional
    threads: 4
    resources:
        mem_mb=1024, cpus=4
    shell:
    	"gatk --java-options '{params.java_opts}' GenotypeGVCFs {params.extra} "
    	"-V gendb://{input.genomicsdb} "
    	"-R {input.ref} "
    	"-O {output.vcf}"

rule select_snps:
    input:
        vcf="VC_calls/output.vcf",
        ref=config["gatk"]["fasta"],
    output:
        vcf="VC_calls/snps.vcf"
    log:
        "logs/gatk/select/snps.log"
    params:
        extra="--select-type-to-include SNP",  # optional filter arguments, see GATK docs
        java_opts="-Xmx4g", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    threads: 4
    resources:
        mem_mb=1024, cpus=4
    shell:
    	"gatk --java-options '{params.java_opts}' SelectVariants -R {input.ref} -V {input.vcf} "
    	"{params.extra} -O {output.vcf}"

rule select_indels:
    input:
        vcf="VC_calls/output.vcf",
        ref=config["gatk"]["fasta"],
    output:
        vcf="VC_calls/indels.vcf"
    log:
        "logs/gatk/select/indels.log"
    params:
        extra="--select-type-to-include INDEL",  # optional filter arguments, see GATK docs
        java_opts="-Xmx4g", # optional
    threads: 4
    resources:
        mem_mb=1024, cpus=4
    shell:
        "gatk --java-options '{params.java_opts}' SelectVariants -R {input.ref} -V {input.vcf} "
        "{params.extra} -O {output.vcf}"

rule gatk_filter_snps:
    input:
        vcf="VC_calls/snps.vcf",
        ref=config["gatk"]["fasta"],
    output:
        vcf="VC_calls/snps.filtered.vcf"
    log:
        "logs/gatk/filter/snps.log"
    params:
        filters="-filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'",
        extra="",  # optional arguments, see GATK docs
        java_opts="-Xmx4g", # optional
    threads: 4
    resources:
        mem_mb=1024, cpus=4
    shell:
        "gatk --java-options '{params.java_opts}' VariantFiltration -R {input.ref} -V {input.vcf} "
    	"{params.extra} {params.filters} -O {output.vcf}"
        
rule gatk_filter_indels:
    input:
        vcf="VC_calls/indels.vcf",
        ref=config["gatk"]["fasta"],
    output:
        vcf="VC_calls/indels.filtered.vcf"
    log:
        "logs/gatk/filter/snps.log"
    params:
        filters="-filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' ",
        extra="",  # optional arguments, see GATK docs
        java_opts="-Xmx4g", # optional
    threads: 4
    resources:
        mem_mb=1024, cpus=4
    shell:
        "gatk --java-options '{params.java_opts}' VariantFiltration -R {input.ref} -V {input.vcf} "
        "{params.extra} {params.filters} -O {output.vcf}"

rule sort_snps:
    input:
        snp="VC_calls/snps.filtered.vcf",
        picard=config["picard"]
    output:
        snp="VC_calls/snps.filtered.sorted.vcf",
    log:
        "logs/gatk/sort_and_merge/sort_snps.log"
    threads: 4
    resources: cpus=4
    shell:
        "java -jar {input.picard}/picard.jar SortVcf -I ./{input.snp} -O ./{output.snp}"

rule sort_indels:
    input:
        indel="VC_calls/indels.filtered.vcf",
        picard=config["picard"]
    output:
        indel="VC_calls/indels.filtered.sorted.vcf",
    log:
        "logs/gatk/sort_and_merge/sort_indels.log"
    threads: 4
    resources: cpus=4
    shell:
        "java -jar {input.picard}/picard.jar SortVcf -I ./{input.indel} -O ./{output.indel}"

rule merge_indels_snps:
    input:
        snp="VC_calls/snps.filtered.sorted.vcf",
        indel="VC_calls/indels.filtered.sorted.vcf",
        picard=config["picard"]
    output:
        merged="VC_calls/merged.vcf"
    log:
        "logs/gatk/sort_and_merge/merge.log"
    threads: 4
    resources: cpus=4
    shell:
        "java -jar {input.picard}/picard.jar MergeVcfs -I ./{input.snp}  -I ./{input.indel} -O ./{output.merged}"
