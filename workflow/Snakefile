import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "../config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
ids= [1,2]
genes= ["CFTR","CYP2C9","CYP4F2","CYP3A5","DPYD","NAT2","NUDT15","SLCO1B1","TPMT","G6PD"]

##### target rules #####

rule all_mapping:
    input:
        expand("mapped/{sample}.bam", sample=samples.index),
        expand("mapped/{sample}.sorted.bam", sample=samples.index)
rule all_optitype:
    input:
      expand("preprocessed_optitype/{sample}.fished_{id}.fastq", sample=samples.index, id=ids),
      expand(["optitype_result/{sample}_coverage_plot.pdf", "optitype_result/{sample}_result.tsv"], sample=samples.index)
rule all_VC_preprocessing:
    input:
      expand(["VC_recal/dedup_{sample}.sorted.bam","VC_recal/{sample}.metrics.txt", "VC_recal/dedup_{sample}.sorted.bam.bai","VC_recal/{sample}.grp", \
      "preprocessed_VC/{sample}.bam"], sample=samples.index)
rule all_VC:
    input:
      expand(["preprocessed_VC/{sample}.g.vcf.gz","VC_interm_dir/{sample}.annotated.g.vcf"] , sample=samples.index),
      "VC_calls/output.vcf",
      "VC_calls/snps.vcf",
      "VC_calls/indels.vcf",
      "VC_calls/snps.filtered.vcf",
      "VC_calls/indels.filtered.vcf",
      "VC_calls/snps.filtered.sorted.vcf",
      "VC_calls/indels.filtered.sorted.vcf",
      "VC_calls/merged.vcf"
rule all_aldy:
    input:
      expand("aldy_result/{sample}.{gene}.aldy", sample=samples.index, gene=genes)
rule all_mosdepth:
    input:
      expand("mosdepth_bed/{sample}.mosdepth.global.dist.txt", sample=samples.index),
      expand("mosdepth_bed/{sample}.mosdepth.region.dist.txt", sample=samples.index),
      expand("mosdepth_bed/{sample}.regions.bed.gz", sample=samples.index)

rule all_report:
    input:
      expand(["reports/{sample}_results.tsv", "reports/{sample}_multipleSNP_results.tsv"], sample=samples.index)
##### load rules #####

include: "rules/mapping.smk"
include: "rules/mosdepth.smk"
include: "rules/optitype.smk"
include: "rules/VC_preprocessing.smk"
include: "rules/variant_calling.smk"
include: "rules/aldy.smk"
include: "rules/report.smk"
