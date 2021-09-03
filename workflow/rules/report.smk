rule final_report:
    input:
        "optitype_result/{sample}_result.tsv",
        "VC_calls/output.vcf",
        "VC_calls/merged.vcf",
        "mosdepth_bed/{sample}.regions.bed.gz",
        list=expand(["aldy_result/{{sample}}.{gene}.aldy"], gene=genes)
    output:
    	"reports/{sample}_results.tsv",
    	"reports/{sample}_multipleSNP_results.tsv"
    params:
    	sample="{sample}",
    	variants=config["report"]["variants"]
    log:
       "logs/report/{sample}.log"
    conda:
       "../envs/report_results.yml"
    resources: mem_mb=100000
    script:
    	"../scripts/report_results_snake.py"
