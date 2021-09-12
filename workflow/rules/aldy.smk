rule aldy:
   input: 
       bam="preprocessed_VC/{sample}.bam"
   output:
       aldy_out="aldy_result/{sample}.{gene}.aldy"
   log:
       "logs/aldy/{sample}.{gene}.log"
   conda:
       "../envs/aldy_env.yaml"
   resources: mem_mb=100000
   shell:
   	"aldy genotype -p illumina -g {wildcards.gene} -o {output.aldy_out} {input.bam} &> {log}"
