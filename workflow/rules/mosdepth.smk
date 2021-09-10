rule mosdepth_bed:
    input:
        bam="mapped/{sample}.sorted.bam",
        #bai="mapped/{sample}.sorted.bam.bai",
        bed=config["mosdepth"]["bed"],
    output:
        "mosdepth_bed/{sample}.mosdepth.global.dist.txt",
        "mosdepth_bed/{sample}.mosdepth.region.dist.txt",
        "mosdepth_bed/{sample}.regions.bed.gz",
        summary="mosdepth_bed/{sample}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth_bed/{sample}.log",
    params:
        extra="--no-per-base --fast-mode",  # optional
    # additional decompression threads through `--threads`
    threads: 2  # This value - 1 will be sent to `--threads`
    resources: cpus=2, mem_mb=20000
    wrapper:
        "0.77.0/bio/mosdepth"
