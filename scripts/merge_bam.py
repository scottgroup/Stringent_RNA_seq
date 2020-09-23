from snakemake.shell import shell
id = snakemake.wildcards.id


shell(f"touch {snakemake.params.outdir}Merge_bam.bam ")
print(f"samtools view -h {snakemake.input.bam[0]} >> {snakemake.output.bam_merge}")
#shell(f"samtools view -h {snakemake.input.bam[0]} >> {snakemake.output.bam_merge}")
for i in range(1,4):
    shell(f"samtools view {snakemake.input.bam[i]} >> {snakemake.output.bam_merge}")

