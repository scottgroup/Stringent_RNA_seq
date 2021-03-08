from snakemake.shell import shell
id = snakemake.wildcards.id


shell(f"touch {snakemake.output.bam_merge} ")
shell(f'cp {snakemake.input.bam[0]} {snakemake.output.bam_merge}')

for i in range(1,4):
    shell(f"samtools view -h -b {snakemake.input.bam[i]}  >> {snakemake.output.bam_merge} ")
    

