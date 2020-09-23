import subprocess
input = snakemake.input
output = snakemake.output
params = snakemake.params
threads = snakemake.threads
log = snakemake.log
id = snakemake.wildcards.id

def run_star(index_dir,input1,input2,threads,outdir,stringence,log,zcat):

  
  Star_list = ['STAR', '--runMode', 'alignReads',
          '--genomeDir', f'{index_dir}',
          '--readFilesIn', f'{input1}', f'{input2}',
          '--runThreadN', f'{threads}',
          '--readFilesCommand', 'zcat',
          '--outReadsUnmapped', 'Fastx',
          '--outFilterType', 'BySJout',
          '--outStd', 'Log',
          '--outSAMunmapped', 'None',
          '--outSAMtype', 'BAM', 'SortedByCoordinate',
          
          '--outFileNamePrefix', f'{outdir}',
          '--outFilterScoreMinOverLread', '0.3',
          '--outFilterMatchNminOverLread', '0.3',
          '--outFilterMultimapNmax', '100',
          '--winAnchorMultimapNmax', '100',
          '--alignEndsProtrude', f'{stringence}', 'ConcordantPair',
          '&>', f'{log}']
  if zcat :

    del Star_list[10]
    del Star_list[11]

  
 
  subprocess.call(Star_list)
  #subprocess.call(f'STAR --runMode alignReads --genomeDir {index_dir} --readFilesIn {input1} {input2}')

run_star(f"{params.index_dir}",f"{input.fastq1}",f"{input.fastq2}",f"{threads}",f"{params.outdir0}",0,f"{log}",False)
#subprocess.call(f'mv "{params.outdir2}"0"/{id}/" {params.outdir}Merge_Bam.bam ',shell=True)

for i in range(1,4):
  run_star(f"{params.index_dir}",f"{params.outdir2}"+str(i-1)+f"/{id}/"+"Unmapped.out.mate1",f"{params.outdir2}"+str(i-1)+f"/{id}/"+"Unmapped.out.mate2",f"{threads}",f"{params.outdir2}"+str(i)+f"/{id}/",i,f"{log}"+str(i)+f"/{id}/",True)
  #subprocess.call(f' samtools view {output.bam_merge} {output.bam} > {output.bam_merge} ',shell=True)


"""
Star_list = ['STAR --runMode alignReads',
          f'--genomeDir {index_dir}',
          f'--readFilesIn {input1} {input2} ',
          f'--runThreadN {threads} ',
          '--readFilesCommand zcat ',
          '--outReadsUnmapped Fastx ',
          '--outFilterType BySJout ',
          '--outStd Log ',
          '--outSAMunmapped None ',
          '--outSAMtype BAM SortedByCoordinate ',
          '--outSAMtype BAM Unsorted ',
          f'--outFileNamePrefix {outdir} ',
          '--outFilterScoreMinOverLread 0.3 ',
          '--outFilterMatchNminOverLread 0.3 ',
          '--outFilterMultimapNmax 100 ',
          '--winAnchorMultimapNmax 100 ',
          '--alignEndsProtrude {stringence} ConcordantPair ',
          f'&> {log}']
"""