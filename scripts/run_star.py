import subprocess

input = snakemake.input
output = snakemake.output
params = snakemake.params
threads = snakemake.threads
log = snakemake.log
id = snakemake.wildcards.id


def run_star(index_dir,input1,input2,threads,outdir,stringence,log,zcat):

  Star_list = [
    'STAR', '--runMode', 'alignReads',
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
    '--alignEndsProtrude', '5', 'ConcordantPair',
    '--outFilterMismatchNmax', f'{stringence}',
    '&>', f'{log}'
  ]

  if zcat :
    del Star_list[10]
    del Star_list[11]

  subprocess.call(Star_list)

# Run the first time
run_star(
  f"{params.index_dir}",
  f"{input.fastq1}",
  f"{input.fastq2}",
  f"{threads}",
  f"{params.outdir0}",
  0,
  f"{log}",
  False
)

# Take the Unmapped for other rounds
for i in range(1,4):

  run_star(
    f"{params.index_dir}",
    f"{params.outdir2}{i-1}/{id}/Unmapped.out.mate1",
    f"{params.outdir2}{i-1}/{id}/Unmapped.out.mate2",
    f"{threads}",
    f"{params.outdir2}{i}/{id}/",
    i,
    f"{log}{i}/{id}/",
    True
  )
