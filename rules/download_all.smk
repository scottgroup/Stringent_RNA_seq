import os

rule download_sample_fastq:
    """Download expression datasets of all tissue samples from GEO. This step
        might take long, depending on your downloading speed."""
    output:
        samples_fastq_1_gz = "data/references/fastq/{id}_1.fastq.gz",
        samples_fastq_2_gz = "data/references/fastq/{id}_2.fastq.gz"
    params:
        samples_fastq_1 = "data/references/fastq/{id}_1.fastq",
        samples_fastq_2 = "data/references/fastq/{id}_2.fastq"
    conda:
        "../envs/geo_download.yaml"
    shell:
        "mkdir -p data/references/fastq/ && "
        "fasterq-dump --split-files {wildcards.id} "
        "-O data/references/fastq/ && "
        "gzip -f {params.samples_fastq_1} > {output.samples_fastq_1_gz} && "
        "gzip -f {params.samples_fastq_2} > {output.samples_fastq_2_gz}"

rule download_genome:
    """Download the reference genome (fasta file) used for this analysis from
        ENSEMBL ftp servers."""
    output:
        genome = config['path']['reference_genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gunzip {output.genome}.gz"

rule download_annotation:
    """Download the annotation (gtf file) used for this analysis."""
    output:
        gtf = config['path']['reference_annotation'],
    params:
        link_annotation = config['download']['annotation'],
    shell:
        "wget --quiet -O {output.gtf} {params.link_annotation}"

rule download_coco_git:
    """Download git repository of CoCo."""
    output:
        git_coco_folder = directory('git_repos/coco')
    params:
        git_coco_link = config['path']['coco_git_link']
    conda:
        '../envs/git.yaml'
    shell:
        'mkdir -p {output.git_coco_folder} '
        '&& git clone {params.git_coco_link} {output.git_coco_folder}'
