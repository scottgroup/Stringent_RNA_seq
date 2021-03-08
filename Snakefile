import os
from pathlib import Path


configfile: "config.json"

original_id = list(config['dataset'].values())
simple_id = list(config['dataset'].keys())


rule all:
    input:
        coco_cc = expand('results/coco/{id}.tsv', id=simple_id),
        bam_merge =expand("results/star/{id}/Merge_bam.bam",id=simple_id),
        merge_tpm = "results/coco/merge_tpm_final.tsv",
        bam = expand("results/star/0/{id}/Aligned.sortedByCoord.out.bam", id=simple_id)

rule all_downloads:
    input:
        samples = expand('data/references/fastq/{id}_{pair}.fastq.gz',
            id=original_id, pair=[1, 2]),
        reference_genome = config['path']['reference_genome'],
        reference_gtf = config['path']['reference_annotation'],
        coco_git = 'git_repos/coco'


rule rename_samples:
    """Rename samples with a nicer understandable name"""
    input:
        fastq = expand(
            "data/references/fastq/{id}_{pair}.fastq.gz",
            id=original_id, pair=[1, 2])
    output:
        renamed_fastq = expand(
            "data/references/fastq/{id}_R{pair}.fastq.gz",
            id=simple_id, pair=[1, 2])
    run:
        for new_name, old_name in config['dataset'].items():
            for num in [1, 2]:
                old = "data/references/fastq/{}_{}.fastq.gz".format(old_name, num)
                new = "data/references/fastq/{}_R{}.fastq.gz".format(new_name, num)
                print(old, new)
                os.rename(old, new)


rule trimming:
    """Trims the input FASTQ files using Trimmomatic"""
    input:
        fastq1 = "data/references/fastq/{id}_R1.fastq.gz",
        fastq2 = "data/references/fastq/{id}_R2.fastq.gz"
    output:
        fastq1 = "data/Trimmomatic/trimmed_reads/{id}_R1.fastq.gz",
        fastq2 = "data/Trimmomatic/trimmed_reads/{id}_R2.fastq.gz",
        unpaired_fastq1 = "data/Trimmomatic/trimmed_reads/{id}_R1.unpaired.fastq.gz",
        unpaired_fastq2 = "data/Trimmomatic/trimmed_reads/{id}_R2.unpaired.fastq.gz"
    threads:
        32
    params:
        options = [
            "ILLUMINACLIP:data/Trimmomatic/Adapters-PE_NextSeq.fa:2:12:10:8:true",
            "TRAILING:30", "LEADING:30", "MINLEN:20"
        ]
    log:
        "logs/trimmomatic/{id}.log"
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fastq1} {input.fastq2} "
        "{output.fastq1} {output.unpaired_fastq1} "
        "{output.fastq2} {output.unpaired_fastq2} "
        "{params.options} "
        "&> {log}"

rule coco_ca:
    """ Correct gtf annotation with CoCo correct_annotation for embedded genes"""
    input:
        gtf_ca = config['path']['reference_annotation']
    output:
        gtf = config['path']['reference_annotation_ca'],
        intron_gtf = config['path']['reference_intron_annotation']
    conda:
        "envs/coco.yaml"
    params:
        coco_path = "git_repos/coco/bin"
    shell:
        "python {params.coco_path}/correct_annotation.py "
        "{input.gtf_ca}"

rule star_index:
    """Generate the genome index needed for STAR alignment"""
    input:
        fasta = config['path']['reference_genome'],
        gtf = config['path']['reference_annotation']
    output:
        chrNameLength = "data/star_index/chrNameLength.txt"
    threads:
        32
    params:
        index_dir = "data/star_index/"
    log:
        "logs/star/star_index.log"
    conda:
        "envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 74 "
        "&> {log}"

rule star_align:
    """Align reads to reference genome using STAR"""
    input:
        fastq1 = rules.trimming.output.fastq1,
        fastq2 = rules.trimming.output.fastq2,
        idx = rules.star_index.output
    output:
        bam = "results/star/0/{id}/Aligned.sortedByCoord.out.bam",
        bam1 = "results/star/1/{id}/Aligned.sortedByCoord.out.bam",
        bam2= "results/star/2/{id}/Aligned.sortedByCoord.out.bam",
        bam3 = "results/star/3/{id}/Aligned.sortedByCoord.out.bam",
    threads:
        16
    params:
        index_dir = "data/star_index/",
        outdir2 = "results/star/",
        outdir0 = "results/star/0/{id}/"
    log:
        "logs/star/0/star_align_{id}.log"
    conda:
        "envs/star.yaml"
    script:
        "scripts/run_star.py"

rule Merge_BAM:
    """ Merge the BAM files before going to coco"""
    input:
        bam = expand("results/star/{wild}/{{id}}/Aligned.sortedByCoord.out.bam",wild=[0,1,2,3])
    output:
        bam_merge = "results/star/{id}/Merge_bam.bam"
    threads:
        4
    params:
        outdir = "results/star/{id}/",
        outdir_star = "results/star/0/{id}/"
    conda:
        "envs/star.yaml"
    shell:
        "samtools merge {output.bam_merge} {input.bam} "

rule coco_cc:
    """ Quantify the number of counts, counts per million (CPM) and transcript
        per million (TPM) for each gene using CoCo correct_count (cc). If you
        want to use your own annotation file, you must implement in this
        workflow the correct_annotation function of CoCo in order to use CoCo's
        correct_count function. Otherwise, by default, you will use the given
        annotation file used in this analysis """
    input:
        gtf = config['path']['reference_annotation_ca'],
        intron_gtf = config['path']['reference_intron_annotation'],
        bam = rules.Merge_BAM.output.bam_merge
    output:
        counts = Path("results/coco/", "{id}.tsv")
    threads:
        16
    params:
        coco_path = "git_repos/coco/bin"
    log:
        "logs/coco/coco_{id}.log"
    conda:
        "envs/coco.yaml"
    shell:
        "python {params.coco_path}/coco.py cc "
        "--countType both "
        "--thread {threads} "
        "--strand 1 "
        "--paired "
        "{input.gtf} "
        "{input.bam} "
        "{output.counts} "
        "&> {log}"

rule merge:
    input:
        counts = expand("results/coco/{id}.tsv",id=simple_id)
    output:
        merge_tpm = "results/coco/merge_tpm_final.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/merge_tpm.py"


#Include download_all rules to download necessary datasets and references
include: "rules/download_all.smk"
