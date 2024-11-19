#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakefiles.utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
#samples = pd.read_csv("atac_samplesheet2.txt",sep='\t')
## Convert all columns to strings
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ['Proj','Donor','Condition','Tissue','Protocol_notes']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

## Group by mn and extract Read1 & Read2
read1 = samples.groupby('mn')['Read1'].apply(list).to_dict()
read2 = samples.groupby('mn')['Read2'].apply(list).to_dict()

## Set run summary name using helper script
runName = namer(samples, config['mergeBy'])
#runName = namer(samples, samples[mergeBy])




## Define actions on success
onsuccess:
	## Success message
	print("ATAC_pre_process completed successfully! Wahoo!")

rule all:
    input:
        # Explicitly list all expected output files
        fastq = [expand("output/fastq/{sampleName}_{read}.fastq.gz",
                       sampleName=key, read=['R1','R2'])
                for key in read1.keys()],
        qc = [expand("output/QC/{sampleName}_{read}_fastqc.{ext}",
                    sampleName=key, read=['R1', 'R2'], ext=['zip', 'html'])
              for key in read1.keys()],
        trim_reports = [expand('output/trim/{sampleName}_R{read}.fastq.gz_trimming_report.txt',
                             sampleName=key, read=['1', '2'])
                       for key in read1.keys()],
        trimmed_fastq = [expand('output/trim/{sampleName}_R{read}_val_{read}.fq.gz',
                              sampleName=key, read=['1', '2'])
                        for key in read1.keys()],
        bam = [expand("output/align/{sampleName}_sorted.bam",sampleName=key) for key in read1.keys()],
        #bam_stats = [expand("output/align/{sampleName}_stats.txt",sampleName=key) for key in read1.keys()],
        ataqv = [expand("output/ataqv/{sampleName}.ataqv.json", sampleName=key) for key in read1.keys()]




rule catReads:
    input:
        R1 = lambda wildcards: read1.get(wildcards.sampleName),
        R2 = lambda wildcards: read2.get(wildcards.sampleName)
    output:
        R1 = 'output/fastq/{sampleName}_R1.fastq.gz',
        R2 = 'output/fastq/{sampleName}_R2.fastq.gz'
    benchmark:
        'output/benchmarks/{sampleName}_catReads.tsv'
    log:
        errR1 = 'output/logs/{sampleName}_R1_catReads.err',
        errR2 = 'output/logs/{sampleName}_R2_catReads.err'
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log.errR1}
        cat {input.R2} > {output.R2} 2> {log.errR2}
        """

rule fastqc:
        input:
                QC1 = rules.catReads.output.R1,
                QC2 = rules.catReads.output.R2
        output:
                zip1 = "output/QC/{sampleName}_R1_fastqc.zip", # temp
                zip2 = "output/QC/{sampleName}_R2_fastqc.zip", # temp
                html1 = "output/QC/{sampleName}_R1_fastqc.html", # temp
                html2 = "output/QC/{sampleName}_R2_fastqc.html" # temp
        log:
                err = 'output/logs/fastqc_{sampleName}.err',
                out = 'output/logs/fastqc_{sampleName}.out'
        params:
                dir = "output/QC",
                version = config['fastqcVers']
        benchmark:
                'output/benchmarks/fastqc_{sampleName}.tsv'
        shell:
            """
            module load fastqc/{params.version}
            mkdir -p {params.dir}
            fastqc -o {params.dir} {input.QC1} {input.QC2} 1> {log.out} 2> {log.err};
            """

rule trim:
    input:
        R1 = rules.catReads.output.R1,
        R2 = rules.catReads.output.R2
    output:
        trim1 = 'output/trim/{sampleName}_R1_val_1.fq.gz',
        trim2 = 'output/trim/{sampleName}_R2_val_2.fq.gz',
        report1 = 'output/trim/{sampleName}_R1.fastq.gz_trimming_report.txt',
        report2 = 'output/trim/{sampleName}_R2.fastq.gz_trimming_report.txt'
    threads: 4
    params:
        version = config['trim_galore']
    log:
        err = 'output/logs/trim_{sampleName}.err',
        out = 'output/logs/trim_{sampleName}.out'
    shell:
        """
        module load trim_galore/{params.version}
        module load python/3.9.6
        module load pigz
        mkdir -p output/trim
        trim_galore -o output/trim --cores {threads} --paired {input.R1} {input.R2} 2> {log.err}
        """

rule align:
    input:
        trim1 = rules.trim.output.trim1,
        trim2 = rules.trim.output.trim2
    output:
        sortedBam = "output/align/{sampleName}_sorted.bam",
        stats = "output/align/{sampleName}_stats.txt"
    params:
        bwa_index = config['bwa_index'],
        bwa_version = config['bwaVers'],
        samtools_version = config['samtoolsVers']
    threads: 8
    benchmark:
        'output/benchmarks/align_{sampleName}.tsv'
    log:
        err = 'output/logs/align_{sampleName}.err',
        out = 'output/logs/align_{sampleName}.out'
    shell:
        """
        module load bwa/{params.bwa_version}
        module load samtools/{params.samtools_version}
        mkdir -p output/align

        bwa mem -t {threads} -M {params.bwa_index} {input.trim1} {input.trim2} | \
        samtools view -q 30 -b | \
        samtools sort -o {output.sortedBam} 1> {log.out} 2> {log.err}

        samtools flagstat {output.sortedBam} > {output.stats} 2>> {log.err}
        """

rule mark_duplicates:
    input:
        bam="output/align/{sampleName}_sorted.bam"
    output:
        dedup_bam="output/dedup/{sampleName}_dedup.bam",
        metrics="output/dedup/{sampleName}_dup_metrics.txt",
        index="output/dedup/{sampleName}_dedup.bai"
    params:
        picard_version=config["picardVers"],
        java_version=config['javaVers']
    log:
        err="output/logs/mark_duplicates_{sampleName}.err",
        out="output/logs/mark_duplicates_{sampleName}.out"
    shell:
        """
        module load java/{params.java_version}
        module load picard/{params.picard_version}
        mkdir -p output/dedup

        java -Xmx16g -jar /nas/longleaf/apps/picard/{params.java_version}/picard-{params.java_version}/picard.jar MarkDuplicates \
            I={input.bam} \
            O={output.dedup_bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            1> {log.out} 2> {log.err}

        samtools index {output.dedup_bam} 1>> {log.out} 2>> {log.err}
        """

rule filter_mitochondrial_reads:
    input:
        dedup_bam=rules.mark_duplicates.output.dedup_bam,
        index=rules.mark_duplicates.output.index
    output:
        filtered_bam="output/filtered/{sampleName}_filtered.bam",
        index="output/filtered/{sampleName}_filtered.bam.bai"
    params:
        samtools_version=config['samtoolsVers']
    log:
        err="output/logs/filter_mitochondrial_{sampleName}.err",
        out="output/logs/filter_mitochondrial_{sampleName}.out"
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p output/filtered

        samtools view -bh {input.dedup_bam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
        -o {output.filtered_bam} 1> {log.out} 2> {log.err}

        samtools index {output.filtered_bam} 1>> {log.out} 2>> {log.err}
        """

rule collect_insert_size_metrics:
    input:
        bam=rules.filter_mitochondrial_reads.output.filtered_bam
    output:
        metrics="output/metrics/{sampleName}_insert_size_metrics.txt",
        histogram="output/metrics/{sampleName}_insert_size_histogram.pdf"
    params:
        picard_version=config["picardVers"],
        java_version=config['javaVers']
    log:
        err="output/logs/insert_size_metrics_{sampleName}.err",
        out="output/logs/insert_size_metrics_{sampleName}.out"
    shell:
        """
        module load java/{params.java_version}
        module load picard/{params.picard_version}
        mkdir -p output/metrics

        java -jar /nas/longleaf/apps/picard/{params.picard_version}/picard.jar CollectInsertSizeMetrics \
            I={input.bam} \
            O={output.metrics} \
            H={output.histogram} \
            M=0.05 \
            ASSUME_SORTED=true \
            1> {log.out} 2> {log.err}
        """

rule add_read_groups:
    input:
        bam="output/filtered/{sampleName}_filtered.bam",
        bai="output/filtered/{sampleName}_filtered.bam.bai"
    output:
        bam="output/align/{sampleName}_RG.bam",
        bai="output/align/{sampleName}_RG.bai"
    params:
        picard_version=config["picardVers"],
        samtools_version=config["samtoolsVers"]
    log:
        err="output/logs/add_read_groups_{sampleName}.err",
        out="output/logs/add_read_groups_{sampleName}.out"
    shell:
        """
        module load picard/{params.picard_version}
        module load samtools/{params.samtools_version}
        mkdir -p output/align

        picard AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.bam} \
            RGSM={wildcards.sampleName} \
            RGPL=ILLUMINA \
            RGLB=lib1 \
            RGPU=unit1 \
            1> {log.out} 2> {log.err}

        samtools index {output.bam} 1>> {log.out} 2>> {log.err}
        """

rule run_ataqv:
    input:
        bam="output/align/{sampleName}_RG.bam"
    output:
        json="output/ataqv/{sampleName}.ataqv.json"
    params:
        tss_file="/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.tss.bed.gz",
        ataqv_version=config["ataqvVers"]
    log:
        err="output/logs/ataqv_{sampleName}.err",
        out="output/logs/ataqv_{sampleName}.out"
    shell:
        """
        module load ataqv/{params.ataqv_version}
        mkdir -p output/ataqv

        ataqv \
            --tss-file {params.tss_file} \
            --ignore-read-groups human \
            {input.bam} > {output.json} \
            2> {log.err}
        """
























 

    
      
    



