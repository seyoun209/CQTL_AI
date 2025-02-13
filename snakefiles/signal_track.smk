#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakefiles.utils.namer import namer

# Read in samplesheet
samples = pd.read_csv(config["bam_samplesheet"], sep='\t')
#samples = pd.read_csv("bam_atac_samplsheet.txt",sep='\t')
samples = samples.astype(str)

# Create mn column
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ['Proj','Donor','Condition','Tissue','Protocol_notes']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

# Group by mn and bam
samples['bam'] = samples[['Bam_directory', 'Bam_file']].apply(lambda row: os.path.join(*row), axis=1)
bamfile = samples.groupby('mn')['bam'].apply(list).to_dict()

# Build dictionary of merged samples by condition
mergeSample = samples.groupby('Condition')['mn'].apply(list).to_dict()

# Build dictionary for tissue grouping
# Get unique tissues
tissue = samples['Tissue'].iloc[0]

##### Define rules #####

rule all:
    input:
        expand("output/signals/{tissue}/{sampleName}.bw",tissue=samples['Tissue'].unique(),sampleName=bamfile.keys()),
        expand('output/signals/{tissue}/merged_signal/{cond}_sorted_bw',cond=mergeSample.keys(),tissue=tissue)

rule signal:
    input:
        bam = lambda wildcards: bamfile.get(wildcards.sampleName)
    output:
        signal = "output/signals/{tissue}/{sampleName}.bw"
    log:
        err = 'output/logs/signal_{tissue}_{sampleName}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        binSize=config['bin_size'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/{wildcards.tissue}

        bamCoverage \
            --bam {input.bam} \
            -o {output.signal} \
            --binSize {params.binSize} \
            --normalizeUsing {params.NormOption} \
            --effectiveGenomeSize {params.effective_genomeSize} \
            --extendReads \
            > {log.err} 2>&1

        """

rule mergeAlign:
    input:
        bam_files = lambda wildcards: [bamfile[sample][0] for sample in mergeSample[wildcards.cond]]
    output:
        bam = "output/signals/{tissue}/mergedAlign/{cond}_sorted.bam",
        bai = "output/signals/{tissue}/mergedAlign/{cond}_sorted.bam.bai",
        stats = "output/signals/{tissue}/mergedAlign/{cond}_stats.txt"
    log:
        err='output/logs/mergedAligned_{tissue}_{cond}.err'
    params:
        samtools_version = config['samtoolsVers']
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p output/signals/{wildcards.tissue}/mergedAlign

        samtools merge {output.bam} {input.bam_files} >> {log.err} 2>&1
        samtools flagstat {output.bam} > {output.stats} >> {log.err} 2>&1
        samtools index {output.bam} >> {log.err} 2>&1
        """

rule mergeSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal='output/signals/{tissue}/merged_signal/{cond}_sorted_bw'
    log:
        err='output/logs/mergedSignal_{tissue}_{cond}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        binSize=config['bin_size'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/{wildcards.tissue}/merged_signal

        bamCoverage \
            --bam {input.bam} \
            -o {output.signal} \
            --binSize {params.binSize} \
            --normalizeUsing {params.NormOption} \
            --effectiveGenomeSize {params.effective_genomeSize} \
            --extendReads \
            > {log.err} 2>&1

        """

#rule mergeForwardSigInal:
#    input:
#        bam= rules.mergeAlign.output.bam
#    output:
#        signal="output/signals/{tissue}/strand/{cond}_fwd.bw"
#    log:
#        err="output/logs/mergedSignal_{tissue}_{cond}_forward.err"
#    params:
#        deeptools_ver=config['deeptoolsVers']
#    shell:
#        """
#        module load deeptools/{params.deeptools_ver}
#        mkdir -p output/signals/{wildcards.tissue}/strand
#
#        bamCoverage 
#        """
