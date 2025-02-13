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

# Group samples by condition
condition_samples = samples.groupby('Condition')['mn'].apply(list).to_dict()

tissue = samples['Tissue'].iloc[0]

##### Define rules #####

rule all:
    input:
        expand("output/peaks/{tissue}/macs2/{sampleName}_peaks.narrowPeak",sampleName=bamfile.keys(),tissue=tissue),
        expand("output/peaks/{tissue}/macs3/{sampleName}_peaks.narrowPeak",sampleName=bamfile.keys(),tissue=tissue),
        expand("output/peaks/{tissue}/hmmratac/{sampleName}_accessible_regions.narrowPeak",sampleName=bamfile.keys(),tissue=tissue),
        #expand("output/peaks/{tissue}/{cond}_consensus_macs3_peaks.bed", cond=condition_samples.keys(),tissue=tissue),
        #expand("output/peaks/{tissue}/{cond}_consensus_hmmratac_peaks.bed", cond=condition_samples.keys(),tissue=tissue),
        #expand("output/peaks/{tissue}/{cond}_consensus_macs2_peaks.bed", cond=condition_samples.keys(),tissue=tissue),
        #expand("output/peaks/{tissue}/merged/all_macs3_merged.bed", tissue=tissue),
        #expand("output/peaks/{tissue}/merged/all_macs2_merged.bed", tissue=tissue),
        #expand("output/peaks/{tissue}/merged/all_hmmratac_merged.bed", tissue=tissue),
        #expand("output/peaks/{tissue}/counts/{tooltype}.saf",tooltype=['macs3', 'macs2','hhmratac'],tissue=tissue),
        #expand("output/peaks/{tissue}/counts/{tooltype}_counts.txt",tooltype=['macs3', 'macs2','hmmratac'],tissue=tissue)

rule macs3_peaks:
    input:
        bam=lambda wildcards: bamfile.get(wildcards.sampleName)
    output:
        peak="output/peaks/{tissue}/macs3/{sampleName}_peaks.narrowPeak"
    log:
        err = 'output/logs/{tissue}_macs3_{sampleName}.err'
    threads:2
    params:
        python_ver = config['python']
    shell:
        """
        module load python/{params.python_ver}
        source /users/s/e/seyoun/tools/MACS3env/bin/activate

        mkdir -p output/peaks/{wildcards.tissue}/macs3/

        macs3 callpeak -t {input.bam} \
                -f BAMPE \
                -g hs \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                --keep-dup all \
                -q 0.01 \
                --call-summits \
                -n {wildcards.sampleName} \
                --outdir output/peaks/{wildcards.tissue}/macs3 2> {log.err}
        """
rule hmmratac_peaks:
    input:
        bam=lambda wildcards: bamfile.get(wildcards.sampleName)
    output:
        peak="output/peaks/{tissue}/hmmratac/{sampleName}_accessible_regions.narrowPeak"
    log:
        err='output/logs/{tissue}_hmmratac_{sampleName}.err'
    params:
        python_ver=config['python']
    threads:2
    shell:
        """
        module load python/{params.python_ver}
        source /users/s/e/seyoun/tools/MACS3env/bin/activate

        mkdir -p output/peaks/{wildcards.tissue}/hmmratac

         macs3 hmmratac \
            -i {input.bam} \
            -f BAMPE \
            --keep-duplicates \
            --minlen 100 \
            --save-states \
            --binsize 50 \
            -n {wildcards.sampleName} \
            --outdir output/peaks/{wildcards.tissue}/hmmratac 2> {log.err}

        """

rule macs2_peaks:
    input:
        bam=lambda wildcards: bamfile.get(wildcards.sampleName)
    output:
        peak="output/peaks/{tissue}/macs2/{sampleName}_peaks.narrowPeak"
    log:
        err='output/logs/{tissue}_macs2_{sampleName}.err'
    threads:2
    params:
        python_ver=config['python'],
        macs2_ver=config['macs2']
    shell:
        """
        module load python/{params.python_ver}
        module load macs/{params.macs2_ver}
        
        mkdir -p output/peaks/{wildcards.tissue}/macs2

        macs2 callpeak -t {input.bam} \
            -f BAMPE \
            -g hs \
            --nomodel \
            --shift -75 \
            --extsize 150 \
            --keep-dup all \
            -p 0.01 \
            --call-summits \
            -n {wildcards.sampleName} \
            --outdir output/peaks/{wildcards.tissue}/macs2 2> {log.err}
        """
