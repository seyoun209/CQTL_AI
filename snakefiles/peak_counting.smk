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
        expand("output/peaks/{tissue}/{cond}_consensus_macs3_peaks.bed", cond=condition_samples.keys(),tissue=tissue),
        expand("output/peaks/{tissue}/{cond}_consensus_hmmratac_peaks.bed", cond=condition_samples.keys(),tissue=tissue),
        expand("output/peaks/{tissue}/{cond}_consensus_macs2_peaks.bed", cond=condition_samples.keys(),tissue=tissue),
        expand("output/peaks/{tissue}/merged/all_macs3_merged.bed", tissue=tissue),
        expand("output/peaks/{tissue}/merged/all_macs2_merged.bed", tissue=tissue),
        expand("output/peaks/{tissue}/merged/all_hmmratac_merged.bed", tissue=tissue),
        expand("output/peaks/{tissue}/counts/{tooltype}.saf",tooltype=['macs3', 'macs2','hhmratac'],tissue=tissue),
        expand("output/peaks/{tissue}/counts/{tooltype}_counts.txt",tooltype=['macs3', 'macs2','hmmratac'],tissue=tissue)

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

rule macs3_merge_peaks_by_condition:
    input:
        peaks=lambda wildcards: expand("output/peaks/{tissue}/macs3/{sampleName}_peaks.narrowPeak", 
                tissue=wildcards.tissue,                       
                sampleName=condition_samples[wildcards.cond])
    output:
        consensus_peaks="output/peaks/{tissue}/{cond}_consensus_macs3_peaks.bed"
    log:
        err="output/logs/{tissue}_macs3_merge_bedtools_{cond}.err"
    params:
        min_samples=lambda wildcards: max(1, len(condition_samples[wildcards.cond]) // 2),
        bedtoolsVer=config['bedtools']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}

        mkdir -p output/peaks/{wildcards.tissue}

        # Merge peaks that appear in at least 50% of samples
        multiIntersectBed -i {input.peaks} | \
        awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
        bedtools merge > {output.consensus_peaks} 2> {log.err}

        """

rule macs2_merge_peaks_by_condition:
    input:
        peaks=lambda wildcards: expand("output/peaks/{tissue}/macs2/{sampleName}_peaks.narrowPeak",
                tissue=wildcards.tissue,
                sampleName=condition_samples[wildcards.cond])
    output:
        consensus_peaks="output/peaks/{tissue}/{cond}_consensus_macs2_peaks.bed"
    log:
        err="output/logs/{tissue}_macs2_merge_bedtools_{cond}.err"
    params:
        min_samples=lambda wildcards: max(1, len(condition_samples[wildcards.cond]) // 2),
        bedtoolsVer=config['bedtools']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}

        mkdir -p output/peaks/{wildcards.tissue}

        # Merge peaks that appear in at least 50% of samples
        multiIntersectBed -i {input.peaks} | \
        awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
        bedtools merge > {output.consensus_peaks} 2> {log.err}

        """
rule hmmratac_merge_peaks_by_condition:
    input:
        peaks=lambda wildcards: expand("output/peaks/{tissue}/hmmratac/{sampleName}_accessible_regions.narrowPeak",
                tissue=wildcards.tissue,
                sampleName=condition_samples[wildcards.cond])
    output:
        consensus_peaks="output/peaks/{tissue}/{cond}_consensus_hmmratac_peaks.bed"
    log:
        err="output/logs/{tissue}_hmmratac_merge_bedtools_{cond}.err"
    params:
        min_samples=lambda wildcards: max(1, len(condition_samples[wildcards.cond]) // 2),
        bedtoolsVer=config['bedtools']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}

        mkdir -p output/peaks/{wildcards.tissue}

        # Merge peaks that appear in at least 50% of samples
        multiIntersectBed -i {input.peaks} | \
        awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
        bedtools merge > {output.consensus_peaks} 2> {log.err}
        """
rule merge_all_peaks_by_caller:
    input:
        macs3_peaks=expand("output/peaks/{tissue}/macs3/{sampleName}_peaks.narrowPeak", tissue=tissue, sampleName=bamfile.keys()),
        macs2_peaks=expand("output/peaks/{tissue}/macs2/{sampleName}_peaks.narrowPeak", tissue=tissue, sampleName=bamfile.keys()),
        hmmratac_peaks=expand("output/peaks/{tissue}/hmmratac/{sampleName}_accessible_regions.narrowPeak", tissue=tissue, sampleName=bamfile.keys())
    output:
        macs3_merged="output/peaks/{tissue}/merged/all_macs3_merged.bed",
        macs2_merged="output/peaks/{tissue}/merged/all_macs2_merged.bed",
        hmmratac_merged="output/peaks/{tissue}/merged/all_hmmratac_merged.bed"
    log:
        err="output/logs/{tissue}_merge_peaks.err",
        out="output/logs/{tissue}_merge_peaks.out"
    params:
        bedtoolsVer=config['bedtools']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        mkdir -p output/peaks/{wildcards.tissue}/merged

        # MACS3
        cat {input.macs3_peaks} | awk '{{ OFS="\t"}} {{ print $1, $2, $3 }}' | \
        sort -k1,1 -k2,2n | bedtools merge > {output.macs3_merged} 2>> {log.err}

        # MACS2
        cat {input.macs2_peaks} | awk '{{ OFS="\t"}} {{ print $1, $2, $3 }}' | \
        sort -k1,1 -k2,2n | bedtools merge > {output.macs2_merged} 2>> {log.err}

        # HMMRATAC
        cat {input.hmmratac_peaks} | awk '{{ OFS="\t"}} {{ print $1, $2, $3 }}' | \
        sort -k1,1 -k2,2n | bedtools merge > {output.hmmratac_merged} 2>> {log.err}

        # Log peak counts
        wc -l {output.macs3_merged} {output.macs2_merged} {output.hmmratac_merged} > {log.out}
        """

rule create_count_matrices:
    input:
        macs3_peaks=rules.merge_all_peaks_by_caller.output.macs3_merged,
        macs2_peaks=rules.merge_all_peaks_by_caller.output.macs2_merged,
        hmmratac_peaks=rules.merge_all_peaks_by_caller.output.hmmratac_merged,
        bams=expand("{sample}", sample=[bam for sample in bamfile.keys() 
                                      for bam in bamfile.get(sample)])
    output:
        saf_macs3="output/peaks/{tissue}/counts/macs3.saf",
        saf_macs2="output/peaks/{tissue}/counts/macs2.saf",
        saf_hhmratac="output/peaks/{tissue}/counts/hhmratac.saf",
        macs3_counts="output/peaks/{tissue}/counts/macs3_counts.txt",
        macs2_counts="output/peaks/{tissue}/counts/macs2_counts.txt",
        hmmratac_counts="output/peaks/{tissue}/counts/hmmratac_counts.txt"
    log:
        err="output/logs/{tissue}_count_matrices.err",
        out="output/logs/{tissue}_count_matrices.out"
    threads: 8
    params:
        subreadVer=config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{wildcards.tissue}/counts

        # MACS3
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}}
             {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {input.macs3_peaks} > {output.saf_macs3}
        featureCounts -p -T {threads} -F SAF -a {output.saf_macs3} -o {output.macs3_counts} {input.bams} 2>> {log.err}

        # MACS2
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}}
             {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {input.macs2_peaks} > {output.saf_macs2}
        featureCounts -p -T {threads} -F SAF -a {output.saf_macs2} -o {output.macs2_counts} {input.bams} 2>> {log.err}

        # HMMRATAC
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}}
             {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {input.hmmratac_peaks} > {output.saf_hhmratac}
        featureCounts -p -T {threads} -F SAF -a {output.saf_hhmratac} -o {output.hmmratac_counts} {input.bams} 2>> {log.err}
        """
