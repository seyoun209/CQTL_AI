#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os

# Read in samplesheet
samples = pd.read_csv(config["bam_samplesheet"], sep='\t')
samples = samples.astype(str)

# Create mn column using mergeBy from config
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

# Group by mn and bam
samples['bam'] = samples[['Bam_directory', 'Bam_file']].apply(lambda row: os.path.join(*row), axis=1)
bamfile = samples.groupby('mn')['bam'].apply(list).to_dict()

# Group samples by condition
condition_samples = samples.groupby('Condition')['mn'].apply(list).to_dict()

# Set tissue explicitly to "Ankle"
tissue = "Ankle"

rule all:
    input:
        # Use the actual output filenames produced by the consensus rules (merged_cond outputs)
        expand("output/peaks/{tissue}/merged_cond/{cond}_macs3_merged.bed",
               cond=list(condition_samples.keys()), tissue=tissue),
        expand("output/peaks/{tissue}/merged_cond/{cond}_macs2_merged.bed",
               cond=list(condition_samples.keys()), tissue=tissue),
        expand("output/peaks/{tissue}/merged_cond/{cond}_hmmratac_merged.bed",
               cond=list(condition_samples.keys()), tissue=tissue),
        # Merged peaks across all samples for each caller 
        expand("output/peaks/{tissue}/merged/{merged}.bed",
               tissue=tissue, merged=['allsamples_macs3_merged','allsamples_macs2_merged','allsamples_hmmratac_merged']),
        expand("output/peaks/{tissue}/merged/{merged}.saf",
               tissue=tissue, merged=['allsamples_macs3_merged','allsamples_macs2_merged','allsamples_hmmratac_merged']),
        expand("output/peaks/{tissue}/merged/{merged}_counts.txt",
               tissue=tissue, merged=['allsamples_macs3_merged','allsamples_macs2_merged','allsamples_hmmratac_merged']),
        # FRiP outputs
        expand("output/peaks/{tissue}/frip/{sampleName}_macs3_frip.txt",
               sampleName=list(bamfile.keys()), tissue=tissue),
        expand("output/peaks/{tissue}/frip/{sampleName}_macs2_frip.txt",
               sampleName=list(bamfile.keys()), tissue=tissue),
        expand("output/peaks/{tissue}/frip/{sampleName}_hmmratac_frip.txt",
               sampleName=list(bamfile.keys()), tissue=tissue)

##############################################
# MACS3 COUNTING RULES
##############################################

rule macs3_peak_count:
    input:
        macs3_peak = lambda wildcards: expand(
            "output/peaks/{tissue}/macs3/{sp_nm_cond}_peaks.narrowPeak",
            sp_nm_cond=condition_samples[wildcards.cond], tissue=tissue),
        bam = lambda wildcards: expand(
            "output/filtered/blk_filter/{sp_nm_cond}.sorted_final.bam",
            sp_nm_cond=condition_samples[wildcards.cond])
    output:
        macs3_bed = "output/peaks/{tissue}/merged_cond/{cond}_macs3_merged.bed",
        macs3_saf = "output/peaks/{tissue}/merged_cond/{cond}_macs3_merged.saf",
        count    = "output/peaks/{tissue}/merged_cond/{cond}_macs3_merged_counts.txt"
    log:
        err = "output/logs/{tissue}_{cond}_macs3_merge_peakcount.err",
        out = "output/logs/{tissue}_{cond}_macs3_merge_peakcount.out"
    threads: 8
    params:
        bedtoolsVer = config['bedtools'],
        subreadVer  = config['subread'],
        min_samples = 6
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/merged_cond

        multiIntersectBed -i {input.macs3_peak} | \
            awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
            awk '($3-$2) >= 50' | \
            bedtools merge > {output.macs3_bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
            {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.macs3_bed} > {output.macs3_saf}

        # Count reads in the region
        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {output.macs3_saf} \
            -o {output.count} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

rule macs3_merge_conditions:
    input:
        bed = expand("output/peaks/{tissue}/merged_cond/{cond}_macs3_merged.bed",
                     cond=['CTL', 'FNF'], tissue=tissue),
        bam = expand("output/filtered/blk_filter/{sampleName}.sorted_final.bam",
                     sampleName=bamfile.keys())
    output:
        bed    = "output/peaks/{tissue}/merged/allsamples_macs3_merged.bed",
        saf    = "output/peaks/{tissue}/merged/allsamples_macs3_merged.saf",
        counts = "output/peaks/{tissue}/merged/allsamples_macs3_merged_counts.txt"
    log:
        err = "output/logs/{tissue}_macs3_merged_counts_Allsamples.err",
        out = "output/logs/{tissue}_macs3_merged_counts_Allsamples.out"
    params:
        bedtoolsVer = config['bedtools'],
        subreadVer  = config['subread']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/merged

        cat {input.bed} | sort -k1,1 -k2,2n | bedtools merge > {output.bed} 2>> {log.err}

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
            {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {output.saf} \
            -o {output.counts} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

rule macs3_frip_calc:
    input:
        mergedSAF = rules.macs3_merge_conditions.output.saf,
        bam = lambda wildcards: "output/filtered/blk_filter/{}.sorted_final.bam".format(wildcards.sampleName)
    output:
        count = "output/peaks/{tissue}/frip/{sampleName}_macs3_frip.txt"
    log:
        err = "output/logs/{tissue}_{sampleName}_macs3_frip.err",
        out = "output/logs/{tissue}_{sampleName}_macs3_frip.out"
    threads: 4
    params:
        subreadVer = config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/frip

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {input.mergedSAF} \
            -o {output.count} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

##############################################
# MACS2 COUNTING RULES
##############################################

rule macs2_peak_count:
    input:
        peak = lambda wildcards: expand(
            "output/peaks/{tissue}/macs2/{sp_nm_cond}_peaks.narrowPeak",
            sp_nm_cond=condition_samples[wildcards.cond], tissue=tissue),
        bam = lambda wildcards: expand(
            "output/filtered/blk_filter/{sp_nm_cond}.sorted_final.bam",
            sp_nm_cond=condition_samples[wildcards.cond])
    output:
        bed   = "output/peaks/{tissue}/merged_cond/{cond}_macs2_merged.bed",
        saf   = "output/peaks/{tissue}/merged_cond/{cond}_macs2_merged.saf",
        count = "output/peaks/{tissue}/merged_cond/{cond}_macs2_merged_counts.txt"
    log:
        err = "output/logs/{tissue}_{cond}_macs2_merge_peakcount.err",
        out = "output/logs/{tissue}_{cond}_macs2_merge_peakcount.out"
    threads: 8
    params:
        bedtoolsVer = config['bedtools'],
        subreadVer  = config['subread'],
        min_samples = 6
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/merged_cond

        multiIntersectBed -i {input.peak} | \
            awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
            awk '($3-$2) >= 50' | \
            bedtools merge > {output.bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
            {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {output.saf} \
            -o {output.count} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

rule macs2_merge_conditions:
    input:
        bed = expand("output/peaks/{tissue}/merged_cond/{cond}_macs2_merged.bed", cond=['CTL', 'FNF'], tissue=tissue),
        bam = expand("output/filtered/blk_filter/{sampleName}.sorted_final.bam", sampleName=bamfile.keys())
    output:
        bed    = "output/peaks/{tissue}/merged/allsamples_macs2_merged.bed",
        saf    = "output/peaks/{tissue}/merged/allsamples_macs2_merged.saf",
        counts = "output/peaks/{tissue}/merged/allsamples_macs2_merged_counts.txt"
    log:
        err = "output/logs/{tissue}_macs2_merged_counts_Allsamples.err",
        out = "output/logs/{tissue}_macs2_merged_counts_Allsamples.out"
    params:
        bedtoolsVer = config['bedtools'],
        subreadVer  = config['subread']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/merged

        cat {input.bed} | sort -k1,1 -k2,2n | bedtools merge > {output.bed} 2>> {log.err}

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
            {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {output.saf} \
            -o {output.counts} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

rule macs2_frip_calc:
    input:
        mergedSAF = rules.macs2_merge_conditions.output.saf,
        bam = lambda wildcards: "output/filtered/blk_filter/{}.sorted_final.bam".format(wildcards.sampleName)
    output:
        count = "output/peaks/{tissue}/frip/{sampleName}_macs2_frip.txt"
    log:
        err = "output/logs/{tissue}_{sampleName}_macs2_frip.err",
        out = "output/logs/{tissue}_{sampleName}_macs2_frip.out"
    threads: 4
    params:
        subreadVer = config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/frip

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {input.mergedSAF} \
            -o {output.count} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

##############################################
# HMMRATAC COUNTING RULES
##############################################

rule hmmratac_peak_count:
    input:
        peak = lambda wildcards: expand(
            "output/peaks/{tissue}/hmmratac/{sp_nm_cond}_accessible_regions.narrowPeak",
            sp_nm_cond=condition_samples[wildcards.cond], tissue=tissue),
        bam = lambda wildcards: expand(
            "output/filtered/blk_filter/{sp_nm_cond}.sorted_final.bam",
            sp_nm_cond=condition_samples[wildcards.cond])
    output:
        bed   = "output/peaks/{tissue}/merged_cond/{cond}_hmmratac_merged.bed",
        saf   = "output/peaks/{tissue}/merged_cond/{cond}_hmmratac_merged.saf",
        count = "output/peaks/{tissue}/merged_cond/{cond}_hmmratac_merged_counts.txt"
    log:
        err = "output/logs/{tissue}_{cond}_hmmratac_merge_peakcount.err",
        out = "output/logs/{tissue}_{cond}_hmmratac_merge_peakcount.out"
    threads: 8
    params:
        bedtoolsVer = config['bedtools'],
        subreadVer  = config['subread'],
        min_samples = 6
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/merged_cond

        multiIntersectBed -i {input.peak} | \
            awk -v min_samples={params.min_samples} '$4 >= min_samples' | \
            awk '($3-$2) >= 40' | \
            bedtools merge > {output.bed} 2> {log.err}

        # Convert to SAF
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
            {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {output.saf} \
            -o {output.count} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

rule hmmratac_merge_conditions:
    input:
        bed = expand("output/peaks/{tissue}/merged_cond/{cond}_hmmratac_merged.bed", cond=['CTL', 'FNF'], tissue=tissue),
        bam = expand("output/filtered/blk_filter/{sampleName}.sorted_final.bam", sampleName=bamfile.keys())
    output:
        bed    = "output/peaks/{tissue}/merged/allsamples_hmmratac_merged.bed",
        saf    = "output/peaks/{tissue}/merged/allsamples_hmmratac_merged.saf",
        counts = "output/peaks/{tissue}/merged/allsamples_hmmratac_merged_counts.txt"
    log:
        err = "output/logs/{tissue}_hmmratac_merged_counts_Allsamples.err",
        out = "output/logs/{tissue}_hmmratac_merged_counts_Allsamples.out"
    params:
        bedtoolsVer = config['bedtools'],
        subreadVer  = config['subread']
    shell:
        """
        module load bedtools/{params.bedtoolsVer}
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/merged

        cat {input.bed} | sort -k1,1 -k2,2n | bedtools merge > {output.bed} 2>> {log.err}

        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} \
            {{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t."}}' {output.bed} > {output.saf}

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {output.saf} \
            -o {output.counts} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """

rule hmmratac_frip_calc:
    input:
        mergedSAF = rules.hmmratac_merge_conditions.output.saf,
        bam = lambda wildcards: "output/filtered/blk_filter/{}.sorted_final.bam".format(wildcards.sampleName)
    output:
        count = "output/peaks/{tissue}/frip/{sampleName}_hmmratac_frip.txt"
    log:
        err = "output/logs/{tissue}_{sampleName}_hmmratac_frip.err",
        out = "output/logs/{tissue}_{sampleName}_hmmratac_frip.out"
    threads: 4
    params:
        subreadVer = config['subread']
    shell:
        """
        module load subread/{params.subreadVer}
        mkdir -p output/peaks/{tissue}/frip

        featureCounts -p -F SAF \
            -T {threads} \
            --primary \
            -C \
            -a {input.mergedSAF} \
            -o {output.count} \
            {input.bam} 2>> {log.err} 1>> {log.out}
        """