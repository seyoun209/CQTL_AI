#!/usr/bin/env python3

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakefiles.utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
samples = samples.astype(str)

# Exclusion criteria
if config['exclude_options']['exclude_replicates']:
    samples = samples[~samples['Protocol_notes'].str.contains('replicate', case=False, na=False)]
if config['exclude_options']['exclude_synovium']:
    samples = samples[~samples['Tissue'].str.contains('synovium', case=False, na=False)]
if config['exclude_options']['exclude_femur']:
    samples = samples[~samples['Tissue'].str.contains('Femur', case=False, na=False)]

## Create merged names
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

# Separate samples by condition
ctl_samples = samples[samples['Condition'] == 'CTL']
fnf_samples = samples[samples['Condition'] == 'FNF']

# Create BAM and BAI file mappings
bam_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/output/filtered/{x.name}_filtered.bam").to_dict()
bai_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/output/filtered/{x.name}_filtered.bam.bai").to_dict()

rule all:
    input:
        expand("output/wasp/filter_remapped_reads/{sampleName}_filtered.{ext}",
               sampleName=bam_files.keys(),
               ext=['remap.fq1.gz', 'remap.fq2.gz','to.remap.bam','keep.bam']), 
        expand("output/wasp/filtered_bam/{sampleName}.keep.bam",
                sampleName=bam_files.keys()),
        expand("output/wasp/rmdup/{sampleName}{ext}",
                sampleName=bam_files.keys(),
                ext=['.rmdup.bam','.sorted_wasp.bam','.sorted_wasp.bam.bai','_stats.txt'])



rule extract_sample_ids:
    input:
        pbs_matched="output/geno/pbs_geno/pbs_matched.txt",
        fnf_matched="output/geno/fnf_geno/fnf_matched.txt"
    output:
        sample_ids_pbs="output/wasp/pbs_sample_ids.txt",
        sample_ids_fnf="output/wasp/fnf_sample_ids.txt"
    log:
        pbs_log="output/logs/extract_sample_ids_pbs.log",
        fnf_log="output/logs/extract_sample_ids_fnf.log"
    shell:
        """

        # Extract the first column for PBS sample IDs
        cut -d' ' -f1 {input.pbs_matched} > {output.sample_ids_pbs} 2> {log.pbs_log}

        # Extract the first column for FNF sample IDs
        cut -d' ' -f1 {input.fnf_matched} > {output.sample_ids_fnf} 2> {log.fnf_log}
        """

rule snp2h5:
    input:
        pbs_vcf="output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz",
        fnf_vcf="output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz",
        sample_ids_pbs=rules.extract_sample_ids.output.sample_ids_pbs,
        sample_ids_fnf=rules.extract_sample_ids.output.sample_ids_fnf
    output:
        pbs_haplotypes="output/wasp/pbs/haplotypes.h5",
        pbs_snp_index="output/wasp/pbs/snp_index.h5",
        pbs_snp_tab="output/wasp/pbs/snp_tab.h5",
        fnf_haplotypes="output/wasp/fnf/haplotypes.h5",
        fnf_snp_index="output/wasp/fnf/snp_index.h5",
        fnf_snp_tab="output/wasp/fnf/snp_tab.h5"
    params:
        wasp_path=config['wasp_dir'],
        chromsize=config['chrom_size'],
        wasp_ver=config['waspVer']
    log:
        pbs_err="output/logs/pbs_snp2h5.err",
        fnf_err="output/logs/fnf_snp2h5.err"
    shell:
        """
        module add wasp/{params.wasp_ver}
        mkdir -p output/wasp/pbs output/wasp/fnf

        # Process PBS VCF
        {params.wasp_path}/snp2h5/snp2h5 \
        --chrom {params.chromsize} \
        --format vcf \
        --haplotype {output.pbs_haplotypes} \
        --snp_index {output.pbs_snp_index} \
        --snp_tab {output.pbs_snp_tab} \
        {input.pbs_vcf} 2> {log.pbs_err}

        # Process FNF VCF
        {params.wasp_path}/snp2h5/snp2h5 \
        --chrom {params.chromsize} \
        --format vcf \
        --haplotype {output.fnf_haplotypes} \
        --snp_index {output.fnf_snp_index} \
        --snp_tab {output.fnf_snp_tab} \
        {input.fnf_vcf} 2> {log.fnf_err}
        """

rule find_intersecting_snps:
    input:
        bam=lambda wildcards: bam_files[wildcards.sampleName],
        pbs_snp_tab=rules.snp2h5.output.pbs_snp_tab,
        pbs_snp_index=rules.snp2h5.output.pbs_snp_index,
        pbs_haplotypes=rules.snp2h5.output.pbs_haplotypes,
        fnf_snp_tab=rules.snp2h5.output.fnf_snp_tab,
        fnf_snp_index=rules.snp2h5.output.fnf_snp_index,
        fnf_haplotypes=rules.snp2h5.output.fnf_haplotypes,
        sample_ids_pbs=rules.extract_sample_ids.output.sample_ids_pbs,
        sample_ids_fnf=rules.extract_sample_ids.output.sample_ids_fnf
    output:
        pbs_fq1="output/wasp/filter_remapped_reads/{sampleName}_filtered.remap.fq1.gz",
        pbs_fq2="output/wasp/filter_remapped_reads/{sampleName}_filtered.remap.fq2.gz",
        to_remap_bam="output/wasp/filter_remapped_reads/{sampleName}_filtered.to.remap.bam",
        keep_bam="output/wasp/filter_remapped_reads/{sampleName}_filtered.keep.bam"
    threads: 4
    params:
        wasp_path=config['wasp_dir'],
        wasp_ver=config['waspVer']
    log:
        err="output/logs/find_intersecting_snps_pbs_{sampleName}.err"
    shell:
        """
        module add wasp/{params.wasp_ver}

        mkdir -p output/wasp/filter_remapped_reads

        if [[ "{wildcards.sampleName}" == *"CTL"* ]]; then
            snp_tab={input.pbs_snp_tab}
            snp_index={input.pbs_snp_index}
            haplotypes={input.pbs_haplotypes}
            sample_id={input.sample_ids_pbs}
        elif [[ "{wildcards.sampleName}" == *"FNF"* ]]; then
            snp_tab={input.fnf_snp_tab}
            snp_index={input.fnf_snp_index}
            haplotypes={input.fnf_haplotypes}
            sample_id={input.sample_ids_fnf}
        else
            echo "Error: Sample type not recognized in {wildcards.sampleName}"
            exit 1
        fi

        # Add error checking here
        if [[ -z "$snp_tab" ]]; then
            echo "Error: snp_tab not set"
            exit 1
        fi
        if [[ -z "$snp_index" ]]; then
            echo "Error: snp_index not set"
            exit 1
        fi
        if [[ -z "$haplotypes" ]]; then
            echo "Error: haplotypes not set"
            exit 1
        fi
        if [[ -z "$sample_id" ]]; then
            echo "Error: sample_id not set"
            exit 1
        fi 

        # Then run WASP command
        python {params.wasp_path}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir output/wasp/filter_remapped_reads \
        --snp_tab $snp_tab \
        --snp_index $snp_index \
        --haplotype $haplotypes \
        --samples $sample_id \
        {input.bam} 2> {log.err}

        """


rule remap_wasp_reads:
    input:
        fq1=rules.find_intersecting_snps.output.pbs_fq1,
        fq2=rules.find_intersecting_snps.output.pbs_fq2,
    output:
        sorted_bam="output/wasp/remapped/{sampleName}_remapped_sorted.bam",
        bai="output/wasp/remapped/{sampleName}_remapped_sorted.bam.bai",
        stats="output/wasp/remapped/{sampleName}_stats.txt" 
    params:
        bwa_index=config['bwa_index'],
        bwa_version=config['bwaVers'],
        samtools_version=config['samtoolsVers']
    threads: 8
    log:
        err="output/logs/remap_wasp_{sampleName}.err",
        out="output/logs/remap_wasp_{sampleName}.out"
    shell:
        """
        module load bwa/{params.bwa_version}
        module load samtools/{params.samtools_version}
        mkdir -p output/wasp/remapped
        
        # Remap with BWA and pipe to sorting
        bwa mem -t {threads} -M {params.bwa_index} {input.fq1} {input.fq2} | \
        samtools view -b | \
        samtools sort -o {output.sorted_bam} \
        1> {log.out} 2> {log.err}

        samtools flagstat {output.sorted_bam} > {output.stats} 2>> {log.err}

        # Index the sorted BAM
        samtools index {output.sorted_bam} 1>> {log.out} 2>> {log.err}
        """

rule filter_remapped_reads:
    input:
        to_remap_bam=rules.find_intersecting_snps.output.to_remap_bam,
        remapped_bam=rules.remap_wasp_reads.output.sorted_bam,
        remapped_bai=rules.remap_wasp_reads.output.bai
    output:
        keep_bam="output/wasp/filtered_bam/{sampleName}.keep.bam"
    params:
        wasp_path=config['wasp_dir'],
        wasp_version=config['waspVer']
    log:
        "output/logs/filter_remapped_{sampleName}.log"
    shell:
        """
        module add wasp/{params.wasp_version}
        
        mkdir -p output/wasp/filtered_bam
        
        python {params.wasp_path}/mapping/filter_remapped_reads.py \
            {input.to_remap_bam} \
            {input.remapped_bam} \
            {output.keep_bam} \
            2> {log}
        """

rule merge:
    input:
        keepbam=rules.filter_remapped_reads.output.keep_bam,
        interbam=rules.find_intersecting_snps.output.keep_bam
    output:
        merged_bam="output/wasp/merge/{sampleName}.keep.merge.bam",
        sort_merged_bam="output/wasp/merge/{sampleName}.sort_keep.merge.bam",
        sort_merged_bai="output/wasp/merge/{sampleName}.sort_keep.merge.bam.bai"
    params:
        samtools_version=config['samtoolsVers']
    threads: 4
    log:
        "output/logs/merge_{sampleName}.log"
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p output/wasp/merge

        samtools merge {output.merged_bam} {input.keepbam} {input.interbam}
        samtools sort -o {output.sort_merged_bam} {output.merged_bam}
        samtools index {output.sort_merged_bam}

        """

rule rmDupWASP:
    input:
        merged_bam=rules.merge.output.sort_merged_bam
    output:
        rmdup_bam="output/wasp/rmdup/{sampleName}.rmdup.bam",
        sort_rmdup_bam="output/wasp/rmdup/{sampleName}.sorted_wasp.bam",
        index_bam="output/wasp/rmdup/{sampleName}.sorted_wasp.bam.bai",
        stats="output/wasp/rmdup/{sampleName}_stats.txt"
    params:
        pythonVer=config['python'],
        wasp_path=config['wasp_dir'],
        wasp_version=config['waspVer'],
        samtools_version=config['samtoolsVers']
    threads:8
    log:
        err="output/logs/rmdup_wasp_{sampleName}.err",
        out="output/logs/rmdup_wasp_{sampleName}.out"
    shell:
        """
        module add wasp/{params.wasp_version}
        module load samtools/{params.samtools_version}
        module load python/{params.pythonVer}
        mkdir -p output/wasp/rmdup

        python {params.wasp_path}/mapping/rmdup_pe.py {input.merged_bam} {output.rmdup_bam}

        samtools sort -o {output.sort_rmdup_bam} {output.rmdup_bam} 1> {log.out} 2> {log.err}
        samtools index {output.sort_rmdup_bam} 1>> {log.out} 2>> {log.err}
        samtools flagstat {output.sort_rmdup_bam} > {output.stats} 2>> {log.err}

        """

