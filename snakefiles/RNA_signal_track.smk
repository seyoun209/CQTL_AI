#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakemake.utils.namer import namer

# Read in samplesheet (use bam)
samples = pd.read_csv(config["bam_samplesheet"], sep='\t')
#samples = pd.read_csv("rna_samplesheet.txt",sep='\t')
samples = samples.astype(str)


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
#tissue = samples['Tissue'].iloc[0]
tissues = samples['Tissue'].unique().tolist()

ankle_samples = [s for s in bamfile.keys() if "_Ankle_" in s]
femur_samples = [s for s in bamfile.keys() if "_Femur_" in s]
synovium_samples = [s for s in bamfile.keys() if "_synovium_" in s]

# Create condition groups for each tissue
ankle_ctl = [s for s in ankle_samples if "_CTL_" in s]
ankle_fnf = [s for s in ankle_samples if "_FNF_" in s]

femur_ctl = [s for s in femur_samples if "_CTL_" in s]
femur_fnf = [s for s in femur_samples if "_FNF_" in s]

synovium_ctl = [s for s in synovium_samples if "_CTL_" in s]
synovium_fnf = [s for s in synovium_samples if "_FNF_" in s]

# Create mergeSample dictionary
mergeSample = {
    "Ankle_CTL": ankle_ctl,
    "Ankle_FNF": ankle_fnf,
    "Femur_CTL": femur_ctl,
    "Femur_FNF": femur_fnf,
    "synovium_CTL": synovium_ctl,
    "synovium_FNF": synovium_fnf
}

##### Define rules #####

rule all:
    input:
        # Signal files
        expand("rna_output/signals/Ankle/{sampleName}.bw", sampleName=ankle_samples),
        expand("rna_output/signals/Ankle/norm/{sampleName}.bw", sampleName=ankle_samples),
        expand("rna_output/signals/Femur/{sampleName}.bw", sampleName=femur_samples),
        expand("rna_output/signals/Femur/norm/{sampleName}.bw", sampleName=femur_samples),
        expand("rna_output/signals/synovium/{sampleName}.bw", sampleName=synovium_samples),
        expand("rna_output/signals/synovium/norm/{sampleName}.bw", sampleName=synovium_samples),
        
        # Merged files for Ankle
        expand("rna_output/signals/Ankle/merged_signal/{cond}_{ext}",
               cond=["CTL", "FNF"],
               ext=['sorted.bw', 'norm_sorted.bw']),
        expand("rna_output/signals/Ankle/strand/{cond}_{ext}",
               cond=["CTL", "FNF"],
               ext=['fwd.bw','rev.bw','norm_fwd.bw','norm_rev.bw']),
               
        # Merged files for Femur
        expand("rna_output/signals/Femur/merged_signal/{cond}_{ext}",
               cond=["CTL", "FNF"],
               ext=['sorted.bw', 'norm_sorted.bw']),
        expand("rna_output/signals/Femur/strand/{cond}_{ext}",
               cond=["CTL", "FNF"],
               ext=['fwd.bw','rev.bw','norm_fwd.bw','norm_rev.bw']),
               
        # Merged files for synovium
        expand("rna_output/signals/synovium/merged_signal/{cond}_{ext}",
               cond=["CTL", "FNF"],
               ext=['sorted.bw', 'norm_sorted.bw']),
        expand("rna_output/signals/synovium/strand/{cond}_{ext}",
               cond=["CTL", "FNF"],
               ext=['fwd.bw','rev.bw','norm_fwd.bw','norm_rev.bw'])

rule signal:
    input:
        bam= lambda wildcards: bamfile.get(wildcards.sampleName)[0]
    output:
        signal = "rna_output/signals/{tissue}/{sampleName}.bw",
        norm_signal = "rna_output/signals/{tissue}/norm/{sampleName}.bw"
    log:
        err = 'rna_output/logs/signal_{tissue}_{sampleName}.err',
        err_norm = 'rna_output/logs/Normsignal_{tissue}_{sampleName}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/{wildcards.tissue}
        mkdir -p rna_output/signals/{wildcards.tissue}/norm

        bamCoverage \
            --bam {input.bam} \
            -o {output.signal} \
            > {log.err} 2>&1

        bamCoverage --normalizeUsing {params.NormOption} \
                --effectiveGenomeSize {params.effective_genomeSize} \
                --bam {input.bam} -o {output.norm_signal} > {log.err_norm} 2>&1
        """
rule mergeAlign:
    input:
        bam_files = lambda wildcards: [bamfile[sample][0] for sample in mergeSample[f"{wildcards.tissue}_{wildcards.cond}"]]
    output:
        bam = "rna_output/signals/{tissue}/mergedAlign/{cond}_sorted.bam",
        bai = "rna_output/signals/{tissue}/mergedAlign/{cond}_sorted.bam.bai",
        stats = "rna_output/signals/{tissue}/mergedAlign/{cond}_stats.txt"
    log:
        err='rna_output/logs/mergedAligned_{tissue}_{cond}.err'
    params:
        samtools_version = config['samtoolsVers']
    shell:
        """
        module load samtools/{params.samtools_version}
        mkdir -p rna_output/signals/{wildcards.tissue}/mergedAlign

        samtools merge {output.bam} {input.bam_files} >> {log.err} 2>&1
        samtools flagstat {output.bam} > {output.stats} >> {log.err} 2>&1
        samtools index {output.bam} >> {log.err} 2>&1
        """


rule mergeSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal='rna_output/signals/{tissue}/merged_signal/{cond}_sorted.bw',
        norm_signal='rna_output/signals/{tissue}/merged_signal/{cond}_norm_sorted.bw'
    log:
        err='rna_output/logs/mergedSignal_{tissue}_{cond}.err',
        err_norm='rna_output/logs/Norm_mergedSignal_{tissue}_{cond}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/{wildcards.tissue}/merged_signal

        bamCoverage \
            --bam {input.bam} \
            -o {output.signal} \
            > {log.err} 2>&1

        bamCoverage --normalizeUsing {params.NormOption} \
                --effectiveGenomeSize {params.effective_genomeSize} \
                --bam {input.bam} -o {output.norm_signal} > {log.err_norm} 2>&1
        """

rule mergeForwardSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal_fw = "rna_output/signals/{tissue}/strand/{cond}_fwd.bw",
        norm_signal_fw = "rna_output/signals/{tissue}/strand/{cond}_norm_fwd.bw"
    log:
        err = 'rna_output/logs/mergeSignal_forward_{tissue}_{cond}.err',
        err_norm = 'rna_output/logs/Norm_mergeSignal_forward_{tissue}_{cond}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/{wildcards.tissue}/strand

        bamCoverage --filterRNAstrand forward -b {input.bam} -o {output.signal_fw}  > {log.err} 2>&1
        
        bamCoverage --filterRNAstrand forward \
                --normalizeUsing {params.NormOption} \
                --effectiveGenomeSize {params.effective_genomeSize} \
                --bam {input.bam} -o {output.norm_signal_fw} > {log.err_norm} 2>&1

        """

rule mergeReverseSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal_rv = "rna_output/signals/{tissue}/strand/{cond}_rev.bw",
        norm_signal_rv = "rna_output/signals/{tissue}/strand/{cond}_norm_rev.bw"
    log:
        err = 'rna_output/logs/mergeSignal_reverse_{tissue}_{cond}.err',
        err_norm = 'rna_output/logs/Norm_mergeSignal_reverse_{tissue}_{cond}.err'
    params:
        deeptools_ver=config['deeptoolsVers'],
        effective_genomeSize=config['effective_genome_size'],
        NormOption=config['normalize_option']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p rna_output/signals/{wildcards.tissue}/strand

        bamCoverage --filterRNAstrand reverse -b {input.bam} -o {output.signal_rv}  > {log.err} 2>&1
        
        bamCoverage --filterRNAstrand reverse \
                --normalizeUsing {params.NormOption} \
                --effectiveGenomeSize {params.effective_genomeSize} \
                --bam {input.bam} -o {output.norm_signal_rv} > {log.err_norm} 2>&1

        """


