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

# Add this after reading the samplesheet (exclude if it's needed
#if config['exclude_options']['exclude_replicates']:
#    samples = samples[~samples['Protocol_notes'].str.contains('replicate', case=False, na=False)]
if config['exclude_options']['exclude_synovium']:
    samples = samples[~samples['Tissue'].str.contains('synovium', case=False, na=False)]
if config['exclude_options']['exclude_femur']:
    samples = samples[~samples['Tissue'].str.contains('Femur', case=False, na=False)]


## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ['Proj','Donor','Condition','Tissue','Protocol_notes']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

# Create BAM and BAI file mappings with base directory
bam_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/rna_output/align/{x.name}_RG.bam").to_dict()
bai_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/rna_output/align/{x.name}_RG.bam.bai").to_dict()

## Define actions on success
onsuccess:
        ## Success message
        print("VerifybamID  completed successfully! Wahoo!")

rule all:
    input:
        expand("rna_output/QC/{sampleName}_verifyBamID.{ext}",sampleName=bam_files.keys(),ext=['selfSM','selfRG','bestRG','bestSM','depthRG','depthSM'])

rule verifybamid:
    input:
        bam = lambda wildcards: bam_files.get(wildcards.sampleName),
        bai = lambda wildcards: bai_files.get(wildcards.sampleName),
        pbs_vcf = "output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz",
        fnf_vcf = "output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz"
    output:
        selfSM = 'rna_output/QC/{sampleName}_verifyBamID.selfSM',
        selfRG = 'rna_output/QC/{sampleName}_verifyBamID.selfRG',
        bestRG = 'rna_output/QC/{sampleName}_verifyBamID.bestRG',
        bestSM = 'rna_output/QC/{sampleName}_verifyBamID.bestSM',
        depthRG = 'rna_output/QC/{sampleName}_verifyBamID.depthRG',
        depthSM = 'rna_output/QC/{sampleName}_verifyBamID.depthSM'
    params:
        verifybamid = config['verifybamid'],
        out_prefix = "rna_output/QC/{sampleName}_verifyBamID"
    log:
        err = "rna_output/logs/verifyBamID_{sampleName}.err"
    shell:
        """

        if [[ "{wildcards.sampleName}" == *"CTL"* ]]; then
            vcf_file={input.pbs_vcf}
        elif [[ "{wildcards.sampleName}" == *"FNF"* ]]; then
            vcf_file={input.fnf_vcf}
        else
            echo "Error: Sample type not recognized in {wildcards.sampleName}"
            exit 1
        fi

        echo "$vcf_file"

        {params.verifybamid} --vcf $vcf_file --bam {input.bam} --bai {input.bai} --best --out {params.out_prefix} 2> {log.err}
        """

