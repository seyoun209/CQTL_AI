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
if config['exclude_options']['exclude_replicates']:
    samples = samples[~samples['Protocol_notes'].str.contains('replicate', case=False, na=False)]
if config['exclude_options']['exclude_synovium']:
    samples = samples[~samples['Tissue'].str.contains('synovium', case=False, na=False)]
if config['exclude_options']['exclude_femur']:
    samples = samples[~samples['Tissue'].str.contains('Femur', case=False, na=False)]


## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ['Proj','Donor','Condition','Tissue','Protocol_notes']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

# Create BAM and BAI file mappings with base directory
bam_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/output/align/{x.name}_RG.bam").to_dict()
bai_files = samples.groupby('mn').apply(lambda x: f"{config['base_dir']}/output/align/{x.name}_RG.bam.bai").to_dict()

## Define actions on success
onsuccess:
        ## Success message
        print("VCF_preprocess  completed successfully! Wahoo!")

rule all:
    input:
        "output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz",
        "output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz",
        expand("output/QC/{sampleName}_verifyBamID.{ext}",sampleName=bam_files.keys(),ext=['selfSM','selfRG','bestRG','bestSM','depthRG','depthSM'])
      
    

rule categorize_samples:
    output:
        ctl="output/CTL_samples.txt",
        fnf="output/FNF_samples.txt"
    log:
        out="output/logs/categorize_samples.out",
        err="output/logs/categorize_samples.err"
    shell:
        """
        mkdir -p output/logs
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/categorize_samples.sh 1> {log.out} 2> {log.err}
        
        # Move output files to expected locations
        mv /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/CTL_samples.txt {output.ctl} 2>> {log.err}
        mv /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/FNF_samples.txt {output.fnf} 2>> {log.err}
        """

rule process_vcf:
    input:
        PBS=rules.categorize_samples.output.ctl,
        FNF=rules.categorize_samples.output.fnf
    output:
        vcf_filtered="output/geno/cQTL_21_filtered.vcf.gz",
        pbs_matched="output/geno/pbs_geno/pbs_matched.txt",
        fnf_matched="output/geno/fnf_geno/fnf_matched.txt"
    params:
        samtoolsVer = config['samtoolsVers'],
        raw_vcf = config['vcf']
    benchmark:
        "output/benchmarks/process_vcf.txt"
    log:
        out="output/logs/process_vcf.out",
        err="output/logs/process_vcf.err"
    shell:
        """
        mkdir -p output/geno
        mkdir -p output/geno/pbs_geno output/geno/fnf_geno
        mkdir -p output/benchmarks
        mkdir -p output/logs
        
        # Copy sample files
        cp {input.PBS} {output.pbs_matched} 2>> {log.err}
        cp {input.FNF} {output.fnf_matched} 2>> {log.err}
        
        # Process VCF
        ml samtools/{params.samtoolsVer}
        
        cp {params.raw_vcf} {output.vcf_filtered} 2>> {log.err}
        tabix -p vcf {output.vcf_filtered} 2>> {log.err}

        """

rule reheader:
    input:
        VCF_in=rules.process_vcf.output.vcf_filtered,
        matched_pbs=rules.process_vcf.output.pbs_matched,
        matched_fnf=rules.process_vcf.output.fnf_matched
    output:
        pbs_vcf="output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz",
        fnf_vcf="output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz",
        pbs_temp="output/geno/pbs_geno/pbs.rename.vcf.gz",
        fnf_temp="output/geno/fnf_geno/fnf.rename.vcf.gz"
    params:
       samtoolsVer = config['samtoolsVers']
    log:
        pbsout="output/logs/PBS_rename.out",
        pbserr="output/logs/PBS_rename.err",
        fnfout="output/logs/FNF_rename.out",
        fnferr="output/logs/FNF_rename.err"
    shell:
        """
        mkdir -p output/logs
        
       sh /work/users/s/e/seyoun/CQTL_AI/scripts/rename.sh {input.VCF_in} {input.matched_pbs} {output.pbs_temp} {output.pbs_vcf} 1> {log.pbsout} 2> {log.pbserr}
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/rename.sh {input.VCF_in} {input.matched_fnf} {output.fnf_temp} {output.fnf_vcf} 1> {log.fnfout} 2> {log.fnferr}
        """

rule verifybamid:
    input:
        bam = lambda wildcards: bam_files.get(wildcards.sampleName),
        bai = lambda wildcards: bai_files.get(wildcards.sampleName),
        pbs_vcf = rules.reheader.output.pbs_vcf,
        fnf_vcf = rules.reheader.output.fnf_vcf
    output:
        selfSM = 'output/QC/{sampleName}_verifyBamID.selfSM',
        selfRG = 'output/QC/{sampleName}_verifyBamID.selfRG',
        bestRG = 'output/QC/{sampleName}_verifyBamID.bestRG',
        bestSM = 'output/QC/{sampleName}_verifyBamID.bestSM',
        depthRG = 'output/QC/{sampleName}_verifyBamID.depthRG',
        depthSM = 'output/QC/{sampleName}_verifyBamID.depthSM'
    params:
        verifybamid = config['verifybamid'],
        out_prefix = "output/QC/{sampleName}_verifyBamID"
    log:
        err = "output/logs/verifyBamID_{sampleName}.err"
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

