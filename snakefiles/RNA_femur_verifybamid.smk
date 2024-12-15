#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Define the list of all possible sample names and their corresponding VCF files
SAMPLE_VCF_MAP = {
    "CQTL_AM7778FE_FNF_Femur_1": "output/geno/fnf_geno_other/fnf.rename_wCHR.vcf.gz",
    "CQTL_AM7778FE_CTL_Femur_1": "output/geno/pbs_geno_other/pbs.rename_wCHR.vcf.gz"
}

# Print samples to be processed (optional, for verification)
print("Samples to be processed:", list(SAMPLE_VCF_MAP.keys()))

rule all:
    input:
        expand("rna_output/QC/{sampleName}_verifyBamID.{ext}",
               sampleName=SAMPLE_VCF_MAP.keys(),
               ext=['selfSM', 'selfRG', 'bestRG', 'bestSM', 'depthRG', 'depthSM'])

rule verifybamid:
    input:
        bam = "rna_output/align/{sampleName}_RG.bam",
        bai = "rna_output/align/{sampleName}_RG.bam.bai",
        vcf = lambda wildcards: SAMPLE_VCF_MAP[wildcards.sampleName]
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
    wildcard_constraints:
        sampleName = "|".join(SAMPLE_VCF_MAP.keys())
    log:
        err = "rna_output/logs/verifyBamID_{sampleName}.err"
    shell:
        """
        echo "Using VCF file: {input.vcf}"

        {params.verifybamid} --vcf {input.vcf} --bam {input.bam} --bai {input.bai} --best --out {params.out_prefix} 2> {log.err}
        """

onsuccess:
    print("VerifybamID completed successfully! Wahoo!")
