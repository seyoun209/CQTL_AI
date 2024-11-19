#!/usr/bin/env python3

include: '/work/users/s/e/seyoun/CQTL_sQTL/snakefiles/ATAC_preprocess.smk'

rule all:
    input:
        "output/geno/pbs_geno_other/pbs.rename_wCHR.vcf.gz",
        "output/geno/fnf_geno_other/fnf.rename_wCHR.vcf.gz",
        "output/geno/pbs_geno_other/pbs.rename.vcf.gz",
        "output/geno/fnf_geno_other/fnf.rename.vcf.gz",
        expand("output/align/{sampleName}_RG.bam", sampleName=["CQTL_AM7778FE_FNF_Femur_1",
                                                              "CQTL_AM7754_FNF_Ankle_replicate",
                                                              "CQTL_AM7755_FNF_Ankle_replicate",
                                                              "CQTL_AM7778FE_CTL_Femur_1",
                                                              "CQTL_AM7754_CTL_Ankle_replicate",
                                                              "CQTL_AM7755_CTL_Ankle_replicate"]),
        expand("output/ataqv/{sampleName}.ataqv.json", sampleName=["CQTL_AM7778FE_FNF_Femur_1",
                                                                   "CQTL_AM7754_FNF_Ankle_replicate",
                                                                   "CQTL_AM7755_FNF_Ankle_replicate",
                                                                   "CQTL_AM7778FE_CTL_Femur_1",
                                                                   "CQTL_AM7754_CTL_Ankle_replicate",
                                                                   "CQTL_AM7755_CTL_Ankle_replicate"])


rule categorize_femur_replicate_samples:
    output:
        ctl_femur_replicate="output/geno/pbs_geno_other/CTL_samples_with_replacements.txt",
        fnf_femur_replicate="output/geno/fnf_geno_other/FNF_samples_with_replacements.txt"
    log:
        out="output/logs/categorize_femur_replicate_samples.out",
        err="output/logs/categorize_femur_replicate_samples.err"
    shell:
        """
        mkdir -p output/logs
        mkdir -p output/geno/pbs_geno_other
        mkdir -p output/geno/fnf_geno_other
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/categorize_samplesv3.sh 1> {log.out} 2> {log.err}
        
        # Move output files to expected locations
        mv /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/CTL_samples_with_replacements.txt {output.ctl_femur_replicate} 2>> {log.err}
        mv /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/FNF_samples_with_replacements.txt {output.fnf_femur_replicate} 2>> {log.err}
        """


rule reheader_subset:
    input:
        matched_pbs=rules.categorize_femur_replicate_samples.output.ctl_femur_replicate,
        matched_fnf=rules.categorize_femur_replicate_samples.output.fnf_femur_replicate
    output:
        pbs_vcf="output/geno/pbs_geno_other/pbs.rename_wCHR.vcf.gz",
        fnf_vcf="output/geno/fnf_geno_other/fnf.rename_wCHR.vcf.gz",
        pbs_temp="output/geno/pbs_geno_other/pbs.rename.vcf.gz",
        fnf_temp="output/geno/fnf_geno_other/fnf.rename.vcf.gz"
    params:
        vcf_copied=config['vcf_copy']
    log:
        pbsout="output/logs/PBS_rename_other.out",
        pbserr="output/logs/PBS_rename_other.err",
        fnfout="output/logs/FNF_rename_other.out",
        fnferr="output/logs/FNF_rename_other.err"
    shell:
        """
        
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/rename.sh {params.vcf_copied} {input.matched_pbs} {output.pbs_temp} {output.pbs_vcf} 1> {log.pbsout} 2> {log.pbserr}
        sh /work/users/s/e/seyoun/CQTL_AI/scripts/rename.sh {params.vcf_copied} {input.matched_fnf} {output.fnf_temp} {output.fnf_vcf} 1> {log.fnfout} 2> {log.fnferr}
        """

