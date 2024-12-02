#!/usr/bin/env python3

include: '/work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk'

# Define the list of all possible sample names
all_samples = [
    "CQTL_AM7778FE_FNF_Femur_1",
    "CQTL_AM7754_FNF_Ankle_replicate",
    "CQTL_AM7755_FNF_Ankle_replicate",
    "CQTL_AM7778FE_CTL_Femur_1",
    "CQTL_AM7754_CTL_Ankle_replicate",
    "CQTL_AM7755_CTL_Ankle_replicate"
    #"OTHER_SAMPLE_1",
    #"OTHER_SAMPLE_2"
]

# Filter samples to include only those with 'replicate' or 'Femur' in their names
selected_samples = [s for s in all_samples if "replicate" in s or "Femur" in s]

rule all_subset:
    input:
        "output/geno/pbs_geno_other/pbs.rename_wCHR.vcf.gz",
        "output/geno/fnf_geno_other/fnf.rename_wCHR.vcf.gz",
        "output/geno/pbs_geno_other/pbs.rename.vcf.gz",
        "output/geno/fnf_geno_other/fnf.rename.vcf.gz",
        expand("output/QC/{sampleName}_verifyBamID.{ext}", sampleName=selected_samples,
            ext=['selfSM','selfRG','bestRG','bestSM','depthRG','depthSM'])


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

rule verifybamid:
    input:
        bam = rules.add_read_groups.output.bam,
        bai = rules.add_read_groups.output.bai,
        pbs_vcf = rules.reheader_subset.output.pbs_vcf,
        fnf_vcf = rules.reheader_subset.output.fnf_vcf
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
    wildcard_constraints:
        sampleName = "|".join([
            "CQTL_AM7778FE_FNF_Femur_1",
            "CQTL_AM7754_FNF_Ankle_replicate", 
            "CQTL_AM7755_FNF_Ankle_replicate",
            "CQTL_AM7778FE_CTL_Femur_1",
            "CQTL_AM7754_CTL_Ankle_replicate",
            "CQTL_AM7755_CTL_Ankle_replicate"
        ])
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

