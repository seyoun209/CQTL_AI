#!/usr/bin/env python3

rule all:
    input:
        "output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz",
        "output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz"

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
