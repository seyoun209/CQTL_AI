#!/usr/bin/env python3

# Include required preprocessing files
include: '/work/users/s/e/seyoun/CQTL_sQTL/snakefiles/VCFpreprocess.smk'
include: '/work/users/s/e/seyoun/CQTL_sQTL/snakefiles/ATAC_preprocess.smk'

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
        cut -f 1 {input.pbs_matched} > {output.sample_ids_pbs} 2> {log.pbs_log}

        # Extract the first column for FNF sample IDs
        cut -f 1 {input.fnf_matched} > {output.sample_ids_fnf} 2> {log.fnf_log}
        """

# Rule: Convert VCF to HDF5 for WASP
rule snp2h5:
    input:
        pbs_vcf=rules.reheader.output.pbs_vcf,
        fnf_vcf=rules.reheader.output.fnf_vcf,
        sample_ids_pbs=rules.extract_sample_ids.output.sample_ids_pbs,
        sample_ids_fnf=rules.extract_sample_ids.output.sample_ids_fnf
    output:
        # Output files for PBS and FNF
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
        bam="output/align/{sampleName}_sorted.bam",
        pbs_snp_tab="output/wasp/pbs/snp_tab.h5",
        pbs_snp_index="output/wasp/pbs/snp_index.h5",
        pbs_haplotypes="output/wasp/pbs/haplotypes.h5",
        fnf_snp_tab="output/wasp/fnf/snp_tab.h5",
        fnf_snp_index="output/wasp/fnf/snp_index.h5",
        fnf_haplotypes="output/wasp/fnf/haplotypes.h5",
        sample_ids_pbs="output/wasp/pbs_sample_ids.txt",
        sample_ids_fnf="output/wasp/fnf_sample_ids.txt"
    output:
        pbs_filtered_bam="output/wasp/pbs/{sampleName}_filtered.bam",
        fnf_filtered_bam="output/wasp/fnf/{sampleName}_filtered.bam"
    params:
        wasp_path=config['wasp_dir']
    log:
        pbs_err="output/logs/find_intersecting_snps_pbs_{sampleName}.err",
        fnf_err="output/logs/find_intersecting_snps_fnf_{sampleName}.err"
    shell:
        """
        module add wasp

        # PBS SNP Intersection
        python {params.wasp_path}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir output/wasp/pbs/{wildcards.sampleName}/filter_remapped_reads \
        --snp_tab {input.pbs_snp_tab} \
        --snp_index {input.pbs_snp_index} \
        --haplotype {input.pbs_haplotypes} \
        --samples {input.sample_ids_pbs} \
        {input.bam} 2> {log.pbs_err}

        # FNF SNP Intersection
        python {params.wasp_path}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir output/wasp/fnf/{wildcards.sampleName}/filter_remapped_reads \
        --snp_tab {input.fnf_snp_tab} \
        --snp_index {input.fnf_snp_index} \
        --haplotype {input.fnf_haplotypes} \
        --samples {input.sample_ids_fnf} \
        {input.bam} 2> {log.fnf_err}
        """

