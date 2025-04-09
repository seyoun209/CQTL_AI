#!/bin/bash

ml plink/2.00a-20220129
ml samtools/1.21
ml bcftools/1.12  # if needed

# Create output directories for recoded files recursively
mkdir -p /work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/vcf/recode_012/pbs
mkdir -p /work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/vcf/recode_012/fnf

# Set the base directory for the input VCF files
base_dir="/work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/vcf"

# Process each condition separately
for CONDITION in fnf pbs; do
  echo "Processing condition: ${CONDITION}"

  # Find all chromosome-level VCF files for this condition.
  # Files should match the pattern: chr*.${CONDITION}.ASCounts.vcf.gz
  find "${base_dir}/${CONDITION}" -type f -name "chr*.${CONDITION}.ASCounts.vcf.gz" | sort -V | while read vcf; do
    base=$(basename "${vcf}")
    # Extract chromosome from filename (e.g. "chr4" from "chr4.pbs.ASCounts.vcf.gz")
    chr=$(echo "${base}" | cut -d'.' -f1)

    echo "  Processing file: ${base} (chromosome: ${chr})"

    # Run plink2 recode conversion with --const-fid 0 so that FID is set to 0 
    plink2 --vcf "${vcf}" \
           --const-fid 0 \
           --recode A-transpose \
           --keep-allele-order \
           --out /work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/vcf/recode_012/${CONDITION}/recodeA_${chr}_${CONDITION}

    # Remove the unwanted "0_" prefix from sample names in the header.
    # This regex replaces any occurrence of start-of-line or a tab followed by "0_" with just start-of-line or tab.
    sed -i '1s/\(^\|\t\)0_/\1/g' /work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/vcf/recode_012/${CONDITION}/recodeA_${chr}_${CONDITION}.traw
  done
done
