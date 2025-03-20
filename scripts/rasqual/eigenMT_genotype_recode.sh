#!/bin/bash
  
# Load required modules
ml plink/2.00a-20220129
ml samtools/1.21

# Set directories
recode_dir="/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf/recode_012"
base_dir="/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf"

# Create output directories
mkdir -p ${recode_dir}/pbs
mkdir -p ${recode_dir}/fnf

# Define chromosomes in the correct order
CHROMOSOMES=($(echo {1..22} X))

# Process each condition
for CONDITION in fnf pbs; do
    echo "Processing condition: ${CONDITION}"
    
    # Process each chromosome separately
    for CHR in "${CHROMOSOMES[@]}"; do
        echo "Processing chromosome ${CHR} for ${CONDITION}"
        
        # Define VCF file for this chromosome
        VCF_FILE="${base_dir}/${CONDITION}/chr${CHR}.${CONDITION}.ASCounts.vcf.gz"
        
        # Check if VCF file exists
        if [ ! -f "${VCF_FILE}" ]; then
            echo "Warning: VCF file not found: ${VCF_FILE}"
            continue
        fi
        
        # Convert to PLINK additive format for this chromosome
        echo "Converting chr${CHR} to PLINK additive format..."
        plink2 --vcf ${VCF_FILE} \
               --double-id \
               --recode A-transpose \
               --keep-allele-order \
               --out ${recode_dir}/${CONDITION}/recodeA_chr${CHR}_${CONDITION}
        
        echo "Completed processing chromosome ${CHR} for ${CONDITION}"
    done
    
    echo "Completed processing all chromosomes for ${CONDITION}"
done

echo "All processing complete!"
