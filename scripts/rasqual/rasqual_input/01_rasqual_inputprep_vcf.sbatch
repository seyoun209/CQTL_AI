#!/bin/bash
#SBATCH -p general
#SBATCH --mem=5g
#SBATCH -t 1-00:00:00

# Input parameters
input_vcf=$1
chr=$2
bam_list=$3
output_vcf=$4
temp_dir=$5

# Load required modules and set paths
export RASQUALDIR=/users/s/e/seyoun/tools/rasqual/
module load samtools/1.21

# Create temp directory if it doesn't exist
mkdir -p ${temp_dir}

# Create full path for temp file
temp_vcf="${temp_dir}/temp_chr${chr}.vcf.gz"

# Process VCF (using full path)
bcftools view ${input_vcf} ${chr} -Oz -o ${temp_vcf}
tabix -p vcf ${temp_vcf}

# Print debug info
echo "Debug information:"
echo "Current directory: $(pwd)"
echo "Temp VCF path: ${temp_vcf}"
echo "Checking if temp VCF exists: $(ls -l ${temp_vcf})"

# Run RASQUAL AS counts with full path to temp file
bash $RASQUALDIR/src/ASVCF/createASVCF.sh paired_end ${bam_list} ${temp_vcf} ${output_vcf} atac


tabix -p vcf ${output_vcf}
