#!/bin/bash
#SBATCH -p general
#SBATCH --mem=5g
#SBATCH -t 1-00:00:00
#SBATCH -J rasqual_vcf
#SBATCH -o /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A_%a.out
#SBATCH -e /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A_%a.err

# filepath: /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_input_processing/01_rasqual_atac_vcf_prep.sbatch

# Input parameters
input_vcf=$1
chr=$2
bam_list=$3
output_vcf=$4
temp_dir=$5

# Print info for debugging
echo "======== RASQUAL VCF Preparation ========"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Input VCF: $input_vcf"
echo "Chromosome: $chr"
echo "BAM list: $bam_list"
echo "Output VCF: $output_vcf"
echo "Temp directory: $temp_dir"
echo "========================================"

# Load required modules and set paths
export RASQUALDIR=/users/s/e/seyoun/tools/rasqual/
module load samtools/1.21
module load bcftools/1.21
module load htslib/1.21

# Echo loaded modules for debugging
echo "Loaded modules:"
module list

# Echo BAM file count
echo "Number of BAM files: $(wc -l < $bam_list)"
echo "First few BAM files:"
head -3 $bam_list

# Make sure temp directory exists
mkdir -p ${temp_dir}

# Create full path for temp file
temp_vcf="${temp_dir}/temp_${chr}.vcf.gz"

# Process VCF - extract chromosome
echo "Extracting chromosome $chr from VCF..."
bcftools view ${input_vcf} ${chr} -Oz -o ${temp_vcf}

# Check if extraction was successful
if [ $? -ne 0 ] || [ ! -f "${temp_vcf}" ]; then
    echo "ERROR: Failed to extract chromosome from VCF"
    echo "Input VCF exists: $(ls -la ${input_vcf})"
    echo "Temp VCF created: $(ls -la ${temp_vcf} 2>/dev/null || echo "FILE NOT CREATED")"
    exit 1
fi

# Create tabix index
echo "Creating tabix index for temp VCF..."
tabix -p vcf ${temp_vcf}

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create tabix index for temp VCF"
    exit 1
fi

# Print debug info
echo "Debug information:"
echo "Current directory: $(pwd)"
echo "Temp VCF size: $(ls -lh ${temp_vcf} | awk '{print $5}')"
echo "VCF sites: $(bcftools view -H ${temp_vcf} | wc -l)"

# Run RASQUAL AS counts
echo "Running RASQUAL ASCounts..."
bash $RASQUALDIR/src/ASVCF/createASVCF.sh paired_end ${bam_list} ${temp_vcf} ${output_vcf} atac

# Check if ASCounts was successful
if [ $? -ne 0 ] || [ ! -f "${output_vcf}" ]; then
    echo "ERROR: Failed to create ASCounts VCF"
    exit 1
fi

# Create tabix index for output VCF
echo "Creating tabix index for output VCF..."
tabix -p vcf ${output_vcf}

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create tabix index for output VCF"
    exit 1
fi

echo "Output VCF size: $(ls -lh ${output_vcf} | awk '{print $5}')"
echo "Output VCF sites: $(bcftools view -H ${output_vcf} | wc -l)"
echo "Completed successfully at $(date)"

# Optional: Clean up temp files if desired
# rm -f ${temp_vcf} ${temp_vcf}.tbi
