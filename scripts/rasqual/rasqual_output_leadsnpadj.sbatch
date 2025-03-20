#!/bin/bash
#SBATCH -J leadsnp_rasqual
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=25g
#SBATCH -t 3-00:00:00
#SBATCH --output=./logs/leadsnp_rasqual_%j.out
#SBATCH --error=./logs/leadsnp_rasqual_%j.err

module load r/4.4.0

# Create logs directory
mkdir -p ./logs

# Process both fnf and pbs conditions
for tissue in fnf pbs; do
  for pcnum in {1..10}; do
      pc="PC${pcnum}"
      job_name="leadsnp_${tissue}_${pc}"
      
      echo "Submitting job for ${tissue} ${pc}"
      
      sbatch -p general -N 1 --mem=20g -n 1 -t 1-00:00:00 \
        --job-name=${job_name} \
        --output=./logs/${job_name}_%j.out \
        --error=./logs/${job_name}_%j.err \
        --wrap="Rscript /work/users/s/e/seyoun/CQTL_AI/scripts/rasqual/rasqual_combine_and_adjust.R ${tissue} ${pc}"
  done
done
