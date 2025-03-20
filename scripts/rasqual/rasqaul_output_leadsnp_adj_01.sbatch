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

# Process both fnf and pbs conditions
for tissue in fnf pbs; do
  for pc in {1..10}; do
      condition="${tissue}"
      pc="PC${pc}"
      sbatch -p general -N 1 --mem=40g -n 1 -t 2- \
        --wrap="Rscript /work/users/s/e/seyoun/CQTL_AI/scripts/rasqual/rasqaul_output_leadsnp_adj_01.R  $condition $pcs"
done
done
