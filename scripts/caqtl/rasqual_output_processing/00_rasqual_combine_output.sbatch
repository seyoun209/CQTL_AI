#!/bin/bash
#SBATCH -p general
#SBATCH --mem=100g
#SBATCH -t 5-00:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -J rasqual_combine
#SBATCH -o /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A.out
#SBATCH -e /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A.err

ml r/4.4.0

Rscript /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/00_rasqual_combine_output.R $1 $2



#Rscript /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/00_rasqual_combine_output.R pbs pc0
