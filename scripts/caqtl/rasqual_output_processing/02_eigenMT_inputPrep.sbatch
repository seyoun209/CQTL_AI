#!/bin/bash
#SBATCH -p general
#SBATCH --mem=20g
#SBATCH -t 2-00:00:00
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -J eigenMT_inputprep
#SBATCH -o /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A.out
#SBATCH -e /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A.err

ml r/4.4.0

Rscript /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/02_eigenMT_inputPrep.R $1 $2 $3 $4
