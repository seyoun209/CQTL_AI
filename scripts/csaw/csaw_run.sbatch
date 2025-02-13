#!/bin/bash
#SBATCH -J csaw
#SBATCH -p general
#SBATCH -n 3
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -t 3-00:00:00
#SBATCH -e /work/users/s/e/seyoun/CQTL_AI/output/logs/csaw_ankle_e.%j
#SBATCH -o /work/users/s/e/seyoun/CQTL_AI/output/logs/csaw_ankle_o.%j

module load r/4.4.0

Rscript ankle_only_csaw.R


