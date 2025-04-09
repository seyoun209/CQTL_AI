#!/bin/bash
#SBATCH -p general
#SBATCH --mem=50g
#SBATCH -t 5-00:00:00
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -J response_qtl
#SBATCH -o /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A.out
#SBATCH -e /work/users/s/e/seyoun/CQTL_AI/output/logs/%x_%A.err

ml r/4.4.0

Rscript /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/response_caQTL/01_response_qtl_run.R

