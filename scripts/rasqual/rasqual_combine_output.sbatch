#!/bin/bash
#SBATCH -J merge_rasqual
#SBATCH -p general
#SBATCH -n 3
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -t 3-00:00:00
#SBATCH --output=./logs/combine_rasqual_%j.out
#SBATCH --error=./logs/combine_rasqual_%j.err

# Load R module if needed (adjust based on your system)
module load r/4.4.0

# Run the R script
#Rscript rasqualtoEigenMT_input.R

Rscript eigenMT_leadsnpadjust.R
