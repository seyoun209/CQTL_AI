#!/bin/bash
ml python/3.9.6

# Set the directory where eigenMT is located
eigenMT_dir="/users/s/e/seyoun/tools/eigenMT"

# Parse input arguments
QTL_file="${1}"
Gen_file="${2}"
Genpos_file="${3}"
Phepos_file="${4}"
Out_file="${5}"
Chr="${6}"

python ${eigenMT_dir}/eigenMT.py --CHROM ${Chr} \
        --QTL ${QTL_file} \
        --GEN ${Gen_file} \
        --GENPOS ${Genpos_file} \
        --PHEPOS ${Phepos_file} \
        --var_thresh 0.99 \
        --cis_dist 1000000 \
        --OUT ${Out_file} \
        --window 200