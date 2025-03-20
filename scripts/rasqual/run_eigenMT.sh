#!/bin/bash

ml python/3.9.6

dir="/users/s/e/seyoun/tools/eigenMT"

QTL_file="${1}";
Gen_file="${2}";
Genpos_file="${3}";
Phepos_file="${4}";
Out_file="${5}";
Chr="${6}";

python ${dir}/eigenMT.py --CHROM ${Chr} \
        --QTL ${QTL_file} \
        --GEN ${Gen_file} \
        --GENPOS ${Genpos_file} \
        --PHEPOS ${Phepos_file} \
        --var_thresh 0.99 \
        --cis_dist 50000 \
        --OUT ${Out_file} \
        --window 200


#dir_input="/work/users/s/e/seyoun/CQTL_AI/output/eigenMT/pbs"
#python ${dir}/eigenmt_debug.py ${dir_input}/qtl_chr22_pbs_PC1.txt ${dir_input}/gen.positions_chr22_pbs_PC1.txt ${dir_input}/phe.positions_chr22_pbs_PC1.txt ${dir_input}/genotypes_chr22_pbs.txt
