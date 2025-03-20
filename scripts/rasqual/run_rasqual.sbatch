#!/bin/bash

#. /etc/profile.d/modulex.sh
module add tabix/0.2.6
module add gsl/2.7

# I will do 1kb and 100kb both

peak_file="${1}"
peak_counts="${2}"
sample_offsets="${3}"
covariates="${4}"
VCF_dir="${5}"
VCF="${6}"
Output="${8}"

peak_rowNum=1

while read peak; do
    # Skip the header line
    is_header=$(echo -e "${peak}" | grep -cw "PeakID")
    if [ ${is_header} -eq 1 ]
    then
        continue
    fi

    params=(${peak//\t/ })
    peak_ID=${params[0]}
    peak_chr=${params[1]}
    peak_start=${params[2]}
    peak_end=${params[3]}
    cisWindow_start=$(awk -v m=${peak_start} 'BEGIN {print m-25000}')
    cisWindow_end=$(awk -v m=${peak_end} 'BEGIN {print m+25000}')

    num_cisSNPs=$(tabix ${VCF_dir}/${VCF} ${peak_chr}:${cisWindow_start}-${cisWindow_end} | wc -l)
    num_peakSNPs=$(tabix ${VCF_dir}/${VCF} ${peak_chr}:${peak_start}-${peak_end} | wc -l)

    WindowSize=$(($cisWindow_end - $cisWindow_start))

    tabix ${VCF_dir}/${VCF} ${peak_chr}:${cisWindow_start}-${cisWindow_end} \
    | /users/s/e/seyoun/tools/rasqual/bin/rasqual -s ${peak_start} -e ${peak_end} -j ${peak_rowNum} \
    -y ${peak_counts} -k ${sample_offsets} \
    -n ${7} -l ${num_cisSNPs} -m ${num_peakSNPs} \
    -f ${peak_ID} -x ${covariates} -p 10 -w ${WindowSize} \
    -q 0.0 -a 0.00 -h 0.00 --min-coverage-depth 0.0 -z --n-threads 10 \
    >> ${Output}

    peak_rowNum=$(( peak_rowNum+1 ))

done < ${peak_file}
