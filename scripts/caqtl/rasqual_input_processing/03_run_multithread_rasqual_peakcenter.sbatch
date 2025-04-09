#!/bin/bash
# Load required modules
module add tabix/0.2.6
module add gsl/2.7

# Get arguments from command line
peak_file="${1}"
peak_counts="${2}"
sample_offsets="${3}"
covariates="${4}"
VCF_dir="${5}"
VCF="${6}"
samplenumber="${7}"
Output="${8}"
pc_value="${9}" # Parameter for PC value
permutation_flag="${10}"  # Parameter for permutation flag

# Calculate total_covariates based on pc_value
total_covariates=$((pc_value + 5))
echo "total_covariates: ${total_covariates}"


# Print all parameters for debugging
echo "=== RASQUAL Script Parameters ==="
echo "peak_file: ${peak_file}"
echo "peak_counts: ${peak_counts}"
echo "sample_offsets: ${sample_offsets}"
echo "covariates: ${covariates}"
echo "VCF_dir: ${VCF_dir}"
echo "VCF: ${VCF}"
echo "samplenumber: ${samplenumber}"
echo "Output: ${Output}"
echo "pc_value: ${pc_value}"
echo "permutation_flag: ${permutation_flag}"
echo "================================="

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

# Calculate peak center
peak_center=$(awk -v start=${peak_start} -v end=${peak_end} 'BEGIN {print int((start + end) / 2)}')

# Set window to 1kb or 100kb around peak center
cisWindow_start=$(awk -v center=${peak_center} 'BEGIN {print center - 100000}')
cisWindow_end=$(awk -v center=${peak_center} 'BEGIN {print center + 100000}')

num_cisSNPs=$(tabix ${VCF_dir}/${VCF} ${peak_chr}:${cisWindow_start}-${cisWindow_end} | wc -l)
num_peakSNPs=$(tabix ${VCF_dir}/${VCF} ${peak_chr}:${peak_start}-${peak_end} | wc -l)

WindowSize=$(($cisWindow_end - $cisWindow_start))

# Print peak info for debugging (for all peaks)
echo "Processing peak #${peak_rowNum}: ${peak_ID} (${peak_chr}:${peak_start}-${peak_end})"
echo "  Peak center: ${peak_center}"
echo "  1kb Window: ${cisWindow_start}-${cisWindow_end} (size: ${WindowSize})"
echo "  SNPs in window: ${num_cisSNPs}, SNPs in peak: ${num_peakSNPs}"
if [ "$permutation_flag" == "--random-permutation" ]; then
echo "  Using random permutation for this peak"
fi

# Added permutation flag to RASQUAL command if specified
if [ "$permutation_flag" == "--random-permutation" ]; then
tabix ${VCF_dir}/${VCF} ${peak_chr}:${cisWindow_start}-${cisWindow_end} \
| /users/s/e/seyoun/tools/rasqual/bin/rasqual -s ${peak_start} -e ${peak_end} -j ${peak_rowNum} \
-y ${peak_counts} -k ${sample_offsets} \
-n ${samplenumber} -l ${num_cisSNPs} -m ${num_peakSNPs} \
-f ${peak_ID} -x ${covariates} -p ${total_covariates} -w ${WindowSize} \
-q 0.0 -a 0.00 -h 0.00 --min-coverage-depth 0.0 -z --n-threads 10 --random-permutation \
>> ${Output}
else
  tabix ${VCF_dir}/${VCF} ${peak_chr}:${cisWindow_start}-${cisWindow_end} \
| /users/s/e/seyoun/tools/rasqual/bin/rasqual -s ${peak_start} -e ${peak_end} -j ${peak_rowNum} \
-y ${peak_counts} -k ${sample_offsets} \
-n ${samplenumber} -l ${num_cisSNPs} -m ${num_peakSNPs} \
-f ${peak_ID} -x ${covariates} -p ${total_covariates} -w ${WindowSize} \
-q 0.0 -a 0.00 -h 0.00 --min-coverage-depth 0.0 -z --n-threads 10 \
>> ${Output}
fi

peak_rowNum=$(( peak_rowNum+1 ))
done < ${peak_file}

echo "RASQUAL processing completed. Total peaks processed: $((peak_rowNum-1))"

