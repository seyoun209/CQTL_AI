#!/bin/bash

ml plink/1.90b3
ml samtools/1.21
ml vcftools/0.1.15

input_vcf=$1
sample_list=$2
output_temp=$3
output_final=$4

mapping_dir="$(dirname "$output_final")/mappings"
mkdir -p "$mapping_dir"
mapping_file="${mapping_dir}/sample_name_mapping.txt"

# Get all VCF samples
vcf_samples=$(bcftools query -l "$input_vcf")
echo "Found VCF samples:"
echo "$vcf_samples"

# First pass: Match exact 4-digit AM numbers
while IFS=' ' read -r sample_name category; do
    if [[ $sample_name =~ AM([0-9]{4}) ]]; then
        am_num=${BASH_REMATCH[1]}
        # Adjust pattern: use any leading number instead of literal "0_"
        vcf_name=$(echo "$vcf_samples" | grep -P "^[0-9]+_AM${am_num}$")

        if [ ! -z "$vcf_name" ]; then
            echo "Exact match found: $vcf_name -> $sample_name"
            echo -e "$vcf_name\t$sample_name" >> "$mapping_file"
            # Remove matched sample
            vcf_samples=$(echo "$vcf_samples" | grep -v "^$vcf_name$")
        fi
    fi
done < "$sample_list"

# Second pass: Match Phanstiel format only for unmatched samples
while IFS=' ' read -r sample_name category; do
    # Only process if not already matched
    if ! grep -q "$sample_name" "$mapping_file"; then
        if [[ $sample_name =~ AM([0-9]{4}) ]]; then
            am_num=${BASH_REMATCH[1]}
            last_two=${am_num: -2}
            # Adjust pattern similarly here
            vcf_name=$(echo "$vcf_samples" | grep -P "^[0-9]+_Phanstiel_${last_two}$")

            if [ ! -z "$vcf_name" ]; then
                echo "Phanstiel match found: $vcf_name -> $sample_name"
                echo -e "$vcf_name\t$sample_name" >> "$mapping_file"
                vcf_samples=$(echo "$vcf_samples" | grep -v "^$vcf_name$")
            fi
        fi
    fi
done < "$sample_list"

# Display final mapping
echo "Final mappings:"
cat "$mapping_file"

if [ ! -s "$mapping_file" ]; then
    echo "ERROR: Mapping file is empty. Check sample names and VCF content." >&2
    exit 1
fi

bcftools reheader --samples "$mapping_file" "$input_vcf" -o "$output_temp"
tabix -p vcf "$output_temp"
bcftools annotate "$output_temp" --rename-chrs /work/users/s/e/seyoun/CQTL_AI/scripts/chrnm.txt -Oz -o "$output_final"
tabix -p vcf "$output_final"
