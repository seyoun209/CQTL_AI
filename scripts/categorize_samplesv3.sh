#!/bin/bash

# Define the directory containing the sample files
sample_directory="/work/users/s/e/seyoun/CQTL_AI/output/fastq"

# Define categories
declare -A category_mapping=(
    ["CTL"]="CTL"
    ["FNF"]="FNF"
)

# Replacement mapping for specific samples
declare -A replacement_mapping=(
    ["CQTL_AM7754_CTL_Ankle_1"]="CQTL_AM7754_CTL_Ankle_replicate"
    ["CQTL_AM7755_CTL_Ankle_1"]="CQTL_AM7755_CTL_Ankle_replicate"
    ["CQTL_AM7778_CTL_Ankle_1"]="CQTL_AM7778FE_CTL_Femur_1"
    ["CQTL_AM7754_FNF_Ankle_1"]="CQTL_AM7754_FNF_Ankle_replicate"
    ["CQTL_AM7755_FNF_Ankle_1"]="CQTL_AM7755_FNF_Ankle_replicate"
    ["CQTL_AM7778_FNF_Ankle_1"]="CQTL_AM7778FE_FNF_Femur_1"
)

# Get all fastq files
sample_files=("$sample_directory"/*_R1.fastq.gz)

# Initialize arrays for categorized samples
declare -A categorized_samples
synovium_samples=()

# Process each sample file
for file in "${sample_files[@]}"; do
    # Extract base sample name
    sample_name=$(basename "$file")
    sample_name=${sample_name%%_R1.fastq.gz}

    # Handle synovium samples
    if [[ $sample_name == *"synovium"* ]]; then
        synovium_samples+=("$sample_name")
        continue
    fi

    # Check if the sample is in the replacement mapping
    if [[ -n ${replacement_mapping[$sample_name]} ]]; then
        sample_name=${replacement_mapping[$sample_name]}
    fi

    # Categorize samples into CTL or FNF
    for category in "${!category_mapping[@]}"; do
        if [[ $sample_name == *"${category_mapping[$category]}"* ]]; then
            categorized_samples[$category]+="$sample_name $category"$'\n'
            break
        fi
    done
done

# Save categorized samples into separate files, removing duplicates
for category in "${!categorized_samples[@]}"; do
    output_file="/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/${category}_samples_with_replacements.txt"
    printf "%s" "${categorized_samples[$category]}" | sort | uniq > "$output_file"
done

# Save synovium samples into a separate file
synovium_file="/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/Synovium_samples.txt"
printf "%s\n" "${synovium_samples[@]}" > "$synovium_file"

echo "Sample processing complete. Files saved in /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/"

