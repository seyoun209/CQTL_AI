#!/bin/bash

# Define the directory containing the sample files
sample_directory="/work/users/s/e/seyoun/CQTL_AI/output/fastq"

# Define categories
declare -A category_mapping=(
    ["CTL"]="CTL"
    ["FNF"]="FNF"
)

# Get all fastq files
sample_files=("$sample_directory"/*_R1.fastq.gz)

# Initialize arrays for categorized samples
declare -A categorized_samples
declare -A categorized_samples_v2
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

    # Handle Femur and replicate samples (v2)
    if [[ $sample_name == *"Femur"* ]] || [[ $sample_name == *"replicate"* ]]; then
        for category in "${!category_mapping[@]}"; do
            if [[ $sample_name == *"${category_mapping[$category]}"* ]]; then
                categorized_samples_v2[$category]+="$sample_name"$'\n'
                break
            fi
        done
        continue
    fi

    # Categorize remaining samples (CTL and FNF)
    for category in "${!category_mapping[@]}"; do
        if [[ $sample_name == *"${category_mapping[$category]}"* ]]; then
            categorized_samples[$category]+="$sample_name"$'\n'
            break
        fi
    done
done

# Save original categorized samples into separate files
for category in "${!categorized_samples[@]}"; do
    output_file="/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/${category}_samples_v2.txt"
    printf "%s" "${categorized_samples[$category]}" > "$output_file"
done

# Save Femur and replicate samples (v2) into separate files
for category in "${!categorized_samples_v2[@]}"; do
    output_file="/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/${category}_femur_replicate_samples_v2.txt"
    printf "%s" "${categorized_samples_v2[$category]}" > "$output_file"
done

# Save synovium samples into a separate file
synovium_file="/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/Synovium_samples.txt"
printf "%s\n" "${synovium_samples[@]}" > "$synovium_file"

