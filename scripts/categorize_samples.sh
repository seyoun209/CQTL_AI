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

# Process each sample file
for file in "${sample_files[@]}"; do
    # Extract base sample name
    sample_name=$(basename "$file")
    sample_name=${sample_name%%_R1.fastq.gz}
    
    # Skip if contains Femur, synovium, or replicate
    if [[ $sample_name == *"Femur"* ]] || [[ $sample_name == *"synovium"* ]] || \
       [[ $sample_name == *"replicate"* ]]; then
        continue
    fi
    
    # Find category (CTL or FNF)
    for category in "${!category_mapping[@]}"; do
        if [[ $sample_name == *"${category_mapping[$category]}"* ]]; then
            categorized_samples[$category]+="$sample_name $category"$'\n'
            break
        fi
    done
done

# Save categorized samples into separate files
for category in "${!categorized_samples[@]}"; do
    output_file="/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/${category}_samples.txt"
    printf "%s" "${categorized_samples[$category]}" > "$output_file"
done
