#!/usr/bin/env Rscript

# Script to process multiple RASQUAL permutation files and adjust p-values for tied SNPs
library(data.table)

# Define parameters
base_dir <- "/work/users/s/e/seyoun/CQTL_AI"  # Adjust this to your base directory
input_dir <- file.path(base_dir, "output/rasqual_qtl/permutation_results")
output_dir <- file.path(base_dir, "output/rasqual_qtl/permutation_results/combined")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parse command line arguments for specific condition and PC (optional)
args <- commandArgs(trailingOnly = TRUE)
specific_condition <- NULL
specific_pc <- NULL

if (length(args) >= 1) {
  specific_condition <- args[1]
  cat("Will only process condition:", specific_condition, "\n")
}

if (length(args) >= 2) {
  specific_pc <- args[2]
  cat("Will only process PC:", specific_pc, "\n")
}

# Get all permutation files
permutation_files <- list.files(
  path = input_dir, 
  pattern = "_RASQUALResults_permutation_25kb.txt$", 
  recursive = TRUE, 
  full.names = TRUE
)

cat("Found", length(permutation_files), "permutation files\n")

# Extract condition and PC from filenames
file_info <- data.table(
  filepath = permutation_files,
  filename = basename(permutation_files)
)

# Parse filename to get chromosome, condition, and PC
file_info[, c("chr", "condition", "pc", "rest") := tstrsplit(filename, "_", keep = 1:4)]

# Filter for specific condition and PC if provided
if (!is.null(specific_condition)) {
  file_info <- file_info[condition == specific_condition]
}
if (!is.null(specific_pc)) {
  file_info <- file_info[pc == specific_pc]
}

# Get unique combinations of condition and PC
condition_pc_combos <- unique(file_info[, .(condition, pc)])

cat("Processing", nrow(condition_pc_combos), "condition-PC combinations\n")

# Process each condition-PC combination
for (i in 1:nrow(condition_pc_combos)) {
  current_condition <- condition_pc_combos$condition[i]
  current_pc <- condition_pc_combos$pc[i]
  
  cat(sprintf("Processing %s_%s (%d of %d)\n", 
              current_condition, current_pc, i, nrow(condition_pc_combos)))
  
  # Get all files for this condition and PC
  current_files <- file_info[condition == current_condition & pc == current_pc, filepath]
  
  # Initialize combined data table
  combined_data <- NULL
  
  # Process each file
  for (f in current_files) {
    cat("  Reading file:", basename(f), "\n")
    
    # Try to read the file
    tryCatch({
      # Determine if it's the first file (to get column names)
      is_first_file <- is.null(combined_data)
      
      # Read file
      current_data <- fread(f, header = FALSE, stringsAsFactors = FALSE)
      
      # Set column names if this is the first file
      if (is_first_file) {
        col_names <- c("Feature", "rs_ID", "Chromosome", "SNP_Position", "Ref_Allele", 
                       "Alt_allele", "Allele_Frequency", "HWE_Chi_Square", 
                       "Imputation_Quality", "Log10_BH_Qvalue", 
                       "Chi_square_statistic", "Effect_Size", 
                       "Error_Rate", "Ref_Bias", 
                       "Overdispersion", "SNP_ID_Region", "Feature_SNPs", 
                       "Tested_SNPs", "Null_Iterations", "Alt_Iterations", 
                       "Random_Ties", "Null_LogLikelihood", "Convergence", 
                       "FSNP_Correlation", "RSNP_Correlation")
        
        # Ensure we don't try to set more names than columns
        setnames(current_data, col_names[1:min(length(col_names), ncol(current_data))])
      } else {
        # Use same names as previously established
        setnames(current_data, names(combined_data)[1:min(length(names(combined_data)), ncol(current_data))])
      }
      
      # Calculate p-value from chi-square statistic
      current_data[, PValue := pchisq(as.numeric(Chi_square_statistic), df = 1, lower.tail = FALSE)]
      
      # Combine with previous data
      if (is_first_file) {
        combined_data <- current_data
      } else {
        combined_data <- rbindlist(list(combined_data, current_data), use.names = TRUE, fill = TRUE)
      }
      
      # Clean up to save memory
      rm(current_data)
      gc()
      
    }, error = function(e) {
      cat("  Error processing file:", basename(f), ":", e$message, "\n")
    })
  }
  
  # Skip if no combined data
  if (is.null(combined_data) || nrow(combined_data) == 0) {
    cat("  No data for", current_condition, current_pc, "\n")
    next
  }
  
  # Filter out invalid rows
  combined_data <- combined_data[rs_ID != "SKIPPED" & SNP_Position != "-1"]
  
  if (nrow(combined_data) == 0) {
    cat("  No valid data after filtering for", current_condition, current_pc, "\n")
    next
  }
  
  # Calculate distance from peak center for features formatted like "chr1_733775_734064"
  cat("  Calculating distances from peak centers...\n")
  
  # Parse the feature format - extract chromosome, start and end positions
  combined_data[, c("chr", "Peak_Start", "Peak_End") := tstrsplit(Feature, "_", keep = 1:3)]
  
  # Convert to numeric
  combined_data[, `:=`(
    Peak_Start = as.numeric(Peak_Start),
    Peak_End = as.numeric(Peak_End)
  )]
  
  # Calculate peak center and distance
  combined_data[, Peak_Center := round((Peak_End - Peak_Start)/2 + Peak_Start)]
  combined_data[, Distance_From_Peak := abs(as.numeric(SNP_Position) - Peak_Center)]
  
  # Clean up temporary columns
  combined_data[, c("chr", "Peak_Start", "Peak_End", "Peak_Center") := NULL]
  
  # Create combined file output
  combined_file <- file.path(output_dir, 
                             paste0(current_condition, "_", current_pc, "_permutation_combined.txt"))
  
  cat("  Writing combined file with", nrow(combined_data), "entries\n")
  fwrite(combined_data, combined_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # Adjust p-values for tied SNPs
  cat("  Adjusting p-values for tied SNPs...\n")
  
  # Get all unique features
  features <- unique(combined_data$Feature)
  cat("  Processing", length(features), "unique features for tie adjustment\n")
  
  # Create a copy of the data for adjustment
  adjusted_data <- copy(combined_data)
  
  # Track statistics
  features_with_ties <- 0
  total_ties_adjusted <- 0
  
  # Process each feature
  for (j in 1:length(features)) {
    if (j %% 1000 == 0 || j == length(features)) {
      cat("    Processed", j, "of", length(features), "features\n")
    }
    
    current_feature <- features[j]
    feature_data <- adjusted_data[Feature == current_feature]
    
    # Skip if only one SNP
    if (nrow(feature_data) <= 1) next
    
    # First sort by p-value and distance from peak center, exactly like your original code
    # This matches your: RASQUAL_filtered <- RASQUAL_filtered[with(RASQUAL_filtered, order(PValue, Distance_from_PeakCenter)),]
    setorder(feature_data, PValue, Distance_From_Peak)
    
    # Then simply add 1E-15 to all p-values except the first one
    # This matches your: RASQUAL_filtered$PValue[2:nrow(RASQUAL_filtered)] <- RASQUAL_filtered$PValue[2:nrow(RASQUAL_filtered)] + 1E-15
    if (nrow(feature_data) > 1) {
      # Add 1E-15 to all except the first row
      feature_data$PValue[2:nrow(feature_data)] <- feature_data$PValue[2:nrow(feature_data)] + 1E-15
      
      # Update statistics
      features_with_ties <- features_with_ties + 1
      total_ties_adjusted <- total_ties_adjusted + (nrow(feature_data) - 1)
      
      # Properly update in the main dataset using data.table syntax
      adjusted_data[Feature == current_feature, PValue := feature_data$PValue]
    }
  }
  
  cat("  Features with tied p-values:", features_with_ties, "out of", length(features),
      sprintf("(%.2f%%)", 100 * features_with_ties / length(features)), "\n")
  cat("  Total tied p-values adjusted:", total_ties_adjusted, "\n")
  
  # Write adjusted RASQUAL results
  adjusted_file <- file.path(output_dir, 
                             paste0(current_condition, "_", current_pc, "_permutation_adjusted.txt"))
  
  cat("  Writing adjusted RASQUAL results\n")
  fwrite(adjusted_data, adjusted_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # Clean up
  rm(combined_data, adjusted_data)
  gc()
}

cat("Processing completed!\n")