#### filepath: /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/00_rasqual_combine_output.R
#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(stats)
library(data.table)

# Optional: Increase printed precision (for debugging ties)
options(digits=20)

### Define base directories and window type
base_dir <- "/work/users/s/e/seyoun/CQTL_AI"
window_type <- "100kb"   # Change to "100kb" when processing 100kb files, otherwise "1kb"
input_dir <- file.path(base_dir, paste0("output/rasqual_output/window_", window_type))
output_dir <- file.path(base_dir, paste0("output/rasqual_output/window_", window_type, "/combined_", window_type))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

### Parse command line arguments for specific condition and PC (optional)
args <- commandArgs(trailingOnly = TRUE)
specific_condition <- if(length(args) >= 1) args[1] else NULL
if (!is.null(specific_condition)) cat("Will only process condition:", specific_condition, "\n")
specific_pc <- if(length(args) >= 2) args[2] else NULL
if (!is.null(specific_pc)) cat("Will only process PC:", specific_pc, "\n")

### Get all RASQUAL result files
pattern_str <- paste0("_", window_type, "\\.txt$")
result_files <- list.files(path = input_dir, pattern = pattern_str, recursive = TRUE, full.names = TRUE)
cat("Found", length(result_files), window_type, "RASQUAL result files\n")

### Create a data.table of file information and parse filenames
file_info <- data.table(
  filepath = result_files,
  filename = basename(result_files)
)
# Expected filename format: <condition>_chr<chr>_<pc>_<window>.txt e.g., "pfnf_chr10_pc0_1kb.txt"
file_info[, c("condition", "chr", "pc", "window") := tstrsplit(filename, "_", keep = 1:4)]
file_info[, chr := str_remove(chr, "^chr")]

### Filter for provided condition/PC, if given
if (!is.null(specific_condition)) file_info <- file_info[condition == specific_condition]
if (!is.null(specific_pc)) file_info <- file_info[pc == specific_pc]

### Get unique combinations of condition and PC
condition_pc_combos <- unique(file_info[, .(condition, pc)])
cat("Processing", nrow(condition_pc_combos), "condition-PC combinations\n")

### Loop over each condition and PC combination
for (i in 1:nrow(condition_pc_combos)) {
  current_condition <- condition_pc_combos$condition[i]
  current_pc <- condition_pc_combos$pc[i]
  cat(sprintf("Processing %s_%s (%d of %d)\n",
              current_condition, current_pc, i, nrow(condition_pc_combos)))
  
  ### Get all files for this condition and PC
  current_files <- file_info[condition == current_condition & pc == current_pc, filepath]
  
  combined_data <- NULL
  
  ### Loop over files (chromosomes) and merge
  for (f in current_files) {
    cat("  Reading file:", basename(f), "\n")
    tryCatch({
      is_first_file <- is.null(combined_data)
      # Read file without header (adjust if headers present)
      current_data <- fread(f, header = FALSE, stringsAsFactors = FALSE)
      
      ### Define column names on the first file (adjust order as needed)
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
        setnames(current_data, col_names[1:min(length(col_names), ncol(current_data))])
      } else {
        setnames(current_data, names(combined_data)[1:min(length(names(combined_data)), ncol(current_data))])
      }
      
      ### Filter out rows with invalid SNP_Position and apply quality filters
      current_data <- current_data[!is.na(SNP_Position) & SNP_Position != "-1"]
      current_data <- current_data[RSNP_Correlation >= 0.8]
      current_data <- current_data[Imputation_Quality >= 0.7]
      
      ### Compute p-value from Chi_square_statistic (df = 1)
      current_data[, PValue := pchisq(as.numeric(Chi_square_statistic), df = 1, lower.tail = FALSE)]
      current_data[, PValue := ifelse(is.na(PValue) | PValue < 0, 1, PValue)] # Ensure no negative p-values
      ### Merge into combined_data
      if (is_first_file) {
        combined_data <- current_data
      } else {
        combined_data <- rbindlist(list(combined_data, current_data), use.names = TRUE, fill = TRUE)
      }
      rm(current_data)
      gc()
      
    }, error = function(e) {
      cat("  Error processing file:", basename(f), ":", e$message, "\n")
    })
  }
  
  ### If no data, skip to next combination
  if (is.null(combined_data) || nrow(combined_data) == 0) {
    cat("  No data for", current_condition, current_pc, "\n")
    next
  }
  
  ### Calculate distance from peak center
  cat("  Calculating distances from peak centers...\n")
  # Expected feature format: e.g. "chr1_733775_734064"
  combined_data[, c("chr_temp", "Peak_Start", "Peak_End") := tstrsplit(Feature, "_", keep = 1:3)]
  combined_data[, `:=`(
    Peak_Start = as.numeric(Peak_Start),
    Peak_End = as.numeric(Peak_End)
  )]
  combined_data[, Peak_Center := round((Peak_End - Peak_Start)/2 + Peak_Start)]
  combined_data[, Distance_From_Peak := abs(as.numeric(SNP_Position) - Peak_Center)]
  combined_data[, c("chr_temp", "Peak_Start", "Peak_End", "Peak_Center") := NULL]
  
  ### Write combined file (with window tag)
  combined_file <- file.path(output_dir,
                             paste0(current_condition, "_", current_pc, "_", window_type, "_combined.txt"))
  cat("  Writing combined file with", nrow(combined_data), "entries\n")
  fwrite(combined_data, combined_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  ### Adjust tied p-values for each unique feature
  cat("  Adjusting p-values for tied SNPs...\n")
  features <- unique(combined_data$Feature)
  cat("  Processing", length(features), "unique features for tie adjustment\n")
  
  adjusted_data <- copy(combined_data)
  features_with_ties <- 0
  total_ties_adjusted <- 0
  
  for (j in 1:length(features)) {
    if (j %% 1000 == 0 || j == length(features)) {
      cat("    Processed", j, "of", length(features), "features\n")
    }
    
    current_feature <- features[j]
    # Subset data for current feature
    feature_data <- adjusted_data[Feature == current_feature]
    
    if (nrow(feature_data) <= 1) next
    
    # Sort feature_data by PValue and Distance_From_Peak
    setorder(feature_data, PValue, Distance_From_Peak)
    
    # Identify duplicated (tied) p-values (beyond the first occurrence)
    dup_indices <- which(duplicated(feature_data$PValue))
    
    if (length(dup_indices) > 0) {
      feature_data$PValue[dup_indices] <- feature_data$PValue[dup_indices] + 1e-15
      features_with_ties <- features_with_ties + length(dup_indices)
      total_ties_adjusted <- total_ties_adjusted + length(dup_indices)
      
      ### Update adjusted_data for this feature
      adjusted_data[Feature == current_feature, PValue := feature_data$PValue]
    }
  }
  
  cat("  Features with tied p-values adjusted:", features_with_ties, "out of", length(features),
      sprintf("(%.2f%%)", 100 * features_with_ties / length(features)), "\n")
  cat("  Total tied p-values adjusted:", total_ties_adjusted, "\n")
  
  ### Write adjusted results (with window tag)
  adjusted_file <- file.path(output_dir,
                             paste0(current_condition, "_", current_pc, "_", window_type, "_adjusted.txt"))
  cat("  Writing adjusted RASQUAL results to file\n")
  fwrite(adjusted_data, adjusted_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  cat("  Successfully wrote adjusted file with", nrow(adjusted_data), "rows\n")
  rm(combined_data, adjusted_data)
  gc()
}
cat("Processing completed!\n")
