#!/usr/bin/env Rscript
#---------------------------------------------------------------
# RASQUAL Multi-threaded Job Submission Script for PBS/FNF
# Runs all PC values (0-20) for both conditions
#---------------------------------------------------------------
library(dplyr)
library(stringr)

# Clear environment
rm(list=ls())
options(stringsAsFactors=FALSE)

#---------------------------------------------------------------
# Configuration parameters
#---------------------------------------------------------------
# Set your conditions
conditions <- c("pbs", "fnf")

# Set chromosomes
chromosomes <- c(1:22)
#chromosomes <- c(22)
# Set PC values to test
pc_values <- 0:20  # Run all PCs from 0 to 20
#pc_values <- c(1)
# Set window size (in kb)
# Options: 1 = 1kb window, 100 = 100kb window
window_kb <- 100  # Change to 1 for 1kb window

# Base directories
base_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input"
output_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_output"
script_dir <- "/work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_input_processing"

# Resources for SLURM
slurm_params <- list(
  partition = "general",
  n_cores = 10,
  mem_gb = 5,
  time = "3-00:00:00"  # 3 days
)

# RASQUAL parameters
rasqual_params <- list(
  p = 10,  # Initial value, will be updated in loop
  a = 0.00,  # Alpha parameter
  h = 0.00,  # Lambda parameter 
  q = 0.0,  # Min. MAF
  min_cov_depth = 0.0  # Min. coverage depth
)

# Flag for permutations (set to TRUE to enable permutations)
run_permutations <- FALSE

#---------------------------------------------------------------
# Job submission function
#---------------------------------------------------------------
submit_rasqual_jobs <- function(
    condition,
    chromosome,
    window_kb,
    pc_value,
    run_permutation = FALSE
) {
  # Convert chromosome name to format used in files
  chr_name <- paste0("chr", chromosome)
  
  # File paths - UPDATED to use input_count folder
  peak_info <- file.path(base_dir, "input_count", paste0("peak_info_", chr_name, ".txt"))
  peak_counts <- file.path(base_dir, "input_count", paste0(chr_name, "_", condition, ".expression.bin"))
  sample_offsets <- file.path(base_dir, "input_count", paste0(chr_name, "_", condition, ".size_factors.bin"))
  covariates <- file.path(base_dir, "covar", paste0(condition, "_covariates_pc", pc_value, ".bin"))
  
  # VCF paths
  vcf_dir <- file.path(base_dir, "vcf", condition)
  vcf_file <- paste0(chr_name, ".", condition, ".ASCounts.vcf.gz")
  
  # Check if all files exist
  required_files <- c(peak_info, peak_counts, sample_offsets, covariates, file.path(vcf_dir, vcf_file))
  missing_files <- required_files[!file.exists(required_files)]
  
  if(length(missing_files) > 0) {
    warning("Missing files for ", condition, " ", chr_name, " PC", pc_value, ":")
    for(file in missing_files) {
      warning("  ", file)
    }
    return(FALSE)
  }
  
  # Get sample number from expression file
  txt_file <- gsub("\\.bin$", ".txt", peak_counts)
  if(file.exists(txt_file)) {
    sample_number <- ncol(read.table(txt_file, header=TRUE)) - 1
  } else {
    warning("Could not find text version of counts file. Using default sample number.")
    sample_number <- 21
  }
  
  # Set output directory with window size and PC value
  pc_results_dir <- file.path(output_dir, paste0("window_", window_kb, "kb"), paste0("pc", pc_value))
  pc_log_dir <- file.path(output_dir, "logs", paste0("window_", window_kb, "kb"), paste0("pc", pc_value))
  
  # Create PC-specific directories
  dir.create(pc_results_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(pc_log_dir, recursive=TRUE, showWarnings=FALSE)
  
  # Set output file name with PC value in filename
  output_suffix <- ifelse(run_permutation, "_permuted", "")
  output_file <- file.path(
    pc_results_dir,
    paste0(condition, "_", chr_name, "_pc", pc_value, "_", window_kb, "kb", output_suffix, ".txt")
  )
  
  # Set permutation flag
  perm_flag <- ifelse(run_permutation, "--random-permutation", "")
  
  # Create SBATCH command
  job_name <- paste0("RASQUAL_", condition, "_", chr_name, "_pc", pc_value)
  log_file <- file.path(pc_log_dir, paste0(job_name, output_suffix, ".log"))
  error_file <- file.path(pc_log_dir, paste0(job_name, output_suffix, ".err"))
  
sbatch_cmd <- paste0(
  "sbatch",
  " -p ", slurm_params$partition,
  " -n ", slurm_params$n_cores,
  " --mem ", slurm_params$mem_gb * 1000,
  " -t ", slurm_params$time,
  " -o ", log_file,
  " -e ", error_file,
  " -J ", job_name,
  " ", script_dir, "/03_run_multithread_rasqual_peakcenter.sbatch",
  " ", peak_info,
  " ", peak_counts,
  " ", sample_offsets,
  " ", covariates,
  " ", vcf_dir,
  " ", vcf_file,
  " ", sample_number,
  " ", output_file,
  " ", pc_value,
  " ", perm_flag  
)
  
  # Print command for debugging
  message("Submitting job for ", condition, " ", chr_name, " PC", pc_value, ":")
  message("  ", sbatch_cmd)
  
  # Submit job
  system(sbatch_cmd)
  
  # Small delay to prevent overwhelming scheduler
  Sys.sleep(1)
  
  return(TRUE)
}

#---------------------------------------------------------------
# Main execution
#---------------------------------------------------------------
message("=== RASQUAL Job Submission ===")
message("Window size: ", window_kb, "kb")
message("Conditions: ", paste(conditions, collapse=", "))
message("Chromosomes: ", paste(chromosomes, collapse=", "))
message("PC values: ", paste(pc_values, collapse=", "))
message("Permutations: ", ifelse(run_permutations, "YES", "NO"))

# Submit jobs with PC loop
job_count <- 0
for(pc_value in pc_values) {
  message("\n====== Processing PC", pc_value, " ======")
  
  # Update the PC value in parameters
  rasqual_params$p <- pc_value
  
  for(condition in conditions) {
    message("Processing condition: ", condition)
    
    for(chrom in chromosomes) {
      # Submit standard analysis
      success <- submit_rasqual_jobs(condition, chrom, window_kb, pc_value)
      if(success) job_count <- job_count + 1
      
      # Submit permutation analysis if requested
      if(run_permutations) {
        success <- submit_rasqual_jobs(condition, chrom, window_kb, pc_value, run_permutation = TRUE)
        if(success) job_count <- job_count + 1
      }
    }
  }
}

message("\n==== Job submission summary ====")
message("Total jobs submitted: ", job_count)
message("Results will be saved in: ", file.path(output_dir, paste0("window_", window_kb, "kb"), "pc[0-20]"))
message("Logs will be saved in: ", file.path(output_dir, "logs", paste0("window_", window_kb, "kb"), "pc[0-20]"))
message("Check job status with: squeue -u $USER")
