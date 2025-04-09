#!/usr/bin/env Rscript

# Function to create BAM list file
create_bam_list <- function(condition, bam_dir, output_dir) {
  # Map condition to BAM file pattern
  bam_condition <- if(condition == "pbs") "CTL" else "FNF"
  
  # Get BAM files
  bam_files <- list.files(
    path = bam_dir,
    pattern = paste0(".*_", bam_condition, ".*sorted_final.bam$"),
    full.names = TRUE
  )
  
  # Filter out Femur and replicate samples
  bam_files <- bam_files[!grepl("Femur|replicate", bam_files)]
  
  # Check if any BAM files were found
  if(length(bam_files) == 0) {
    stop(sprintf("No BAM files found for condition %s in directory %s", condition, bam_dir))
  }
  
  # Print summary
  message(sprintf("Found %d BAM files for condition %s", length(bam_files), condition))
  
  # Write BAM list
  bam_list_file <- file.path(output_dir, paste0(condition, "_bamlist.txt"))
  writeLines(bam_files, bam_list_file)
  message(sprintf("BAM list written to %s", bam_list_file))
  
  return(bam_list_file)
}

# Function to create donor list file from BAM filenames
create_donor_list <- function(condition, bam_dir, output_dir) {
  # Map condition to BAM file pattern
  bam_condition <- if(condition == "pbs") "CTL" else "FNF"
  
  # Get BAM files
  bam_files <- list.files(
    path = bam_dir,
    pattern = paste0(".*_", bam_condition, ".*sorted_final.bam$"),
    full.names = TRUE
  )
  
  # Extract donor IDs from BAM filenames
  # Adjust this regex pattern to match how your donor IDs appear in filenames
  # Example pattern: "/path/to/CQTL_AM7717_CTL_Ankle_1.sorted_final.bam"
  donor_pattern <- ".*/(CQTL_[A-Z0-9]+)_.*\\.bam$"
  donor_ids <- unique(gsub(donor_pattern, "\\1", bam_files))
  
  # Write donor list
  donor_list_file <- file.path(output_dir, paste0(condition, "_donor_list.txt"))
  writeLines(donor_ids, donor_list_file)
  
  message(sprintf("Found %d unique donors for condition %s", length(donor_ids), condition))
  message(sprintf("Donor list written to %s", donor_list_file))
  
  return(donor_list_file)
}

# Function to submit jobs
submit_rasqual_jobs <- function(
    conditions = c("pbs", "fnf"),
    chromosomes = c(1:22, "X"),
    bam_dir = "/work/users/s/e/seyoun/CQTL_AI/output/filtered/blk_filter",
    output_dir = "/work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/vcf",
    input_base = "/work/users/s/e/seyoun/CQTL_AI/output/geno",
    sbatch_script = "/work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_input_processing/01_rasqual_atac_vcf_prep.sbatch") {
  
  # Create output directories
  for(condition in conditions) {
    # Create main output directory
    dir.create(file.path(output_dir, condition), recursive = TRUE, showWarnings = FALSE)
    # Create condition-specific temp directory
    dir.create(file.path(output_dir, condition, "temp"), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Process each condition
  for(condition in conditions) {
    message(sprintf("Processing condition: %s", condition))
    
    # Create BAM list file
    bam_list_file <- create_bam_list(condition, bam_dir, output_dir)
    
    # Also create donor list file (useful for other scripts)
    donor_list_file <- create_donor_list(condition, bam_dir, output_dir)
    
    # Verify that input VCF exists before submitting jobs
    input_vcf <- file.path(input_base, paste0(condition, "_geno/04.final/finalFiltered.vcf.gz"))
    if(!file.exists(input_vcf)) {
      warning(sprintf("Input VCF not found: %s", input_vcf))
      warning("Skipping condition: ", condition)
      next
    }
    
    # Process each chromosome
    for(chr in chromosomes) {
      # Add 'chr' prefix to chromosome number
      chr_name <- paste0("chr", chr)
      message(sprintf("  Submitting job for chromosome %s", chr_name))
      
      # Set up paths
      output_vcf <- file.path(output_dir, condition,
                              sprintf("%s.%s.ASCounts.vcf.gz", chr_name, condition))
      # Use condition-specific temp directory
      temp_dir <- file.path(output_dir, condition, "temp")
      
      # Create and execute sbatch command
      cmd <- sprintf("sbatch %s %s %s %s %s %s",
                     sbatch_script,
                     input_vcf,
                     chr_name,  # Using chr_name with prefix
                     bam_list_file,
                     output_vcf,
                     temp_dir)
      
      # Print command for debugging
      message(sprintf("  Running command: %s", cmd))
      
      # Submit job
      system(cmd)
      
      # Small delay to prevent overwhelming the scheduler
      Sys.sleep(1)
    }
  }
}

# ========== MAIN EXECUTION ==========

# Set parameters (modify as needed)
params <- list(
  conditions = c("pbs", "fnf"),
  chromosomes = c(1:22, "X"),
  bam_dir = "/work/users/s/e/seyoun/CQTL_AI/output/filtered/blk_filter",
  output_dir = "/work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/vcf",
  input_base = "/work/users/s/e/seyoun/CQTL_AI/output/geno",
  sbatch_script = "/work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_input_processing/01_rasqual_atac_vcf_prep.sbatch"
)

# Print parameters
message("=== RASQUAL VCF Preparation ===")
message("Using parameters:")
for(param_name in names(params)) {
  param_value <- params[[param_name]]
  if(length(param_value) > 5) {
    message(sprintf("  %s: %s ... (%d total)", param_name, 
                    paste(param_value[1:5], collapse=", "), length(param_value)))
  } else {
    message(sprintf("  %s: %s", param_name, paste(param_value, collapse=", ")))
  }
}

# Execute the job submission
message("Starting job submission...")
do.call(submit_rasqual_jobs, params)

message("All RASQUAL VCF preparation jobs submitted.")
message("Check job status with: squeue -u $USER")
