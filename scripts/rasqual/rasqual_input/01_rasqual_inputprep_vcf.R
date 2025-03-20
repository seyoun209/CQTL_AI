#!/usr/bin/env Rscript



## run this Rscript rasqual_inputprep_vcf.R


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
  
  # Write BAM list
  bam_list_file <- file.path(output_dir, paste0(condition, "_bamlist.txt"))
  writeLines(bam_files, bam_list_file)
  return(bam_list_file)
}

# Function to submit jobs
submit_rasqual_jobs <- function(
    conditions = c("pbs", "fnf"),
    chromosomes = c(1:22, "X"),
    bam_dir = "/work/users/s/e/seyoun/CQTL_AI/output/wasp/blk_filter",
    output_dir = "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf",
    input_base = "/work/users/s/e/seyoun/CQTL_AI/output/geno",
    sbatch_script = "/work/users/s/e/seyoun/CQTL_AI/scripts/rasqual/rasqual_inputprep_vcf.sbatch") {
  
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
    
    # Process each chromosome
    for(chr in chromosomes) {
      # Add 'chr' prefix to chromosome number
      chr_name <- paste0("chr", chr)
      message(sprintf("  Submitting job for chromosome %s", chr_name))
      
      # Set up paths
      input_vcf <- file.path(input_base, 
                             paste0(condition, "_geno/04.final/finalFiltered_wCHR.vcf.gz"))
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
      
      # Submit job
      system(cmd)
      
      # Small delay to prevent overwhelming the scheduler
      Sys.sleep(1)
    }
  }
}

# Run the submission
submit_rasqual_jobs(
  conditions = c("pbs", "fnf"),
  chromosomes = c(1:22, "X"),
  bam_dir = "/work/users/s/e/seyoun/CQTL_AI/output/wasp/blk_filter",
  output_dir = "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf",
  input_base = "/work/users/s/e/seyoun/CQTL_AI/output/geno",
  sbatch_script = "/work/users/s/e/seyoun/CQTL_AI/scripts/rasqual/rasqual_inputprep_vcf.sbatch"
)