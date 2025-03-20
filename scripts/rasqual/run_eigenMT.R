#Submitting eigenMT Jobs

library(data.table)
library(dplyr)

# Define parameters
conditions <- c("fnf", "pbs")
chromosomes <- c(1:22)
PCs <- c(1:10)

# Set base directories
base_dir <- "/work/users/s/e/seyoun/CQTL_AI"
input_dir <- file.path(base_dir, "output/eigenMT/input")
output_dir <- file.path(base_dir, "output/eigenMT/results")
logs_dir <- file.path(output_dir, "logs")
script_path <- "/work/users/s/e/seyoun/CQTL_AI/scripts/rasqual/run_eigenMT.sh"

# Create logs directory
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

# Create subdirectories for each condition and PC
for (condition in conditions) {
  for (pc in PCs) {
    pc_dir <- file.path(output_dir, condition, paste0("PC", pc))
    dir.create(pc_dir, showWarnings = FALSE, recursive = TRUE)
  }
}

# Job counter
job_count <- 0

# Process each condition
for (condition in conditions) {
  cat("Processing condition:", condition, "\n")
  
  # Process each PC covariate
  for (pc in PCs) {
    pc_cov <- paste0("PC", pc)
    cat("Processing PC covariate:", pc_cov, "\n")
    
    # Process each chromosome
    for (chr in chromosomes) {
      cat("Submitting job for chromosome", chr, "\n")
      
      # Define input files
      qtl_file <- file.path(input_dir, condition, paste0("qtl_chr", chr, "_", condition, "_", pc_cov, ".txt"))
      gen_file <- file.path(input_dir, condition, paste0("genotypes_chr", chr, "_", condition, "_", pc_cov, ".txt"))
      genpos_file <- file.path(input_dir, condition, paste0("gen.positions_chr", chr, "_", condition, "_", pc_cov, ".txt"))
      phepos_file <- file.path(input_dir, condition, paste0("phe.positions_chr", chr, "_", condition, "_", pc_cov, ".txt"))
      
      # Define output file
      output_file <- file.path(output_dir, condition, pc_cov, paste0("eigenMT_chr", chr, "_", condition, "_", pc_cov, ".txt"))
      
      # Check if all required files exist
      if (!all(file.exists(qtl_file, gen_file, genpos_file, phepos_file))) {
        missing_files <- c()
        if (!file.exists(qtl_file)) missing_files <- c(missing_files, qtl_file)
        if (!file.exists(gen_file)) missing_files <- c(missing_files, gen_file)
        if (!file.exists(genpos_file)) missing_files <- c(missing_files, genpos_file)
        if (!file.exists(phepos_file)) missing_files <- c(missing_files, phepos_file)
        
        cat("Warning: Missing input files for chr", chr, "with", pc_cov, ":\n")
        cat(paste(" -", missing_files), sep="\n")
        cat("Skipping this job.\n")
        next
      }
      
      # Create job name and log files
      job_name <- paste0("eigen_", condition, "_chr", chr, "_", pc_cov)
      log_out <- file.path(logs_dir, paste0(job_name, "_%j.out"))
      log_err <- file.path(logs_dir, paste0(job_name, "_%j.err"))
      
      # Create sbatch command
      sbatch_cmd <- paste(
        "sbatch",
        "-p general",  # Partition (adjust as needed)
        "-c 1",
        "--mem 5G",  # Memory allocation
        "-t 12:00:00",  # Time limit (12 hours)
        paste0("-o ", log_out),  # Standard output file
        paste0("-e ", log_err),  # Standard error file
        paste0("-J ", job_name),  # Job name
        script_path,  # Your existing script
        qtl_file,     # Arguments
        gen_file,
        genpos_file,
        phepos_file,
        output_file,
        chr
      )
      
      # Submit job
      cat("Submitting job:", job_name, "\n")
      system(sbatch_cmd)
      
      # Increment job counter
      job_count <- job_count + 1
      
      # Add a small delay to avoid overwhelming the scheduler
      Sys.sleep(0.1)
    }
  }
}

cat("Total jobs submitted:", job_count, "\n")
cat("All eigenMT jobs have been submitted. Check job status with 'squeue'.\n")
