#RASQUAL permutation

library(dplyr)
library(stringr)
rm(list=ls())
options(stringsAsFactors=FALSE)

# Setting up variables for all conditions and chromosomes
path2vcf <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf"
conditions <- c("pbs", "fnf")
chromosomes <- c(1:22, "X")
#conditions <- "fnf"
#chromosomes <- 11

# Create a log file to track job submissions
log_file <- file.path("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/logs", 
                      paste0("rasqual_permutation_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

# Open the log file connection
log_conn <- file(log_file, open = "w")

# Write header to the log file
writeLines(c(
  "RASQUAL Permutation Job Submission Log",
  "====================================",
  paste("Log generated on:", Sys.time()),
  "Conditions: ", paste(conditions, collapse=", "),
  "Chromosomes: ", paste(chromosomes, collapse=", "),
  "----------------------------"
), log_conn)

# Loop through all conditions, chromosomes, and PCs
for(condition in conditions) {
  vcf_dir <- file.path(path2vcf, condition)
  
  for(i in chromosomes) {
    chr <- paste0("chr", i)
    
    for(pc in 1:10) {
      # Set up file paths
      peakinfo <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs/peak_info_", chr, ".txt")
      
      # Set counts and offsets based on condition (pbs->CTL, fnf->FNF)
      condition_label <- ifelse(condition == "pbs", "CTL", "FNF")
      peakcounts <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs/", chr, "_", condition_label, ".expression.bin")
      sampleoffsets <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs/", chr, "_", condition_label, ".size_factors.bin")
      
      # Covariate and VCF files
      covariates <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/covariates/", condition, "_covariates_pc", pc, ".bin")
      InputVCF <- paste0(chr, ".", condition, ".ASCounts.vcf.gz")
      
      # Get sample number
      samplenumber <- ncol(read.table(paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs/", chr, "_", condition_label, ".expression.txt"))) - 1
      
      # Create output directory for permutation results if it doesn't exist
      permutation_dir <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/permutation_results/", condition)
      dir.create(permutation_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Set up output file with permutation specification
      Output_File <- paste0(permutation_dir, "/", chr, "_", condition, "_PC", pc, "_RASQUALResults_permutation_25kb.txt")
      
      # IMPORTANT CHANGE: You need to pass the permutation flag as the 9th parameter
      # Create sbatch command with random-permutation flag
      commandline <- paste0("sbatch -p general -n 10 --mem 5000 -t 7-00:00:00 ",
                            "-o RASQUAL_perm_", condition, "_", chr, "_PC", pc, ".out ",
                            "-e RASQUAL_perm_", condition, "_", chr, "_PC", pc, ".err ",
                            "-J RASQUAL_perm_", condition, "_", chr, "_PC", pc, " ",
                            "/work/users/s/e/seyoun/CQTL_AI/scripts/rasqual/run_rasqual_perm.sbatch ",  # Use the new script
                            peakinfo, " ",
                            peakcounts, " ",
                            sampleoffsets, " ",
                            covariates, " ",
                            vcf_dir, " ",
                            InputVCF, " ",
                            samplenumber, " ",
                            Output_File, " ",
                            "--random-permutation")  # This is now properly passed as the 9th parameter
      
      # Execute the command
      job_submission_result <- system(commandline, intern = TRUE)
      
      # Write detailed log information
      log_entry <- paste(
        "----------------------------",
        paste("Timestamp:", Sys.time()),
        paste("Condition:", condition),
        paste("Chromosome:", chr),
        paste("PC:", pc),
        paste("Permutation: YES"),
        "Command Details:",
        commandline,
        "Job Submission Result:",
        paste(job_submission_result, collapse = "\n"),
        sep = "\n"
      )
      
      writeLines(log_entry, log_conn)
      
      # Print to console as well
      cat(log_entry, "\n")
    }
  }
}

# Close the log file connection
close(log_conn)

# Print the location of the log file
cat("Permutation job submission log saved to:", log_file, "\n")