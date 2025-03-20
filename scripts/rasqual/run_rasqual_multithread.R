library(dplyr)
library(stringr)

rm(list=ls())
options(stringsAsFactors=FALSE)

# Setting up variables for all conditions and chromosomes
path2vcf <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf"
conditions <- c("pbs", "fnf")
chromosomes <- c(1:22, "X")
#conditions <- "pbs"
#chromosomes <- 11

# Create a log file to track job submissions
log_file <- file.path("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/logs", paste0("rasqual100kb_job_submission_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

# Open the log file connection
log_conn <- file(log_file, open = "w")

# Write header to the log file
writeLines(c(
  "RASQUAL Job Submission Log",
  "=====================",
  paste("Log generated on:", Sys.time()),
  "Conditions: ", conditions,
  "Chromosomes: ", as.character(chromosomes),
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
      outdir <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results_100kb/", condition)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      # Set up output file with 1kb or 100kb or whatever you want  specification
      Output_File <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results_100kb/", 
                            condition, "/", chr, "_", condition, "_PC", pc, "_RASQUALResults_100kb.txt")
      
      
      # Create sbatch command with your specific structure
      commandline <- paste0("sbatch -p general -n 10 --mem 5000 -t 3-00:00:00 ",
                            "-o RASQUAL_", condition, "_", chr, ".out ",
                            "-e RASQUAL_", condition, "_", chr, ".err ",
                            "-J RASQUAL_", condition, "_", chr, " ",
                            "/work/users/s/e/seyoun/CQTL_AI/scripts/rasqual/run_rasqual_multithread_peakcenter.sbatch ",
                            peakinfo, " ",
                            peakcounts, " ",
                            sampleoffsets, " ",
                            covariates, " ",
                            vcf_dir, " ",
                            InputVCF, " ",
                            samplenumber, " ",
                            Output_File)
      
      # Execute the command
      job_submission_result <- system(commandline, intern = TRUE)
      
      # Write detailed log information
      log_entry <- paste(
        "----------------------------",
        paste("Timestamp:", Sys.time()),
        paste("Condition:", condition),
        paste("Chromosome:", chr),
        paste("PC:", pc),
        "Command Details:",
        commandline,
        "Job Submission Result:",
        paste(job_submission_result, collapse = "\n"),
        sep = "\n"
      )
      
      writeLines(log_entry, log_conn)
      
      # Optional: Print to console as well
      cat(log_entry, "\n")
    }
  }
}

# Close the log file connection
close(log_conn)

# Print the location of the log file
cat("Job submission log saved to:", log_file, "\n")
