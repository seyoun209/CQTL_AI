#!/usr/bin/env Rscript
# Script for adjusting tied p-values in combined RASQUAL results (optimized version)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript eigenMT_leadsnpadjust.R <condition> <pc>")
}
condition <- args[1]
pc <- args[2]
print(paste("Processing condition:", condition, "PC:", pc))

base_dir <- "/work/users/s/e/seyoun/CQTL_AI"

# Set input and output directories
combined_dir <- file.path(base_dir, "output/rasqual_qtl/results/combined")
tied_adjusted_dir <- file.path(base_dir, "output/rasqual_qtl/results/tied_adjusted")
dir.create(tied_adjusted_dir, recursive = TRUE, showWarnings = FALSE)

# Read combined RASQUAL file
combined_file <- file.path(combined_dir, paste0(condition, "_", pc, "_combined.txt"))
if (!file.exists(combined_file)) {
  stop(paste("Combined file not found:", combined_file))
}
print(paste("Reading combined file:", combined_file))

# Read file with data.table for better performance
# Set stringsAsFactors=FALSE to maintain character strings
RASQUAL <- fread(combined_file, stringsAsFactors = FALSE)

# Get unique features
Unique_Peak <- unique(RASQUAL$Feature)
print(paste("Total unique peaks to process:", length(Unique_Peak)))

# Create a new empty data.table with the same structure as RASQUAL
# This ensures correct column types
Lead_table <- RASQUAL[0]

# Create a progress tracker
total_peaks <- length(Unique_Peak)
progress_step <- max(1, floor(total_peaks / 20))  # Report progress every 5%
start_time <- Sys.time()
last_report_time <- start_time

# Process peaks in batches to manage memory
batch_size <- 1000
for (batch_start in seq(1, total_peaks, by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, total_peaks)
  batch_peaks <- Unique_Peak[batch_start:batch_end]
  
  # Create a batch result table
  batch_results <- list()
  
  # Process this batch of peaks
  for (i in 1:length(batch_peaks)) {
    current_peak <- batch_peaks[i]
    global_i <- batch_start + i - 1
    
    # Create filtered subset using data.table syntax
    RASQUAL_filtered <- RASQUAL[Feature == current_peak, .SD]
    
    # Apply the adjustment logic
    if (nrow(RASQUAL_filtered) == 1) {
      # If only one row, no need to adjust
      Lead_SNP <- RASQUAL_filtered
    } else {
      # Order by PValue and Distance_From_Peak, then adjust p-values
      setorder(RASQUAL_filtered, PValue, Distance_From_Peak)
      RASQUAL_filtered$PValue[2:nrow(RASQUAL_filtered)] <- 
        RASQUAL_filtered$PValue[2:nrow(RASQUAL_filtered)] + 1E-15
      Lead_SNP <- RASQUAL_filtered
    }
    
    # Add to batch results
    batch_results[[i]] <- Lead_SNP
    
    # Report progress periodically
    if (global_i %% progress_step == 0 || global_i == total_peaks) {
      current_time <- Sys.time()
      elapsed <- difftime(current_time, start_time, units = "mins")
      time_since_last_report <- difftime(current_time, last_report_time, units = "secs")
      
      # Only show progress updates if at least 30 seconds have passed
      if (time_since_last_report > 30 || global_i == total_peaks) {
        percent_done <- (global_i / total_peaks) * 100
        est_total_time <- elapsed / (percent_done / 100)
        est_remaining <- est_total_time - elapsed
        
        print(sprintf("Processed %d/%d peaks (%.1f%%) - Elapsed: %.1f min, Est. remaining: %.1f min", 
                      global_i, total_peaks, percent_done, 
                      as.numeric(elapsed), as.numeric(est_remaining)))
        
        last_report_time <- current_time
      }
    }
  }
  
  # Combine batch results and append to main results
  if (length(batch_results) > 0) {
    batch_table <- rbindlist(batch_results, use.names = TRUE)
    Lead_table <- rbindlist(list(Lead_table, batch_table), use.names = TRUE)
  }
  
  # Print a sample of the data to check
  if (batch_start == 1) {
    print("Sample of processed data (first 3 rows):")
    print(head(Lead_table, 3))
  }
  
  # Garbage collection to free memory after each batch
  rm(batch_results, batch_table)
  gc()
}

# Write adjusted results
adjusted_file <- file.path(tied_adjusted_dir, paste0(condition, "_", pc, "_tied_adjusted.txt"))
print(paste("Writing results to:", adjusted_file))
fwrite(Lead_table, adjusted_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Create eigenMT input file
eigenmt_input_dir <- file.path(base_dir, "output/eigenMT/tiedSNPs_adjusted_input")
dir.create(eigenmt_input_dir, recursive = TRUE, showWarnings = FALSE)

# Use data.table syntax for better performance
eigenmt_input <- Lead_table[, .(
  snp = rs_ID,
  peak = Feature,
  statistic = as.numeric(Chi_square_statistic),
  `p-value` = PValue,
  FDR = PValue,
  beta = as.numeric(Effect_Size)
)]

eigenmt_file <- file.path(eigenmt_input_dir, paste0(condition, "_", pc, "_eigenMT_input.txt"))
print(paste("Writing eigenMT input to:", eigenmt_file))
fwrite(eigenmt_input, eigenmt_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")
print(paste("Processing completed! Total time:", round(total_time, 2), "minutes"))