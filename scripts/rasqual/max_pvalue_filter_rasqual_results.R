# Script for filtering RASQUAL results by eigenMT-derived p-value cutoffs
library(dplyr)
library(stringr)
library(data.table)


fnf_pc1 <- fread("/work/users/s/e/seyoun/CQTL_AI/output/eigenMT/input/fnf/gen.positions_chr1_fnf_PC1.txt")
              
rm(list=ls())
options(stringsAsFactors=FALSE)

# Define parameters
conditions <- c("fnf", "pbs")
pcs <- paste0("PC", 1:10)
base_dir <- "/work/users/s/e/seyoun/CQTL_AI"

# Create output directory
output_dir <- file.path(base_dir, "output/rasqual_qtl/results/pvalue_filtered")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Process each condition and PC
for (condition in conditions) {
  for (pc in pcs) {
    print(paste("Processing", condition, pc))
    
    # Load RASQUAL results
      rasqual_file <- file.path(base_dir, "output/rasqual_qtl/results/combined", 
                              paste0(condition, "_", pc, "_combined.txt"))
    
    if (!file.exists(rasqual_file)) {
      warning(paste("RASQUAL file not found:", rasqual_file))
      next
    }
    
    rasqual <- fread(rasqual_file)
    
    # Load eigenMT p-value cutoff
    cutoff_file <- file.path(base_dir, "output/eigenMT/combined_results", 
                             condition, pc,
                             paste0(condition, "_", pc, "_MaxPvalueCutoff_05.txt"))
    
    if (!file.exists(cutoff_file)) {
      warning(paste("Cutoff file not found:", cutoff_file))
      next
    }
    
    cutoff <- fread(cutoff_file)
    cutoff_value <- as.numeric(cutoff[1,1])
    
    if (is.na(cutoff_value)) {
      warning(paste("No significant hits for", condition, pc, "- skipping"))
      next
    }
    
    # Filter RASQUAL results by p-value cutoff
    pvalue_filtered <- filter(rasqual, PValue <= cutoff_value)
    
    # Write filtered results
    filtered_output_file <- file.path(output_dir, 
                                      paste0(condition, "_", pc, "_PValueFiltered_Results.txt"))
    
    write.csv(pvalue_filtered, filtered_output_file, row.names = FALSE)
    
    # Print summary
    print(paste("  Original rows:", nrow(rasqual), 
                "Filtered rows:", nrow(pvalue_filtered),
                "Cutoff p-value:", cutoff_value))
  }
}