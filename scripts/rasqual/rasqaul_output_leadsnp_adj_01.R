#Script for production of lead SNP (per unique peak) table for EigenMT file prep and run
#This input will be compared against raw EigenMT runs and should produce similar outcomes with tied pvalues being potentially different
##- multiple SNPs can be associated with the same genomic feature
##- finding the most significant one for each feature
##- when multiple SNPs have the same signifcant (p-value), I am planning to use their disance from the peak center as a tiebreaker
  
library(dplyr)
library(stringr)
library(data.table)

#-------------------------------------------------------------------------------
calculate_peak_distances <- function(feature_ids, snp_positions) {
  # Extract start and end positions using regular expressions
  positions <- stringr::str_match(feature_ids, ".*_(\\d+)_(\\d+)")[, 2:3]
  
  # Convert to numeric once
  starts <- as.numeric(positions[,1])
  ends <- as.numeric(positions[,2])
  
  # Vectorized calculation of peak centers and distances
  centers <- (starts + ends) %/% 2
  distances <- abs(centers - as.numeric(snp_positions))  # Convert SNP positions to numeric
  
  return(distances)
}

#------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
condition <- args[1]  # PC number and tissue type (e.g., fnf_PC1, pbs_PC1)

print(condition)


# Read input file
rasqaul_combined <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results/combined/", 
                           condition, "_combined.txt")
rasqual_dt <- fread(rasqaul_combined, header = TRUE)

# Calculate distances
rasqual_dt[, Distance_from_PeakCenter := calculate_peak_distances(Feature, `SNP Position`)]

# Select lead SNPs by:
# 1. Group by Feature
# 2. For each Feature: 
#    - Find SNP with lowest p-value
#    - If tie in p-value, use distance as tiebreaker
lead_snps <- rasqual_dt[, {
  if(.N == 1) {
    # If only one SNP, it's automatically the lead
    .SD
  } else {
    # For features with multiple SNPs:
    # First get the SNP(s) with minimum p-value
    min_pval <- min(PValue)
    min_pval_snps <- .SD[PValue == min_pval]
    
    if(nrow(min_pval_snps) == 1) {
      # If only one SNP has the minimum p-value, it's our lead SNP
      min_pval_snps
    } else {
      # If multiple SNPs have the same p-value,
      # select the one closest to peak center
      min_pval_snps[order(Distance_from_PeakCenter)][1]
    }
  }
}, by = Feature]


# Modified output path
output_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results/lead_snps/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

Outputpath <- paste0(output_dir, condition, "_lead_snps.txt")
write.table(lead_snps, Outputpath, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
