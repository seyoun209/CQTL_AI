# Script for adjusting tied p-values in combined RASQUAL results
library(dplyr)
library(stringr)
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
RASQUAL <- fread(combined_file)

Unique_Peak <- unique(RASQUAL$Feature)

Lead_table <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(RASQUAL)))
colnames(Lead_table) <- colnames(RASQUAL)

for (i in 1:length(Unique_Peak)) {
  print(i)
  
  RASQUAL_filtered <- filter(RASQUAL, Feature == Unique_Peak[i])
  
  if (nrow(RASQUAL_filtered) == 1) {
    Lead_SNP <- RASQUAL_filtered[1,]
  } else {
    RASQUAL_filtered <- RASQUAL_filtered[with(RASQUAL_filtered, order(PValue, Distance_From_Peak)),]
    RASQUAL_filtered$PValue[2:nrow(RASQUAL_filtered)] <- RASQUAL_filtered$PValue[2:nrow(RASQUAL_filtered)] + 1E-15
    Lead_SNP <- RASQUAL_filtered
  }
  
  Lead_table <- rbind(Lead_table, Lead_SNP)
}

Lead_table <- Lead_table[-1,]

# Write adjusted results
adjusted_file <- file.path(tied_adjusted_dir, paste0(condition, "_", pc, "_tied_adjusted.txt"))
fwrite(Lead_table, adjusted_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Also create eigenMT input file
eigenmt_input_dir <- file.path(base_dir, "output/eigenMT/tiedSNPs_adjusted_input")
dir.create(eigenmt_input_dir, recursive = TRUE, showWarnings = FALSE)

eigenmt_input <- Lead_table %>%
  mutate(
    snp = rs_ID,
    peak = Feature,
    statistic = as.numeric(Chi_square_statistic),
    `p-value` = PValue,
    FDR = PValue,
    beta = as.numeric(Effect_Size)
  ) %>%
  select(snp, peak, statistic, `p-value`, FDR, beta)

eigenmt_file <- file.path(eigenmt_input_dir, paste0(condition, "_", pc, "_eigenMT_input.txt"))
fwrite(eigenmt_input, eigenmt_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

print("Processing completed!")