#rasqual count and the offset

library(rasqualTools)
library(dplyr)
library(stringr)
#---------------------------------------------------------------

setwd("/work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/input_count")

# Set global parameters
output_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/input_count"
chromosomes <- paste0("chr", c(1:22, "X"))

message("Loading data...")
# Load CQN normalized ATAC counts
load("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condtion/QC/CQN_results.RData")
load("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condtion/QC/Chon_macs2_se.RData")

# Create peak info with ID
message("Creating peak info with IDs...")
macs2_filtered <- data.frame(
  Chr = as.character(seqnames(macs2_gr)),
  Start = start(macs2_gr),
  End = end(macs2_gr),
  Geneid =macs2_gr$peakID
)

# Create ID format
macs2_filtered$ID <- paste(macs2_filtered$Chr, macs2_filtered$Start, macs2_filtered$End, sep="_")

# Identify CTL and FNF columns
message("Splitting by condition...")
ctl_cols <- grep("CTL", colnames(counts_macs2), value=TRUE)
fnf_cols <- grep("FNF", colnames(counts_macs2), value=TRUE)

message(paste0("Found ", length(ctl_cols), " CTL samples"))
message(paste0("Found ", length(fnf_cols), " FNF samples"))

# Create CTL file
ctl_count <- macs2_filtered 
ctl_count <- cbind(ctl_count, counts_macs2[, ctl_cols])

# Create FNF file
fnf_count <- macs2_filtered
fnf_count <- cbind(fnf_count, counts_macs2[, fnf_cols])

# Write the full count tables for reference
write.table(ctl_count, file.path(output_dir, "peaks_CTL.txt"), 
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(fnf_count, file.path(output_dir, "peaks_FNF.txt"), 
            row.names=FALSE, quote=FALSE, sep="\t")

#---------------------------------------------------------------
# Create expression count files by chromosome
#---------------------------------------------------------------
message("Creating count files by chromosome...")

# Process CTL counts
message("Processing CTL counts...")
for (i in chromosomes) {
  message(paste0("  Processing ", i, " for CTL..."))
  
  # Filter by chromosome
  chr_specific <- subset(ctl_count, ctl_count$Chr == i)
  
  # Extract only count columns (skip the first 5 columns with metadata)
  first_count_col <- which(colnames(chr_specific) == ctl_cols[1])
  counts_data <- chr_specific[, first_count_col:ncol(chr_specific)]
  
  # Set row names to peak IDs
  row.names(counts_data) <- chr_specific$ID
  
  # Ensure all columns are numeric
  counts_data[] <- lapply(counts_data, function(x) as.numeric(x))
  
  # Create identifier and save matrix
  ID <- paste0(i, "_pbs")  # Using pbs instead of CTL for file naming
  countslist <- list(counts_data)
  names(countslist) <- ID
  
  saveRasqualMatrices(countslist, output_dir, file_suffix = "expression")
}

# Process FNF counts
message("Processing FNF counts...")
for (i in chromosomes) {
  message(paste0("  Processing ", i, " for FNF..."))
  
  # Filter by chromosome
  chr_specific <- subset(fnf_count, fnf_count$Chr == i)
  
  # Extract only count columns (skip the first 5 columns with metadata)
  first_count_col <- which(colnames(chr_specific) == fnf_cols[1])
  counts_data <- chr_specific[, first_count_col:ncol(chr_specific)]
  
  # Set row names to peak IDs
  row.names(counts_data) <- chr_specific$ID
  
  # Ensure all columns are numeric
  counts_data[] <- lapply(counts_data, function(x) as.numeric(x))
  
  # Create identifier and save matrix
  ID <- paste0(i, "_fnf")
  countslist <- list(counts_data)
  names(countslist) <- ID
  
  saveRasqualMatrices(countslist, output_dir, file_suffix = "expression")
}

#---------------------------------------------------------------
# Create size factor files for offsets
#---------------------------------------------------------------
message("Creating size factor files...")

# Ensure column ordering matches
message("Checking column ordering...")
ctl_cqn_cols <- intersect(ctl_cols, colnames(cqnNormFactors))
fnf_cqn_cols <- intersect(fnf_cols, colnames(cqnNormFactors))

# Check if all samples have normalization factors
if (length(ctl_cqn_cols) != length(ctl_cols)) {
  warning("Not all CTL samples have normalization factors")
}
if (length(fnf_cqn_cols) != length(fnf_cols)) {
  warning("Not all FNF samples have normalization factors")
}

# Split normalization factors by condition
CQN_CTL <- cqnNormFactors[, ctl_cqn_cols]
CQN_FNF <- cqnNormFactors[, fnf_cqn_cols]

# Process CTL normalization factors
message("Processing CTL normalization factors...")
CQN_CTL_df <- as.data.frame(CQN_CTL)
CQN_CTL_df$Chr <- macs2_filtered$Chr
CQN_CTL_df$PeakID <- macs2_filtered$ID

# Process FNF normalization factors
message("Processing FNF normalization factors...")
CQN_FNF_df <- as.data.frame(CQN_FNF)
CQN_FNF_df$Chr <- macs2_filtered$Chr
CQN_FNF_df$PeakID <- macs2_filtered$ID

# Process by chromosome for CTL
message("Creating size factor files for CTL...")
for (i in chromosomes) {
  message(paste0("  Processing ", i, " for CTL..."))
  
  # Filter by chromosome
  chr_specific <- subset(CQN_CTL_df, CQN_CTL_df$Chr == i)
  
  # Extract only normalization factor columns
  norm_cols <- setdiff(colnames(chr_specific), c("Chr", "PeakID"))
  normfactors <- chr_specific[, norm_cols]
  
  # Set row names to peak IDs
  row.names(normfactors) <- chr_specific$PeakID
  
  # Ensure all columns are numeric
  normfactors[] <- lapply(normfactors, function(x) as.numeric(x))
  
  # Create identifier and save matrix
  ID <- paste0(i, "_pbs")  # Using pbs instead of CTL for file naming
  normlist <- list(normfactors)
  names(normlist) <- ID
  
  saveRasqualMatrices(normlist, output_dir, file_suffix = "size_factors")
}

# Process by chromosome for FNF
message("Creating size factor files for FNF...")
for (i in chromosomes) {
  message(paste0("  Processing ", i, " for FNF..."))
  
  # Filter by chromosome
  chr_specific <- subset(CQN_FNF_df, CQN_FNF_df$Chr == i)
  
  # Extract only normalization factor columns
  norm_cols <- setdiff(colnames(chr_specific), c("Chr", "PeakID"))
  normfactors <- chr_specific[, norm_cols]
  
  # Set row names to peak IDs
  row.names(normfactors) <- chr_specific$PeakID
  
  # Ensure all columns are numeric
  normfactors[] <- lapply(normfactors, function(x) as.numeric(x))
  
  # Create identifier and save matrix
  ID <- paste0(i, "_fnf")
  normlist <- list(normfactors)
  names(normlist) <- ID
  
  saveRasqualMatrices(normlist, output_dir, file_suffix = "size_factors")
}

#---------------------------------------------------------------
# Create peak info files
#---------------------------------------------------------------
message("Creating peak info files...")

# Create standard peak info files
for (i in chromosomes) {
  message(paste0("  Creating peak info for ", i))
  
  # Filter by chromosome
  chr_specific <- subset(macs2_filtered, macs2_filtered$Chr == i)
  
  # Create peak info data frame
  peak_file <- data.frame(
    PeakID = chr_specific$ID,
    Chr = chr_specific$Chr,
    Start = chr_specific$Start,
    End = chr_specific$End
  )
  
  # Write to file
  file_name <- file.path(output_dir, paste0("peak_info_", i, ".txt"))
  write.table(peak_file, file = file_name, 
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Create peak info files with gene IDs
message("Creating peak info files with gene IDs...")
peak_info_dir <- file.path(dirname(output_dir), "peak_info")
dir.create(peak_info_dir, showWarnings = FALSE, recursive = TRUE)

for (i in chromosomes) {
  # Filter by chromosome
  chr_specific <- subset(macs2_filtered, macs2_filtered$Chr == i)
  
  # Create peak info data frame with gene ID
  peak_file <- data.frame(
    GeneID = chr_specific$Geneid,
    PeakID = chr_specific$ID,
    Chr = chr_specific$Chr,
    Start = chr_specific$Start,
    End = chr_specific$End
  )
  
  # Write to file
  file_name <- file.path(peak_info_dir, paste0("peak_info_", i, ".txt"))
  write.table(peak_file, file = file_name, 
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}

message("Done! All files have been created in ", output_dir)