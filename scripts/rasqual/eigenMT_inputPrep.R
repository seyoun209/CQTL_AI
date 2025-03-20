#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(dplyr)
library(stringr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if all required arguments are provided
if (length(args) < 3) {
  cat("Usage: Rscript eigenMT_inputprep.R <condition> <pc> <chromosome>\n")
  cat("Example: Rscript eigenMT_inputprep.R pbs 1 1\n")
  quit(status = 1)
}

# Parse arguments
condition <- args[1]
pc <- as.numeric(args[2])
c <- args[3]  # chromosome

# Set directories
base_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/eigenMT/tiedSNPs_adjusted_input"
output_base_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/eigenMT/input"
genotype_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf/recode_012"

# Create output directories
dir.create(file.path(output_base_dir, condition), showWarnings = FALSE, recursive = TRUE)

# Set PC covariate string
pc_cov <- paste0("PC", pc)

# Detailed logging
cat(sprintf("Processing: Condition=%s, Chromosome=%s, Covariate=%s\n", 
            condition, c, pc_cov))

# Read eigenMT input file
eigenmt_input_file <- file.path(base_dir, paste0(condition, "_", pc_cov, "_eigenMT_input.txt"))

if (!file.exists(eigenmt_input_file)) {
  cat(sprintf("Warning: EigenMT input file not found: %s\n", eigenmt_input_file))
  quit(status = 1)
}

cat(sprintf("Reading eigenMT input file: %s\n", eigenmt_input_file))
eigenmt_data <- fread(eigenmt_input_file, header = TRUE)

# Extract chromosome from peak name and filter for current chromosome
chr_pattern <- paste0("chr", c, "_")
eigenmt_data_chr <- eigenmt_data %>% 
  filter(grepl(chr_pattern, peak))

if (nrow(eigenmt_data_chr) == 0) {
  cat(sprintf("Warning: No data for chromosome %s in file %s\n", c, eigenmt_input_file))
  quit(status = 1)
}

# Write QTL file with only the required columns
qtl_data <- eigenmt_data_chr %>%
  select(snp, peak, `p-value`) 

qtl_output <- file.path(output_base_dir, condition, 
                        paste0("qtl_chr", c, "_", condition, "_", pc_cov, ".txt"))
write.table(qtl_data, qtl_output, quote = FALSE, row.names = FALSE)
cat(sprintf("Saved qtl file: %s\n", qtl_output))


# Extract SNP information for gen.positions.txt
gen_pos <- eigenmt_data_chr %>%
  distinct(snp) %>%
  mutate(
    parts = str_split_fixed(snp, ":", 4),
    chr_snp = str_replace(parts[,1], "chr", ""),
    pos = as.numeric(parts[,2])
  ) %>%
  select(snp, chr_snp, pos)

# Write gen.positions.txt
gen_position_output <- file.path(output_base_dir, condition, 
                                 paste0("gen.positions_chr", c, "_", condition, "_", pc_cov, ".txt"))
write.table(gen_pos, gen_position_output, quote = FALSE, row.names = FALSE)
cat(sprintf("Saved gen.positions file: %s\n", gen_position_output))

# Construction of Phe.position.txt
phe_pos <- eigenmt_data_chr %>% 
  distinct(peak) %>%
  mutate(
    chrom_probe = sub("chr([^_]+)_.*", "\\1", peak),
    coords = sub("chr[^_]+_(.*)", "\\1", peak),
    s1 = as.numeric(gsub("_.*$", "", coords)),  # Get start position
    s2 = as.numeric(gsub("^.*_", "", coords))   # Get end position
  ) %>%
  select(peak, chrom_probe, s1, s2) %>%
  dplyr::rename(peak_id = peak)

# Write phe.positions.txt
phe_position_output <- file.path(output_base_dir, condition, 
                                 paste0("phe.positions_chr", c, "_", condition, "_", pc_cov, ".txt"))
write.table(phe_pos, phe_position_output, quote = FALSE, row.names = FALSE)
cat(sprintf("Saved phe.positions file: %s\n", phe_position_output))

geno_file_path <- file.path(genotype_dir, condition, paste0("recodeA_chr", c, "_", condition, ".traw"))

if (!file.exists(geno_file_path)) {
  cat(sprintf("Warning: Genotype file not found: %s\n", geno_file_path))
  quit(status = 1)
}

cat(sprintf("Reading genotype file: %s\n", geno_file_path))

# Get the list of SNPs from the QTL file
qtl_snps <- unique(eigenmt_data_chr$snp)
cat(sprintf("Found %d unique SNPs in QTL file\n", length(qtl_snps)))

# Read the genotype file
genotype_data <- try(fread(geno_file_path))
if (inherits(genotype_data, "try-error")) {
  cat("Error reading genotype file. Please check the file format and path.\n")
  quit(status = 1)
}

# Create a new dataframe with the format needed for eigenMT
genotypes_out <- data.frame(
  ID = genotype_data$SNP,  # Use SNP column as ID
  stringsAsFactors = FALSE
)

# Add sample data columns (columns 7 and above from traw file)
for (col in names(genotype_data)[7:ncol(genotype_data)]) {
  # Clean up the sample names if needed
  sample_name <- col
  if (grepl("_CQTL_", col)) {
    sample_name <- strsplit(col, "_CQTL_")[[1]][1]
  }
  genotypes_out[[sample_name]] <- genotype_data[[col]]
}

# Filter to keep only SNPs that are in the QTL data
genotypes_out <- genotypes_out[genotypes_out$ID %in% qtl_snps, ]

# Check if we found any matching SNPs
if (nrow(genotypes_out) == 0) {
  cat("WARNING: No matching SNPs found in genotype file! eigenMT will fail.\n")
  cat("Sample QTL SNPs:", head(qtl_snps), "\n")
  cat("Sample genotype SNPs:", head(genotype_data$SNP), "\n")
  quit(status = 1)
}

# Save chromosome-specific genotype file
genotypes_output_file <- file.path(output_base_dir, condition, paste0("genotypes_chr", c, "_", condition, "_", pc_cov, ".txt"))
write.table(genotypes_out, genotypes_output_file, sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Saved matched genotype file: %s with %d rows\n", genotypes_output_file, nrow(genotypes_out)))