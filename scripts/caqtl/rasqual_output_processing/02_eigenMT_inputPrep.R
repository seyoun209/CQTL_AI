#### filepath: /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/02_eigenMT_inputPrep.R
#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(stringr)

# Usage: Rscript 02_eigenMT_inputPrep.R <condition> <pc> <chromosome> [window_type]
if(length(commandArgs(trailingOnly = TRUE)) < 3){
  cat("Usage: Rscript 02_eigenMT_inputPrep.R <condition> <pc> <chromosome> [window_type]\n")
  cat("Example: Rscript 02_eigenMT_inputPrep.R pbs pc0 1\n")
  quit(status = 1)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
condition <- args[1]
pc <- args[2]    # Example: "pc0"
chrom <- args[3] # Chromosome, e.g., "1" or "chr1"
window_type <- ifelse(length(args) >= 4, args[4], "1kb")

# For file naming, ensure that the recode file part uses a "chr" prefix.
chrom_fname <- ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))

# Define base directories
base_dir <- "/work/users/s/e/seyoun/CQTL_AI"
phenotype_dir <- file.path(base_dir, "output/rasqual_output", paste0("window_", window_type), paste0("combined_", window_type))
genotype_dir <- file.path(base_dir, "output/caQTL_rasqual_input/vcf/recode_012")
output_base_dir <- file.path(base_dir, "output/rasqual_output/eigenMT/inputs", window_type)
# Create subdirectory for condition and pc.
output_dir <- file.path(output_base_dir, condition, tolower(pc))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set PC covariate string (if needed downstream)
pc_cov <- pc

cat(sprintf("Processing: Condition=%s, Chromosome=%s, Covariate=%s, Window=%s\n", 
            condition, chrom, pc_cov, window_type))

#############################
# 1. Construct gen.positions #
#############################
phenotype_file <- file.path(phenotype_dir, paste0(condition, "_", pc, "_", window_type, "_adjusted.txt"))
cat("Using phenotype file:", phenotype_file, "\n")
RASQUAL <- fread(phenotype_file)

gen.pos <- RASQUAL %>% 
  select(rs_ID, Chromosome, SNP_Position) %>% 
  filter(rs_ID != "SKIPPED") %>% 
  distinct(rs_ID, .keep_all = TRUE)
gen.pos$Chromosome <- str_remove(gen.pos$Chromosome, "^chr")
gen.pos_chr <- filter(gen.pos, Chromosome == as.character(chrom))
genpos_file <- file.path(output_dir, paste0("gen.positions_chr", chrom, "_", condition, "_", pc, ".txt"))
fwrite(gen.pos_chr, genpos_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Wrote gen.positions file:", genpos_file, "\n")

#############################
# 2. Construct genotypes      #
#############################
# Use chrom_fname to build the recode file name
geno_file <- file.path(genotype_dir, condition, paste0("recodeA_", chrom_fname, "_", condition, ".traw"))
cat("Using genotype file:", geno_file, "\n")
geno_dt <- fread(geno_file)
genotypes <- as.data.frame(geno_dt)[, c("SNP", colnames(geno_dt)[7:ncol(geno_dt)])]
genotypes <- genotypes %>% filter(SNP %in% gen.pos_chr$rs_ID)
geno_out_file <- file.path(output_dir, paste0("genotypes_chr", chrom, "_", condition, "_", pc, ".txt"))
write.table(genotypes, geno_out_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Wrote genotypes file:", geno_out_file, "\n")

#############################
# 3. Construct phe.positions  #
#############################
# Construction of Phe.position.txt for chromosome 3
phe.pos <- RASQUAL %>%
  select(peak_id = Feature, Chromosome) %>%
  mutate(
    chrom_probe = str_remove(Chromosome, "^chr"),
    tmp = str_remove(peak_id, paste0("^chr", chrom, "_")),
    s1 = str_split_fixed(tmp, "_", 2)[, 1],
    s2 = str_split_fixed(tmp, "_", 2)[, 2]
  ) %>%
  select(peak_id, chrom_probe, s1, s2) %>% distinct()
  Phe.pos <- phe.pos %>% filter(chrom_probe == chrom)
phepos_file <- file.path(output_dir, paste0("phe.positions_chr", chrom, "_", condition, "_", pc, ".txt"))
write.table(Phe.pos, phepos_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Wrote phe.positions file:", phepos_file, "\n")

#############################
# 4. Construct qtl file       #
#############################
qtl <- RASQUAL %>% select(rs_ID, Feature, PValue)
colnames(qtl) <- c("snp", "peak", "p-value")
qtl <- qtl %>% filter(snp != "SKIPPED") %>% filter(snp %in% gen.pos_chr$rs_ID)
qtl_file <- file.path(output_dir, paste0("qtl_chr", chrom, "_", condition, "_", pc, ".txt"))
write.table(qtl, qtl_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Wrote qtl file:", qtl_file, "\n")

cat("EigenMT input preparation completed.\n")