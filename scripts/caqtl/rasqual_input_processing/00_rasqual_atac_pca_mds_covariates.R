# Rasqual_prep for covariates
# edited by: Jess
# date: 2023-03-30
setwd("/work/users/s/e/seyoun/CQTL_AI/output/")
#---------------------------------------------------------------
# Load libraries
library(data.table)
library(dplyr)
library(tibble)
library(base)
library(deseq2)
source("/work/users/s/e/seyoun/CQTL_AI/scripts/ATAC_Diff/Correlation_functions.R")

load("diff_deseq2/condtion/QC/Chon_macs2_se.RData")
load("diff_deseq2/QC/CQN_results.RData")

geno_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/geno/"
pbs_pca <- fread(paste0(geno_dir,"pbs_geno/05.pca/cqtl.eigenvec"), data.table = FALSE)
fnf_pca <- fread(paste0(geno_dir,"fnf_geno/05.pca/cqtl.eigenvec"), data.table = FALSE)

output_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/caQTL_rasqual_input/covar"
n_pcs_max <- 20
#----------------------------------------------------------------
# PCA for rasqual covariates
dds <- DESeqDataSetFromMatrix(countData = counts_macs2,
                                  colData = meta_final,
                                  rowRanges=macs2_gr,
                                  design = ~Donor)
normalizationFactors(dds) = cqnNormFactors

# Get VST transformed data 
vsd_donor <- vst(dds, blind=TRUE)
vst_counts = assay(vsd_donor)  

# Calculate PCA - this step was missing!
pca_rasqual <- prcomp(t(vst_counts))

# Correlation plot
# Create PCA data frame
pca_df <- as.data.frame(pca_rasqual$x)
pca_vsd_df <- merge(pca_df, colData(vsd_donor), by = "row.names", all.x = TRUE)
ctl_plot <- create_correlation_plot(pca_vsd_df, condition = "CTL", title = "CTL Correlations")
fnf_plot <- create_correlation_plot(pca_vsd_df, condition = "FNF", title = "FNF Correlations")

# Calculate variance explained
var_explained <- pca_rasqual$sdev^2 / sum(pca_rasqual$sdev^2)
message("Variance explained by first 42 PCs: ", 
        paste(round(var_explained[1:42]*100, 2), collapse="%, "), "%")
pca_points <- as.data.frame(pca_rasqual$x[,1:42])

#---------------------------------------------------------------
# 2. Read genotype PCA files
message("Reading genotype PCA files...")

# Read PBS genotype PCA
pbs_geno_file <- paste0(geno_dir, "pbs_geno/05.pca/cqtl.eigenvec")
pbs_geno_pca <- fread(pbs_geno_file, data.table = FALSE)
pbs_geno_pca <- pbs_geno_pca %>%
  select(-V2) %>%
  rename(sampleID = V1) %>%
  rename_with(~ paste0("genoPCA", 1:20), .cols = V3:V22)

# Read FNF genotype PCA
fnf_geno_file <- paste0(geno_dir, "fnf_geno/05.pca/cqtl.eigenvec")
fnf_geno_pca <- fread(fnf_geno_file, data.table = FALSE)
fnf_geno_pca <- fnf_geno_pca %>%
  select(-V2) %>%
  rename(sampleID = V1) %>%
  rename_with(~ paste0("genoPCA", 1:20), .cols = V3:V22)

#---------------------------------------------------------------
# 3. Process metadata
message("Preparing metadata...")

# Add numeric versions of key covariates
meta_enhanced <- meta_final %>%
  mutate(
    Sex_numeric = as.numeric(Sex == "M"),  # 1=Male, 0=Female
    Protocol_batch = as.numeric(factor(ATACProtocolDate))  # Protocol batch as numeric factor
  )

#---------------------------------------------------------------
# 4. Generate covariates for PBS (CTL condition)
message("Creating PBS covariates...")

# Get CTL samples
ctl_indices <- which(meta_enhanced$Condition == "CTL")
ctl_samples <- meta_enhanced[ctl_indices, ]
message("Number of CTL samples: ", length(ctl_indices))

# Extract sample IDs for CTL
ctl_sample_ids <- ctl_samples$sampleID

# Filter genotype PCA data to only include CTL samples
pbs_geno_subset <- pbs_geno_pca[pbs_geno_pca$sampleID %in% ctl_sample_ids, ]
pbs_geno_subset <- pbs_geno_subset[match(ctl_sample_ids, pbs_geno_subset$sampleID), ]

# Extract ATAC PCs for CTL samples
ctl_pc_data <- pca_points[ctl_indices, ]

# Generate covariates for different PC counts
for(i in 0:n_pcs_max) {
  message("Creating PBS covariates with ", i, " PCs...")
  
  # Start with fixed covariates
  current_covars <- data.frame(
    # Add sex and protocol batch
    Sex = ctl_samples$Sex_numeric,
    Protocol = ctl_samples$Protocol_batch,
    # Add first 3 genotype PCs
    genoPCA1 = pbs_geno_subset$genoPCA1,
    genoPCA2 = pbs_geno_subset$genoPCA2,
    genoPCA3 = pbs_geno_subset$genoPCA3
  )
  
  # Add ATAC PCs if i > 0
  if(i > 0) {
    # Add ATAC PCs (already calculated in pca_points)
    atac_pcs <- ctl_pc_data[, 1:i, drop=FALSE]
    current_covars <- cbind(current_covars, atac_pcs)
  }
  
  # Save as text file (for verification)
  txt_file <- file.path(output_dir, sprintf("pbs_covariates_pc%d.txt", i))
  write.table(current_covars, 
              txt_file, 
              row.names = FALSE, 
              col.names = TRUE,
              quote = FALSE, 
              sep = "\t")
  
  # Save as binary (for RASQUAL)
  bin_file <- file.path(output_dir, sprintf("pbs_covariates_pc%d.bin", i))
  fbin <- file(bin_file, "wb")
  writeBin(as.double(c(as.matrix(current_covars))), fbin)
  close(fbin)
}

#---------------------------------------------------------------
# 5. Generate covariates for FNF
message("Creating FNF covariates...")

# Get FNF samples
fnf_indices <- which(meta_enhanced$Condition == "FNF")
fnf_samples <- meta_enhanced[fnf_indices, ]
message("Number of FNF samples: ", length(fnf_indices))

# Extract sample IDs for FNF
fnf_sample_ids <- fnf_samples$sampleID

# Filter genotype PCA data to only include FNF samples
fnf_geno_subset <- fnf_geno_pca[fnf_geno_pca$sampleID %in% fnf_sample_ids, ]
fnf_geno_subset <- fnf_geno_subset[match(fnf_sample_ids, fnf_geno_subset$sampleID), ]

# Extract ATAC PCs for FNF samples
fnf_pc_data <- pca_points[fnf_indices, ]

# Generate covariates for different PC counts
for(i in 0:n_pcs_max) {
  message("Creating FNF covariates with ", i, " PCs...")
  
  # Start with fixed covariates
  current_covars <- data.frame(
    # Add sex and protocol batch
    Sex = fnf_samples$Sex_numeric,
    Protocol = fnf_samples$Protocol_batch,
    # Add first 3 genotype PCs
    genoPCA1 = fnf_geno_subset$genoPCA1,
    genoPCA2 = fnf_geno_subset$genoPCA2,
    genoPCA3 = fnf_geno_subset$genoPCA3
  )
  
  # Add ATAC PCs if i > 0
  if(i > 0) {
    # Add ATAC PCs (already calculated in pca_points)
    atac_pcs <- fnf_pc_data[, 1:i, drop=FALSE]
    current_covars <- cbind(current_covars, atac_pcs)
  }
  
  # Save as text file (for verification)
  txt_file <- file.path(output_dir, sprintf("fnf_covariates_pc%d.txt", i))
  write.table(current_covars, 
              txt_file, 
              row.names = FALSE, 
              col.names = TRUE,
              quote = FALSE, 
              sep = "\t")
  
  # Save as binary (for RASQUAL)
  bin_file <- file.path(output_dir, sprintf("fnf_covariates_pc%d.bin", i))
  fbin <- file(bin_file, "wb")
  writeBin(as.double(c(as.matrix(current_covars))), fbin)
  close(fbin)
}

# Save mapping information for reference
sample_mapping <- data.frame(
  condition = c(rep("CTL", length(ctl_indices)), rep("FNF", length(fnf_indices))),
  sample_id = c(ctl_sample_ids, fnf_sample_ids),
  sample_index = c(ctl_indices, fnf_indices)
)

write.table(sample_mapping,
            file.path(output_dir, "sample_mapping.txt"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

message("All covariates created successfully. Files saved to: ", output_dir)
message("File naming: pbs_covariates_pc0.bin through pbs_covariates_pc20.bin")
message("File naming: fnf_covariates_pc0.bin through fnf_covariates_pc20.bin")
