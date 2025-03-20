#Rasqual_prep_for_covariates
library(data.table)
library(dplyr)
library(tibble)
library(base)

#-------------------------------------------------------------------------------
# 1. Calculating the PCA
## Building DESeq matrix

ATACdds <- DESeqDataSetFromMatrix(countData = counts_macs2,
                                  colData = meta_final,
                                  rowRanges=macs2_gr,
                                  design = ~Donor)
normalizationFactors(ATACdds) = cqnNormFactors

vsd <- vst(ATACdds)
log2counts = assay(ATACdds)

## Calculating the PCA

pca_logcpm <- prcomp(t(log2counts)) #This transposition is just to have samples as rows and sites as columns, standard from PCA from what I can find reading around
#calculate percent of variance explained by each axis for this new PCA ordination
pca_perc = round((pca_logcpm$sdev)^2 / sum(pca_logcpm$sdev^2) *100,2);
# extract data of first 10 PCs from pca object
pca_points <- as.data.frame(pca_logcpm$x[,1:42]) 

#Cleanup genotype PCA-----------------------------------------------------------

geno_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/geno/"
pbs_pca <- fread(paste0(geno_dir,"pbs_geno/05.pca/cqtl.eigenvec"), data.table = FALSE)
fnf_pca <- fread(paste0(geno_dir,"fnf_geno/05.pca/cqtl.eigenvec"), data.table = FALSE)

output_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/covariates/"

process_rasqual_covariates <- function(pca_points, meta_final, geno_dir, output_dir, n_pcs=42) {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
  })
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  process_genotype_covariates <- function(geno_type) {
    # Read genotype PCA data
    geno_pca <- fread(paste0(geno_dir, geno_type, "_geno/05.pca/cqtl.eigenvec"), 
                      data.table = FALSE)
    
    # Process basic PCA data
    geno_pca_modified <- geno_pca %>%
      select(-V2) %>%
      rename(sampleID = V1) %>%
      rename_with(~ paste0("genopc", seq(1:20)), .cols = V3:V22)
    
    # Convert Sex to numeric (M=0, F=1)
    meta_final_processed <- meta_final %>%
      mutate(
        Sex_numeric = as.numeric(Sex == "F")
        # Convert date to numeric (days since first date)
        #ATACProtocolDate_numeric = as.numeric(as.Date(ATACProtocolDate, format="%m/%d/%y") - 
         #                                       min(as.Date(ATACProtocolDate, format="%m/%d/%y")))
      )
    
    # Add metadata (using numeric versions)
    base_covariates <- geno_pca_modified %>%
      select(sampleID, genopc1, genopc2, genopc3) %>%
      left_join(
        meta_final_processed %>% 
          select(sampleID, Sex_numeric),
        by = "sampleID"
      )
    
    # Process each PC increment
    for(i in 1:n_pcs) {
      message(sprintf("Processing %s PC%d...", geno_type, i))
      
      # Create current covariate set
      current_covariates <- base_covariates
      
      # Add PCs up to current i
      temp_pca <- pca_points[, 1:i, drop=FALSE]
      temp_pca$sampleID <- rownames(pca_points)
      
      current_covariates <- current_covariates %>%
        left_join(temp_pca, by = "sampleID")
      
      # Remove sampleID column for final output
      covariate_matrix <- as.matrix(current_covariates[,-1])
      
      # Save text format (for verification)
      txt_file <- file.path(output_dir, sprintf("%s_covariates_pc%d.txt", geno_type, i))
      write.table(covariate_matrix, 
                  txt_file, 
                  row.names = FALSE, 
                  col.names = FALSE,
                  quote = FALSE,
                  sep = "\t")
      
      # Convert to binary using writeBin
      bin_file <- file.path(output_dir, sprintf("%s_covariates_pc%d.bin", geno_type, i))
      fbin <- file(bin_file, "wb")
      writeBin(as.double(c(as.matrix(covariate_matrix))), fbin)
      close(fbin)
      
      # Print first few rows of current matrix for verification
      if(i == 1) {
        message("First few rows of covariate matrix:")
        print(head(covariate_matrix))
      }
    }
    
    return(base_covariates)
  }
  
  # Process both PBS and FNF
  message("Processing PBS covariates...")
  pbs_covariates <- process_genotype_covariates("pbs")
  
  message("Processing FNF covariates...")
  fnf_covariates <- process_genotype_covariates("fnf")
  
  return(list(
    pbs_base = pbs_covariates,
    fnf_base = fnf_covariates,
    output_dir = output_dir
  ))
}

results <- process_rasqual_covariates(
  pca_points = pca_points,
  meta_final = meta_final,
  geno_dir = "/work/users/s/e/seyoun/CQTL_AI/output/geno/",
  output_dir = "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/covariates",
  n_pcs = 42
)
