#Rasqual script to calculate PCA
library(data.table)
library(PCAForQTL)
library(tidyverse)
library(dplyr)
#-------------------------------------------------------------------------------
# sample and the normfactors
load('/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/QC/QTL_GCContent.RData') # this will give two file: the meta_final and GCcontent
load("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/atac_chon_deseq.RData")
#GC content
load('/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/QC/QTL_GCContent.RData')

#CQN Normalized factors (cqnNormFactors, mapped_reads_ordered, counts_macs2, macs2_gr,macs2_filtered)
load("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/CQTL_cqnNormFactors.RData")


#DESeq2 size factors30 to adjust for library size and the gcCor.R script provided with RASQUAL
normalized_counts <- counts(chonATACdds, normalized=TRUE) # This is already normalized size factor

gc_corrected <- gcCor(Y = normalized_counts, 
                      gcvec = GCcontent, 
                      PLOT = T)

# Overall pattern Looks good to me! Here is the reasons:
## 1. relative flat in the middle GC content range (0.4-0.6) and most correction factors stay within -0.01 to 0.01 range.


#PCA (Use for the PCAforQTL-----------------------------------------------------
#Building DESeq matrix

dds = DESeqDataSetFromMatrix(countData = counts_macs2, colData = meta_final, design = ~DonorID) 
normalizationFactors(dds) = cqnNormFactors

vsd <- vst(dds)
log2counts = assay(vsd)

pca_logcpm <- prcomp(t(log2counts)) #This transposition is just to have samples as rows and sites as columns, standard from PCA from what I can find reading around
#calculate percent of variance explained by each axis for this new PCA ordination
pca_perc = round((pca_logcpm$sdev)^2 / sum(pca_logcpm$sdev^2) *100,2);
# extract data of first 10 PCs from pca object
pca_points <- as.data.frame(pca_logcpm$x[,1:42]) #NOTE: PC1 already explains the 43.41% and maybe PC1-6?


#PCA for the qtl ---------------------------------------------------------------
geno_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/geno/"
pbs_pca <- fread(paste0(geno_dir,"pbs_geno/05.pca/cqtl.eigenvec"), data.table = FALSE)
fnf_pca <- fread(paste0(geno_dir,"fnf_geno/05.pca/cqtl.eigenvec"), data.table = FALSE)

pbs_pca_modified <- pbs_pca %>%
  select(-V2) %>%
  rename(sampleID = V1) %>%
  rename_with(~ paste0("genopc", seq(1:20)), .cols = V3:V22)

#add the meta file (sex and the ATACprotocolDate)
pbs_temp_pca_cov <- pbs_pca_modified %>%
  left_join(
    meta_final %>% select(sampleID, Sex, ATACProtocolDate),
    by = "sampleID"
  )

# Add rownames from pca_points for each PC
for(i in 1:20) {
  pbs_cov_final <- pbs_temp_pca_cov %>%
    left_join(
      rownames_to_column(pca_points[,i, drop=FALSE], var="sampleID") %>%
        rename(!!pc_names[i] := 2),
      by="sampleID"
    )
}

"/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/covariates/"
