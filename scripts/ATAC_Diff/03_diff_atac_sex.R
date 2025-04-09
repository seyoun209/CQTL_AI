# Differential analysis atac-seq - Sex
##  Building deseq2 matrix
# Load libraries
library(data.table)
library(dplyr)
library(tibble)
library(base)
library(cqn)
library(rasqualTools)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DESeq2)
library(httpgd)
library(ggplot2)
library(plotgardener)

setwd("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/sex")
#---------------------------------------------------------------
load("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/Chon_macs2_se.RData") # load the data
load("//work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/CQN_results.RData")
print(ls())

meta_final$Condition <- factor(meta_final$Condition, levels = c("CTL", "FNF"))
meta_final$Sex <- factor(meta_final$Sex, levels = c("M", "F"))


sexdds <- DESeqDataSetFromMatrix(countData = counts_macs2,
                                 colData = meta_final,
                                 rowRanges = macs2_gr,
                                 design = ~ Condition + Sex)

# Apply CQN normalization factors
normalizationFactors(sexdds) <- cqnNormFactors
# Run DESeq2
sexdds <- DESeq(sexdds)

# Get results
sex_res <- results(sexdds)
print(summary(sex_res,alpha=0.05, na.rm=TRUE))
sex_res_Shrink <- lfcShrink(sexdds, coef="Sex_F_vs_M", type="apeglm",format =  "GRanges")
sex_res_Shrink$peakID <- rowData(sexdds)$peakID

save(sexdds, sex_res, sex_res_Shrink,
     file="atac_cqnnorm_sex_deseq2_re.RData")


#----------------------------------------------------------------
# 1. PCA plots

normCountsSex <- counts(sexdds, normalized = TRUE)
rldSex <- rlog(sexdds, blind=TRUE)
vsdSex <- vst(sexdds, blind=TRUE)
save(vsdSex, rldSex, normCountsSex,
 file = "QC/sexATAC_deseq2_normalized.RData")

# Calculate PCA on vst-transformed data
pca_result <- prcomp(t(assay(vsdSex)))
plotPCA(vsdSex, intgroup = c("Sex"), ntop=1000) +
  ggplot2::theme(aspect.ratio = 1)
# Create PCA data frame
pca_df <- as.data.frame(pca_result$x)
variance_explained <- summary(pca_result)$importance[2, ]
pca_df$Donor <- colData(vsdSex)$Donor
pca_df$Sex <-colData(vsdSex)$Sex


# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_var <- cumsum(var_explained)

pca_plot_vsdSex <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Sex)) +
  geom_point(size = 3) +  # Increase point size here
  #geom_text_repel(size = 5, box.padding = unit(1, "lines")) +
  scale_color_manual(values = c("M" = "#c864da", "F" = "#55e0bd")) +
  labs(
    title = "Total peaks (n=281,227) with all samples",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))


# Save PCA plot
save(pca_plot_vsdSex,
file = "plots/PCA_sexdeseq2_plot_vsd.rds")

#------------------------------------------------------

## Finding the significant peaks
atac_res_Shrink_noNA <- atac_res_Shrink[!is.na(atac_res_Shrink$padj)]
atac_res_Shrink_df <- as.data.frame(atac_res_Shrink_noNA)
#Getting significant ATAC-signals
diff_atac_sig <- atac_res_Shrink_df %>%
  dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.05)

#finiding up sig and down sig
gained <- diff_atac_sig[diff_atac_sig$log2FoldChange > 0 & diff_atac_sig$padj < 0.05,]
lost <- diff_atac_sig[diff_atac_sig$log2FoldChange < 0 & diff_atac_sig$padj < 0.05,]
static <- atac_res_Shrink_df[!(atac_res_Shrink_df$peakID %in% c(gained$peakID, lost$peakID)), ]


atac_res_Shrink_df$class <- "static"
atac_res_Shrink_df$class[atac_res_Shrink_df$peakID %in% gained$peakID] <- "gained"
atac_res_Shrink_df$class[atac_res_Shrink_df$peakID %in% lost$peakID] <- "lost"

# Save the significant peaks
save(diff_atac_sig, gained, lost, static,atac_res_Shrink_df,
     file = "atac_diff_significant_peaks.RData")

#--------------------------------------------------------------------------------

