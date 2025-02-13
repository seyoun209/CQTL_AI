library(DESeq2)
library(ggplot2)
library(ggrepel)
library(GenomicFeatures)
library(ChIPseeker)

load('/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_cqnNormFactors.RData') # cqnNormFactors, results
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_CSAW_windowCounts.Rdata") # counts, my.regions, results
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_fragment_analysis.RData")

results$SampleID <- gsub("-", "_", results$SampleID)

results <- data.frame(
  Sample = sample_names,
  SampleID = sapply(strsplit(sample_names, "_"), `[`, 2),
  Condition = sapply(strsplit(sample_names, "_"), `[`, 3),
  Tissue = sapply(strsplit(sample_names, "_"), `[`, 4),
  #replicate = sapply(strsplit(sample_names, "_"), `[`, 5),
  MedianFragmentLength = frag_lens
)

my.regions$peakID <- paste0("peak_", seq_len(length(my.regions)))

# Create DESeq2 object
chonATACdds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(
    DonorID = factor(results$SampleID),
    Condition = factor(results$Condition, levels = c("CTL", "FNF"))
  ),
  rowRanges = my.regions,
  design = ~DonorID + Condition
)

# Apply CQN normalization factors
normalizationFactors(chonATACdds) <- cqnNormFactors

# Run DESeq2
chonATACdds <- DESeq(chonATACdds)


# Get results
chonRes <- results(chonATACdds)
print(summary(chonRes,alpha=0.05, na.rm=TRUE))
chon_res_Shrink <- lfcShrink(chonATACdds, coef='Condition_FNF_vs_CTL', type="apeglm",format =  "GRanges")
chon_res_Shrink$peakID <- rowData(chonATACdds)$peakID

save(chonATACdds,chonRes,chon_res_Shrink,file="/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/atac_chon_deseq.RData")

load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/atac_chon_deseq.RData")
#-------------------------------------------------------------------------------
# Differential atac gained, lost ,  statics Chondrocyte
# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.05)

chon_res_Shrink <- chon_res_Shrink[!is.na(chon_res_Shrink$padj)]

normalized_counts <- counts(chonATACdds, normalized=TRUE)

# Get the condition information
condition <- colData(chonATACdds)$Condition

# Calculate the baseMean for each condition
baseMean_CTL <- rowMeans(normalized_counts[, condition == "CTL"])
baseMean_FNF <- rowMeans(normalized_counts[, condition == "FNF"])

chon_res_Shrink_df <- as.data.frame(chon_res_Shrink)
#Getting significant ATAC-signals
diff_chon_atac_sig <- chon_res_Shrink_df %>%
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)

#finiding up sig and down sig
gained <- diff_chon_atac_sig[diff_chon_atac_sig$log2FoldChange > 0 & diff_chon_atac_sig$padj < 0.05,]
lost <- diff_chon_atac_sig[diff_chon_atac_sig$log2FoldChange < 0 & diff_chon_atac_sig$padj < 0.05,]
static <- diff_chon_atac_sig[!(diff_chon_atac_sig$peakID %in% c(gained$peakID, lost$peakID)), ]


chon_res_Shrink_df$class <- "static"
chon_res_Shrink_df$class[chon_res_Shrink_df$peakID %in% gained$peakID] <- "gained"
chon_res_Shrink_df$class[chon_res_Shrink_df$peakID %in% lost$peakID] <- "lost"


write.table(gained, file = "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/chondrocyte_FNF_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lost, file = "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/chondrocyt_PBS_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(chon_res_Shrink_df, file = "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/chondrocyt_atac_allPeaks.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)



#-------------------------------------------------------------------------------
#PCA plot

rld <- rlog(chonATACdds, blind=TRUE)
plotPCA(rld, intgroup = "Condition",ntop=241209) + ggplot2::theme(aspect.ratio = 1)

# For rlog transformation with PCA visualization:
mat_chon_atac <- assay(rld)
pca_comp_chon <- prcomp(t(mat_chon_atac))

# Create PCA data frame
pca_df <- as.data.frame(pca_comp_chon$x)
variance_explained <- summary(pca_comp_chon)$importance[2, ]

pca_df <- as.data.frame(pca_comp_chon$x)
pca_df$DonorID <- colData(chonATACdds)$DonorID
pca_df$Condition <- colData(chonATACdds)$Condition

chon_atac_PCAplot_all <- ggplot(pca_df, aes(x = PC1, y = PC2, label = DonorID, color = Condition)) +
  geom_point(size = 3) +  # Increase point size here
  #geom_text_repel(size = 5, box.padding = unit(1, "lines")) +
  scale_color_manual(values = c("CTL" = "#1775AE", "FNF" = "#FFC200")) +
  labs(
    title = "Total peaks with all samples-CQTL",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))+
  coord_fixed(ratio = 1.2) 

ggsave("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/plot/atac_pacPLOT.pdf", plot = last_plot(), width = 5, height = 5)
save(chon_atac_PCAplot_all, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/plot/atac_pacPLOT.rda")


