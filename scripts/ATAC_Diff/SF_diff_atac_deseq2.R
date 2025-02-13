library(DESeq2)
library(ggplot2)
library(ggrepel)
library(GenomicFeatures)
library(ChIPseeker)

load("/work/users/s/e/seyoun/CQTL_AI/output/old_peaks/csaw/synovium/synovium_CSAW_filter50_windowCounts.Rdata") # for counts and reg.counts_data
load('/work/users/s/e/seyoun/CQTL_AI/output/old_peaks/csaw/synovium/synovium_cqnNormFactors_filter50.RData') # for cqnNormFactors and results
load("/work/users/s/e/seyoun/CQTL_AI/output/old_peaks/csaw/synovium/synovium_regions_data.RData")
results$SampleID <- gsub("-", "_", results$SampleID)
data.regions$peakID <- paste0("peak_", seq_len(length(data.regions)))
se <- SummarizedExperiment(assays=list(counts = counts),
                           rowRanges=data.regions,  colData = data.frame(
                             DonorID = factor(results$SampleID),
                             Condition = factor(results$Condition, levels = c("CTL", "FNF"))
                           ))
# Create DESeq2 object
atacDDS <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(
    DonorID = factor(results$SampleID),
    Condition = factor(results$Condition, levels = c("CTL", "FNF"))
  ),
  design = ~DonorID + Condition,
  rowRanges = data.regions
)

# Apply CQN normalization factors
normalizationFactors(atacDDS) <- cqnNormFactors

# Run DESeq2
atacDDS <- DESeq(atacDDS)

# Get results
res <- results(atacDDS)
print(summary(res,alpha=0.05, na.rm=TRUE))
res_Shrink <- lfcShrink(atacDDS, coef='Condition_FNF_vs_CTL', type="apeglm",format =  "GRanges")
res_Shrink$peakID <- rowData(atacDDS)$peakID

save(atacDDS,res,res_Shrink,file="/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_deseq.Rdata")

#-------------------------------------------------------------------------------
# Differential atac gained, lost ,  statics
# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.05)

res_Shrink <- res_Shrink[!is.na(res_Shrink$padj)]

normalized_counts <- counts(atacDDS, normalized=TRUE)

# Get the condition information
condition <- colData(atacDDS)$Condition

# Calculate the baseMean for each condition
baseMean_CTL <- rowMeans(normalized_counts[, condition == "CTL"])
baseMean_FNF <- rowMeans(normalized_counts[, condition == "FNF"])

res_Shrink_df <- as.data.frame(res_Shrink)
#Getting significant ATAC-signals
diff_atac_sig <- res_Shrink_df %>%
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)

#finiding up sig and down sig
gained <- diff_atac_sig[diff_atac_sig$log2FoldChange > 0 & diff_atac_sig$padj < 0.05,]
lost <- diff_atac_sig[diff_atac_sig$log2FoldChange < 0 & diff_atac_sig$padj < 0.05,]
static <- diff_atac_sig[!(diff_atac_sig$peakID %in% c(gained$peakID, lost$peakID)), ]


res_Shrink_df$class <- "static"
res_Shrink_df$class[res_Shrink_df$peakID %in% gained$peakID] <- "gained"
res_Shrink_df$class[res_Shrink_df$peakID %in% lost$peakID] <- "lost"


write.table(gained, file = "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/SF_FNF_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lost, file = "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/SF_PBS_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(res_Shrink_df, file = "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/SF_atac_allPeaks.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)