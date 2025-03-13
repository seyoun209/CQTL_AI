#ATAC-seq for the gc fixed and use the cqn normailzied

#-------------------------------------------------------------------------------
library(DESeq2)

#loading data

#GC content
load('/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/QC/QTL_GCContent.RData')

#CQN Normalized factors (cqnNormFactors, mapped_reads_ordered, counts_macs2, macs2_gr,macs2_filtered)
load("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/CQTL_cqnNormFactors.RData")

#Prepare for DESeq2

atacdds <- DESeqDataSetFromMatrix(countData = counts_macs2,
                              colData = meta_final,
                              rowRanges=macs2_gr,
                              design = ~Donor+Condition)

# Apply CQN normalization factors
normalizationFactors(atacdds) <- cqnNormFactors

# Run DESeq2
atacdds <- DESeq(atacdds)

# Get results
atac_res <- results(atacdds)
print(summary(atac_res,alpha=0.05, na.rm=TRUE))
atac_res_Shrink <- lfcShrink(atacdds, coef='Condition_FNF_vs_CTL', type="apeglm",format =  "GRanges")
atac_res_Shrink$peakID <- rowData(atacdds)$peakID

save(atacdds,atac_res,atac_res_Shrink,
     file="/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/atac_cqnnorm_deseq2_re.RData")



#-------------------------------------------------------------------------------
#Differential analysis
load("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/atac_cqnnorm_deseq2_re.RData")

# Differential atac gained, lost ,  statics Chondrocyte
# Get significant differential ATAC peaks (abs log2FC > 1.5 and padj < 0.05)

atac_res_Shrink <- atac_res_Shrink[!is.na(atac_res_Shrink$padj)]

norm_ATACcounts <- counts(atacdds, normalized=TRUE)

# Get the condition information
condition <- colData(atacdds)$Condition

# Calculate the baseMean for each condition
baseMean_ATAC_pbs <- rowMeans(norm_ATACcounts[, condition == "CTL"])
baseMean_ATAC_fnf <- rowMeans(norm_ATACcounts[, condition == "FNF"])

atac_res_ShrinkDF <- as.data.frame(atac_res_Shrink)
#Getting significant ATAC-signals #24560
sig_atac_log2fc15_padj05 <- atac_res_ShrinkDF %>%
  dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.05)

#finiding up sig and down sig
gained <- sig_atac_log2fc15_padj05[sig_atac_log2fc15_padj05$log2FoldChange > 0 & sig_atac_log2fc15_padj05$padj < 0.05,]
lost <- sig_atac_log2fc15_padj05[sig_atac_log2fc15_padj05$log2FoldChange < 0 & sig_atac_log2fc15_padj05$padj < 0.05,]
static <- atac_res_ShrinkDF[!(atac_res_ShrinkDF$peakID %in% c(gained$peakID, lost$peakID)), ]


atac_res_ShrinkDF$class <- "static"
atac_res_ShrinkDF$class[atac_res_ShrinkDF$peakID %in% gained$peakID] <- "gained"
atac_res_ShrinkDF$class[atac_res_ShrinkDF$peakID %in% lost$peakID] <- "lost"


write.table(gained, 
            file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/deseq2_result/chondrocyte_FNF_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lost, 
            file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/deseq2_result/chondrocyt_PBS_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(atac_res_ShrinkDF, 
            file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/deseq2_result/chondrocyt_atac_allPeaks.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
#-------------------------------------------------------------------------------
#PCA plot