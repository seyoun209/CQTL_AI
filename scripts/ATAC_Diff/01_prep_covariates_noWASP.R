# differential atacseq without wasp
# edited by: Jess
# date: 2023-03-30
#---------------------------------------------------------------
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
library(corrplot)

setwd("/work/users/s/e/seyoun/CQTL_AI/output/")
#---------------------------------------------------------------
# 1. GC bias check
load("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/Chon_macs2_se.RData") # load the data
# Calculate peak widths
peakwidths = width(macs2_gr)

# Calculate GC content
message('Calculating GC content...')
peakseqs = getSeq(BSgenome.Hsapiens.UCSC.hg38,
                  seqnames(macs2_gr),
                  start(macs2_gr),
                  end(macs2_gr))
GCcontent = letterFrequency(peakseqs, "GC", as.prob=TRUE) |> as.numeric()

hist(GCcontent,
     breaks=50,
     main="Distribution of GC Content in Peaks",
     xlab="GC Content",
     ylab="Frequency")

save(GCcontent, 
     file = '/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/GCContent.RData')


# Get library size information
message('Reading library size metrics...')
#Shell to run
#ml multiqc/1.25.1
# multiqc * -o /work/users/s/e/seyoun/CQTL_AI/output/filtered/blk_filter/qc_cqtl
flagstat_data <- fread("/work/users/s/e/seyoun/CQTL_AI/output/filtered/blk_filter/qc_cqtl/multiqc_data/multiqc_samtools_flagstat.txt",
                      header = TRUE,
                      stringsAsFactors = FALSE)

mapped_reads <- data.frame(
  sampleID = gsub("_stats$", "", flagstat_data$Sample),
  MappedReads = flagstat_data$mapped_passed
) %>%
  filter(sampleID %in% colnames(counts_macs2)) %>%
  dplyr::slice(match(colnames(counts_macs2), sampleID))

  # Apply CQN normalization
message('Performing CQN normalization...')
cqn.counts <- cqn(counts_macs2,
                  lengths = peakwidths,
                  x = GCcontent,
                  sizeFactors = mapped_reads$MappedReads,
                  sqn = TRUE)

# Calculate normalization factors for RASQUAL
message('Calculating RASQUAL normalization factors...')
cqnOffset <- cqn.counts$offset
cqnNormFactors <- exp(-cqnOffset)  # Convert offsets to factors (note negative sign)

# Center normalization factors around 1.0 (optional but recommended)
cqnNormFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))

# Create normalized counts (adding this missing piece)
message('Calculating normalized expression values...')
normalized_counts <- cqn.counts$y + cqn.counts$offset


# Save the CQN results for inspection/QC
save(cqn.counts, GCcontent, peakwidths,cqnNormFactors,
     file = "/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/CQN_results.RData")

par(mfrow=c(1,2))
# Before correction
smoothScatter(GCcontent, rowMeans(log2(counts_macs2 + 1)),
              main="Before GC Correction",
              xlab="GC Content",
              ylab="Mean log2(Counts + 1)")
abline(lm(rowMeans(log2(counts_macs2 + 1)) ~ GCcontent), col="red", lwd=2)
text(0.2, max(rowMeans(log2(counts_macs2 + 1)))*0.9, 
     paste("Slope =", round(pre_slope, 2)), col="red")

# After correction
smoothScatter(GCcontent, rowMeans(normalized_counts),
              main="After CQN Correction",
              xlab="GC Content",
              ylab="Mean Normalized Expression")
abline(lm(rowMeans(normalized_counts) ~ GCcontent), col="red", lwd=2)
text(0.2, max(rowMeans(normalized_counts))*0.9, 
     paste("Slope =", round(post_slope, 2)), col="red")
#---------------------------------------------------------------
#  Building deseq2 matrix

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
     file="/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/atac_cqnnorm_deseq2_re.RData")


#----------------------------------------------------------------
# 3. PCA plot

message("Calculating PCA for RASQUAL covariates...")
# Calculate variance-stabilized transformation
normCounts <- counts(atacdds, normalized = TRUE)
vsd <- vst(atacdds, blind=TRUE)
rld <- rlog(atacdds, blind=TRUE)

save(vsd, rld,normCounts,
 file = "/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/Chon_macs2_normalized.RData")
# Calculate PCA on vst-transformed data
pca_result <- prcomp(t(assay(vsd)))

# Create PCA data frame
pca_df <- as.data.frame(pca_result$x)
variance_explained <- summary(pca_result)$importance[2, ]
pca_df$Donor <- colData(vsd)$Donor
pca_df$Condition <-colData(vsd)$Condition


# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_var <- cumsum(var_explained)

# Plot variance explained
pdf("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/PCA_variance_explained.pdf", width=10, height=5)
par(mfrow=c(1,2))
barplot(var_explained[1:20]*100, 
        main="Variance Explained by PCs",
        xlab="Principal Component", 
        ylab="Percent Variance Explained",
        names.arg=1:20)
plot(cumulative_var[1:20]*100, 
     type="b", pch=19,
     xlab="Principal Component",
     ylab="Cumulative Percent Variance Explained",
     main="Cumulative Variance")
abline(h=75, col="red", lty=2)
dev.off()

#plotPCA(vsd, intgroup = c("ATACLibrarySubmissionDate"),ntop=281227) + 
ggplot2::theme(aspect.ratio = 1)

pca_plot_vsd <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Condition)) +
  geom_point(size = 3) +  # Increase point size here
  #geom_text_repel(size = 5, box.padding = unit(1, "lines")) +
  scale_color_manual(values = c("CTL" = "#1775AE", "FNF" = "#FFC200")) +
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
save(pca_plot_vsd,
file = "/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/plots/PCA_plot_vsd.rds")

load("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/plots/PCA_plot_vsd.rds")
print(pca_plot_vsd)
#----------------------------------------------------------------

