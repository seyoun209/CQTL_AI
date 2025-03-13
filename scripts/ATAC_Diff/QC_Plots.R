# QC plots
#Libraries
library(corrplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(cqn)
library(rasqualTools)

macs2 <- read_and_clean_peaks("peaks/merged/allsamples_macs2_merged_counts.txt")
# Calculate number of samples and minimum count threshold
n_samples <- ncol(macs2) - 6  # Subtract 6 for the metadata columns (Geneid, Chr, Start, End, Strand, Length)
min_count <- 5 * n_samples

# Calculate row sums only for count columns (excluding metadata columns)
row_sums <- rowSums(macs2[, 7:ncol(macs2)])

# Filter low countds the data frame
macs2_filtered <- macs2[row_sums >= min_count, ]

macs2_gr <- makeGRangesFromDataFrame(macs2_filtered)
macs2_gr$peakID <- macs2_filtered$Geneid
counts_macs2 <- as.matrix(macs2_filtered[,7:ncol(macs2_filtered)])

# GC bias check-----------------------------------------------------------------
#1.  Calculate GC contents

windowsGR <- macs2_gr
peakwidths = width(windowsGR)
message('calculating GC content')
peakseqs = getSeq(Hsapiens, 
                  seqnames(windowsGR),
                  start(windowsGR),
                  end(windowsGR))
GCcontent_nolength = sapply(peakseqs, function(x) letterFrequency(x, "GC", as.prob=TRUE))
GCcontent = sapply(peakseqs, 
                   function(x) (length(gregexpr("[GC]",x)[[1]]))/nchar(x))
# Save GC content
save(GCcontent, meta_final,
     file = '/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/QC/QTL_GCContent.RData')


# Look at GC content distribution
hist(GCcontent_nolength, 
     breaks=50, 
     main="Distribution of GC Content in Peaks",
     xlab="GC Content",
     ylab="Frequency")

hist(GCcontent, 
     breaks=50, 
     main="Distribution of GC Content in Peaks",
     xlab="GC Content",
     ylab="Frequency")


# Check extreme cases
# Very high GC peaks (>70%)
high_gc <- sum(GCcontent_nolength > 0.7)
high_gc_percent <- (high_gc/length(GCcontent_nolength)) * 100

# Very low GC peaks (<30%)
low_gc <- sum(GCcontent_nolength < 0.3)
low_gc_percent <- (low_gc/length(GCcontent_nolength)) * 100

# Based on the value, most of my peaks (about 91%) have normal GC content (30-70%). I am planning to do regular vsd normalization. 
#-------------------------------------------------------------------------------
message('running cqn')


flagstat_data <- fread("/work/users/s/e/seyoun/CQTL_AI/output/wasp/blk_filter/qc_cqtl/multiqc_data/multiqc_samtools_flagstat.txt", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
 mapped_reads <- data.frame(
   File = gsub("_stats$", "", flagstat_data$Sample),  # Remove "_stats" suffix
   MappedReads = flagstat_data$mapped_passed          # Use mapped_passed column
 ) %>%
   rename(sampleID = File)
 
desired_order <- colnames(counts_macs2)
mapped_reads_ordered <- mapped_reads %>%
  filter(sampleID %in% desired_order) %>%
  slice(match(desired_order, sampleID))

cqn.counts.sqn = cqn(counts_macs2,
                     lengths = peakwidths,
                     x = GCcontent,
                     sizeFactors = mapped_reads_ordered$MappedReads,
                     sqn = TRUE)


# Calculate normalization factors
cqnOffset = cqn.counts.sqn$glm.offset
cqnNormFactors = exp(cqnOffset) #exponential function to the offsets to get the raw normalization factors.
cqnNormFactors = cqnNormFactors / exp(rowMeans(log(cqnNormFactors))) #normalization factors around 1.

# Save normalization factors
save(cqnNormFactors, mapped_reads_ordered, counts_macs2, macs2_gr,macs2_filtered,
     file = '/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/CQTL_cqnNormFactors.RData')



#Correlation plot---------------------------------------------------------------

#1. PCA from the ATAC-seq
vsd <- vst(chonATACdds, blind=TRUE)

mat_atac <- assay(vsd)
pca_comp <- prcomp(t(mat_atac))

# Create PCA data frame
pca_df <- as.data.frame(pca_comp$x)
variance_explained <- summary(pca_comp)$importance[2, ]

pca_df <- as.data.frame(pca_comp$x)
pca_vsd_cqtl <- merge(pca_df, colData(vsd), by="row.names")


# PBS

ctl_plot <- create_correlation_plot(pca_vsd_cqtl, condition = "CTL", title = "CTL Correlations")
fnf_plot <- create_correlation_plot(pca_vsd_cqtl, condition = "FNF", title = "FNF Correlations")


all_plot <- create_correlation_plot(pca_vsd_cqtl, condition = NULL)

