# QC plots
#Libraries
library(corrplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(cqn)


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

