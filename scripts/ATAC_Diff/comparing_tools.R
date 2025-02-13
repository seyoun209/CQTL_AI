# Compare the result between MACS3 vs CSAW -Synovial Fibroblast

setwd("/work/users/s/e/seyoun/CQTL_AI/output/peaks/Ankle/counts/")

load('/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_cqnNormFactors.RData') # cqnNormFactors, results
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_CSAW_windowCounts.Rdata") # counts, my.regions, results
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_fragment_analysis.RData")

#load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw//synovium_regions.RData")
#load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_fragment_analysis.RData")
#load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_CSAW_windowCounts.Rdata")

# Convert CSAW regions to data frame format similar to MACS/HMMRATAC
csaw_df <- data.frame(
  Geneid = paste0("peak_", 1:length(my.regions)),
  Chr = seqnames(my.regions),
  Start = start(my.regions),
  End = end(my.regions),
  Strand = ".",
  Length = width(my.regions)
)

# Add count data 
count_cols <- as.data.frame(counts)
names(count_cols) <- results$Sample
csaw_df <- cbind(csaw_df, count_cols)

# Read and clean data
read_and_clean_peaks <- function(file) {
  df <- read.table(file, header=TRUE, skip=1, sep="\t")
  colnames(df) <- gsub("X.work.users.s.e.seyoun.CQTL_AI.output.align.", "", colnames(df))
  colnames(df) <- gsub("_RG.bam", "", colnames(df))
  return(df)
}

macs3 <- read_and_clean_peaks("macs3_counts.txt")
macs2 <- read_and_clean_peaks("macs2_counts.txt") 
hmmratac <- read_and_clean_peaks("hmmratac_counts.txt")

counts_macs3 <- as.matrix(macs3[,7:ncol(macs3)])
counts_macs2 <- as.matrix(macs2[,7:ncol(macs2)])
counts_hmmratac <- as.matrix(hmmratac[,7:ncol(hmmratac)])

macs3_peaks <- makeGRangesFromDataFrame(macs3)
macs2_peaks <- makeGRangesFromDataFrame(macs2)
hmmratac_peaks <- makeGRangesFromDataFrame(hmmratac)

# Calculate SNR for each peak caller signal to noise ratio
calculate_snr <- function(counts_df) {
  signal <- rowMeans(counts_df[,7:ncol(counts_df)])
  background <- apply(counts_df[,7:ncol(counts_df)], 1, sd)
  mean(signal/background, na.rm=TRUE)
}

# peak characteristics
compare_peaks <- function(macs3_df, macs2_df, hmmratac_df) {
  data.frame(
    Method = c("MACS3", "MACS2", "HMMRATA"),
    Peaks = c(nrow(macs3_df), nrow(macs2_df), nrow(hmmratac_df)),
    MedianWidth = c(
      median(macs3_df$End - macs3_df$Start),
      median(macs2_df$End - macs2_df$Start),
      median(hmmratac_df$End - hmmratac_df$Start)
      #median(as.data.frame(my.regions)$end -as.data.frame(my.regions)$start)
    ),
    SNR = c(
      calculate_snr(macs3_df),
      calculate_snr(macs2_df),
      calculate_snr(hmmratac_df)
    )
  )
}


# Calculate metrics
snr_results_macs3 <- calculate_snr(macs3)
snr_results_macs2 <- calculate_snr(macs2)
snr_results_hmmratac <- calculate_snr(hmmratac)
snr_results_cs <- calculate_snr(counts)
peak_comparison <- compare_peaks(macs3, macs2, hmmratac)


#-------------------------------------------------------------------------------
process_peaks <- function(hmmratac_df, csaw_df) {
  # Calculate read coverage and width distribution
  peak_metrics <- data.frame(
    Tool = c("HMMRATAC", "CSAW"),
    MedianWidth = c(
      median(hmmratac_df$End - hmmratac_df$Start),
      median(csaw_df$End - csaw_df$Start)
    ),
    SNR = c(
      calculate_snr(hmmratac_df),
      calculate_snr(csaw_df)
    ),
    TotalPeaks = c(nrow(hmmratac_df), nrow(csaw_df))
  )
  return(peak_metrics)
}

comparison <- process_peaks(hmmratac, csaw_df)
tss <- promoters(genes(txdb), upstream=1500, downstream=500)


csaw_gr <- makeGRangesFromDataFrame(csaw_df)
hmmratac_gr <- makeGRangesFromDataFrame(hmmratac)
macs3_gr <- makeGRangesFromDataFrame(macs3)

# Count overlaps with TSS regions
csaw_tss <- countOverlaps(tss, csaw_gr)
hmmratac_tss <- countOverlaps(tss, hmmratac_gr)
macs3_tss <- countOverlaps(tss, macs3_gr)

# Calculate enrichment scores
csaw_enrich <- mean(csaw_tss > 0)  # Fraction of TSS with peaks
hmmratac_enrich <- mean(hmmratac_tss > 0)
macs3_enrich <- mean(macs3_tss > 0)

print(paste("CSAW TSS enrichment:", round(csaw_enrich, 3)))
print(paste("HMMRATAC TSS enrichment:", round(hmmratac_enrich, 3)))

csaw_norm <- csaw_enrich/nrow(csaw_df)
hmmratac_norm <- hmmratac_enrich/nrow(hmmratac)
macs3_norm <- hmmratac_enrich/nrow(macs3)

# Calculate peak density in TSS regions
csaw_density <- sum(csaw_tss)/length(tss)
hmmratac_density <- sum(hmmratac_tss)/length(tss)
macs3_density <- sum(macs3_tss)/length(tss)

print(paste("CSAW normalized enrichment:", format(csaw_norm, scientific=TRUE)))
print(paste("HMMRATAC normalized enrichment:", format(hmmratac_norm, scientific=TRUE)))
print(paste("CSAW TSS peak density:", round(csaw_density, 2)))
print(paste("HMMRATAC TSS peak density:", round(hmmratac_density, 2)))
print(paste("MACSe TSS peak density:", round(macs3_density, 2)))

#-------------------------------------------------------------------------------
library(DiffBind)
library(ChIPQC)


txdb <- loadDb("/work/users/s/e/seyoun/crispr/02.test_seq/gencode.v45.annotation.TxDb")

# Annotate peaks with genomic features
csaw_anno <- annotatePeak(my.regions, TxDb=txdb)
hmmratac_anno <- annotatePeak(hmmratac_peaks, TxDb=txdb)
macs3_anno <- annotatePeak(macs3_peaks, TxDb=txdb)
macs2_anno <- annotatePeak(macs2_peaks, TxDb=txdb)


# Compare genomic distributions
plotAnnoBar(list(CSAW=csaw_anno, HMMRATAC=hmmratac_anno, MACS3=macs3_anno, MACS2=macs2_anno))

# Get numeric annotation stats
csaw_stats <- as.data.frame(csaw_anno@annoStat)
hmmratac_stats <- as.data.frame(hmmratac_anno@annoStat)
macs3_stats <- as.data.frame(macs3_anno@annoStat)




library(pheatmap)
cor_csaw <- cor(counts)
cor_hmm <- cor(counts_hmmratac)
cor_hmm <- cor(counts_hmmratac)
pheatmap(cor_csaw, main="CSAW replicate correlation")
pheatmap(cor_hmm, main="HMMRATAC replicate correlation")

ggplot(data.frame(width=width(csaw_peaks), 
                  strength=log2(rowMeans(counts))), 
       aes(x=width, y=strength)) +
  geom_hex() +
  labs(title="CSAW Peak Width vs Strength")

ggplot(data.frame(width=width(hmmratac_peaks), 
                  strength=log2(rowMeans(counts_hmmratac))), 
       aes(x=width, y=strength)) +
  geom_hex() +
  labs(title="MACS3 Peak Width vs Strength")


ggplot(data.frame(width=width(macs3_peaks), 
                  strength=log2(rowMeans(counts_macs3))), 
       aes(x=width, y=strength)) +
  geom_hex() +
  labs(title="MACS3 Peak Width vs Strength")


ggplot() +
  geom_density(data=data.frame(width=width(csaw_peaks)), aes(x=width, color="CSAW")) +
  geom_density(data=data.frame(width=width(hmmratac_peaks)), aes(x=width, color="MACS3")) +
  xlim(0, 1000) +
  labs(title="Peak Width Distribution")

hist(width(macs3_peaks), 
     breaks=50, 
     main="Peak Width Distribution", 
     xlab="Peak Width (bp)")

hist(width(csaw_peaks), 
     breaks=50, 
     main="Peak Width Distribution", 
     xlab="Peak Width (bp)")

hist(width(hmmratac_peaks), 
     breaks=50, 
     main="Peak Width Distribution", 
     xlab="Peak Width (bp)")
