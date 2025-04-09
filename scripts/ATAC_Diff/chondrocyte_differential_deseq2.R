#Feb 2nd 
#Author Seyoun Byun
# 02.04.2024
#This is updated to use the no-wasp
setwd("/work/users/s/e/seyoun/CQTL_AI/output/")
# Differential atac-seq -Chondrocyte
#-------------------------------------------------------------------------------
#library
library(DESeq2)
library(Rsubread,lib.loc = "/nas/longleaf/home/seyoun/R/x86_64-pc-linux-gnu-library/4.2")
library(GenomicRanges)
library(base)
library(data.table)
library(tidyverse)


#rocco <- read_and_clean_peaks("peaks/merged/allsamples_rocco_merged_counts.txt")
#macs3 <- read_and_clean_peaks("peaks/merged/allsamples_macs3_merged_counts.txt")
#hmmratac <- read_and_clean_peaks("peaks/merged/allsamples_hmmratac_merged_counts.txt")

macs2 <- read_and_clean_peaks("peaks/merged/allsamples_macs2_merged_counts.txt")

#Differential analysis 

## sample config file generate (Metasamplesheet)
config <- yaml::read_yaml("../config/ATACconfig.yaml")
samplesheet <- fread(paste0("/work/users/s/e/seyoun/CQTL_AI/",config$samplesheet))
donorsheet <- fread("/work/users/s/e/seyoun/CQTL_AI/Chon_DonorInfo.txt") 

#Check
donorsheet[grep("\\(ankle\\)|\\(femur\\)", donorsheet$Donor), ]
donorsheet <- donorsheet[!grepl("\\(femur\\)", Donor)]
donorsheet[, Donor := gsub(" \\(ankle\\)", "", Donor)]

#verify
donorsheet[grep("AM7778", donorsheet$Donor), ]

merged_meta <- left_join(samplesheet, donorsheet, by = "Donor") |> 
  filter(Tissue == "Ankle") |>
  dplyr::select("Proj","Donor","Condition","Tissue","Protocol_notes","Time","Tech_Rep","Seq_Rep","TreatmentDate",
                "ATACProtocolDate","ATACLibrarySubmissionDate","Sex","Age","Race",
                "OAGradeAvg","CauseOfDeath","FragmentBatch")



#Take the replicate out
pairs_count <- merged_meta[, .N, by=.(Donor, Condition)]
duplicated_pairs <- pairs_count[N > 1, .(Donor, Condition)]

meta_data_rep_cleaned <- merged_meta[!(Donor %in% unique(duplicated_pairs$Donor)) | 
                                       (Donor %in% unique(duplicated_pairs$Donor) & Tech_Rep == 1)]

#Add the age group
meta_data_rep_cleaned$AgeGroup <- cut(meta_data_rep_cleaned$Age,
                                      breaks=c(24, 44, 64, 84),
                                      labels=c("25-44", "45-64", "65-84"))



meta_final <- meta_data_rep_cleaned |> 
  mutate(sampleID = paste(Proj, Donor, Condition, Tissue, Protocol_notes, sep="_"))

meta_final$Donor <- factor(meta_final$Donor)
meta_final$Condition <- factor(meta_final$Condition, levels = c("CTL", "FNF"))

#Change to macs2_gr
macs2_gr <- makeGRangesFromDataFrame(macs2)
macs2_gr$peakID <- macs2$Geneid
counts_macs2 <- as.matrix(macs2[,7:ncol(macs2)])


se <- SummarizedExperiment(assays=list(counts = counts_macs2),
                               rowRanges=macs2_gr, colData=meta_final)

# Run DESeq2 -------------------------------------------------------------------
# 1. Create DESeq2 object
chonATACdds <- DESeqDataSet(se,
                            design = ~ Donor + Condition)

# 2. Pre-filtering (remove low count regions)
n_samples <- ncol(chonATACdds)
min_count <- 5 * n_samples
keep <- rowSums(counts(chonATACdds)) >= min_count
chonATACdds <- chonATACdds[keep,]

# 3. Run DESeq2
  # This includes:
  # - estimateSizeFactors
  # - estimateDispersions
  # - nbinomWaldTest
chonATACdds <- DESeq(chonATACdds)

# 4. Get results
chonRes <- results(chonATACdds)
print(summary(chonRes,alpha=0.05, na.rm=TRUE))
chon_res_Shrink <- lfcShrink(chonATACdds, coef='Condition_FNF_vs_CTL', type="apeglm",format =  "GRanges")
chon_res_Shrink$peakID <- rowData(chonATACdds)$peakID

save(chonATACdds,chonRes,chon_res_Shrink,
     file="/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/atac_chon_deseq.RData")
#-------------------------------------------------------------------------------
# Differential atac gained, lost ,  statics Chondrocyte
# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.05)

load("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/atac_chon_deseq.RData")

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
  dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.05)

#finiding up sig and down sig
gained <- diff_chon_atac_sig[diff_chon_atac_sig$log2FoldChange > 0 & diff_chon_atac_sig$padj < 0.05,]
lost <- diff_chon_atac_sig[diff_chon_atac_sig$log2FoldChange < 0 & diff_chon_atac_sig$padj < 0.05,]
static <- chon_res_Shrink_df[!(chon_res_Shrink_df$peakID %in% c(gained$peakID, lost$peakID)), ]


chon_res_Shrink_df$class <- "static"
chon_res_Shrink_df$class[chon_res_Shrink_df$peakID %in% gained$peakID] <- "gained"
chon_res_Shrink_df$class[chon_res_Shrink_df$peakID %in% lost$peakID] <- "lost"


write.table(gained, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/chondrocyte_FNF_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lost, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/chondrocyt_PBS_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(chon_res_Shrink_df, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/chondrocyt_atac_allPeaks.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


#-------------------------------------------------------------------------------
#PCA

rld <- rlog(chonATACdds, blind=TRUE)
vsd <- vst(chonATACdds, blind=TRUE)
plotPCA(rld, intgroup = c("Condition","AgeGroup"),ntop=400) + ggplot2::theme(aspect.ratio = 1)
plotPCA(vsd, intgroup = c("Condition","Sex"),ntop=256743) + ggplot2::theme(aspect.ratio = 1)

mat_atac <- assay(rld)
pca_comp <- prcomp(t(mat_atac))

# Create PCA data frame
pca_df <- as.data.frame(pca_comp$x)
variance_explained <- summary(pca_comp)$importance[2, ]

pca_df <- as.data.frame(pca_comp$x)
pca_df$Donor <- colData(rld)$Donor
pca_df$Condition <-colData(rld)$Condition

 ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Condition)) +
  geom_point(size = 3) +  # Increase point size here
  #geom_text_repel(size = 5, box.padding = unit(1, "lines")) +
  scale_color_manual(values = c("CTL" = "#1775AE", "FNF" = "#FFC200")) +
  labs(
    title = "Total peaks (n=256,743) with all samples",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))
 
 
# MA plot-----------------------------------------------------------------------


ggplot(chonRes, aes(x = log10(baseMean), y = log2FoldChange)) +
   geom_point(data = subset(chonRes, padj > 0.05), 
              color = "grey60", alpha = 0.5 )+
   geom_point(data = subset(chonRes, padj < 0.05 & log2FoldChange > 0), 
              color = "#E07653", alpha = 0.5) +
   geom_point(data = subset(chonRes, padj < 0.05 &  log2FoldChange < 0), 
              color = "#287C6F", alpha = 0.5) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
   theme_classic() +
   theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
         plot.background = element_rect(fill = 'transparent', color = "transparent"),
         text = element_text(family = "Helvetica"),
         legend.position = "None",
         axis.text.y = element_text(color = "black", size = 6),
         axis.title.y = element_text(size = 6),  # Changed from element_markdown
         axis.title.x = element_text(size = 6),  # Changed from element_markdown
         axis.text.x = element_text(color = "black", size = 6),
         strip.background = element_blank(),
         axis.ticks = element_blank(),
         axis.line.x = element_line(linewidth = 0.25),
         strip.text = element_text(size = 8, color = "black"),
         panel.spacing = unit(0, "mm"), 
         plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
   labs(x = expression("Mean of normalized counts (log"[10]*")"),
        y = expression("Log"[2]*" Fold change FN-f / PBS"))+
   coord_cartesian(ylim = c(-8,9))
 
 
#Heatmap------------------------------------------------------------------------
 

