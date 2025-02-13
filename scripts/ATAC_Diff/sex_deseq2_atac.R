# Sex differential chromatin chondrocyte

setwd("/work/users/s/e/seyoun/CQTL_AI/output/")

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


meta_final$Sex <- factor(meta_final$Sex)
meta_final$AgeGroup <- factor(case_when(
  meta_final$Age <= 44 ~ "25_44",
  meta_final$Age <= 64 ~ "45_64",
  TRUE ~ "65_84"
))

#Change to macs2_gr
macs2_gr <- makeGRangesFromDataFrame(macs2)
macs2_gr$peakID <- macs2$Geneid
counts_macs2 <- as.matrix(macs2[,7:ncol(macs2)])


se <- SummarizedExperiment(assays=list(counts = counts_macs2),
                           rowRanges=macs2_gr, colData=meta_final)



# Run DESeq2 -------------------------------------------------------------------
# 1. Create DESeq2 object
sexATACdds <- DESeqDataSet(se,
                            design = ~ Condition + Sex )

# 2. Pre-filtering (remove low count regions)
n_samples <- ncol(sexATACdds)
min_count <- 5 * n_samples
keep <- rowSums(counts(sexATACdds)) >= min_count
sexATACdds <- sexATACdds[keep,]

# 3. Run DESeq2
# This includes:
# - estimateSizeFactors
# - estimateDispersions
# - nbinomWaldTest
sexATACdds <- DESeq(sexATACdds)

# 4. Get results
sex_Res <- results(sexATACdds)
print(summary(sex_Res,alpha=0.05, na.rm=TRUE))
sex_res_Shrink <- lfcShrink(sexATACdds, coef='Sex_M_vs_F', type="apeglm",format =  "GRanges")
sex_res_Shrink$peakID <- rowData(sexATACdds)$peakID

save(sexATACdds,sex_Res,sex_res_Shrink,
     file="/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/atac_chonSex_deseq.RData")
#-------------------------------------------------------------------------------
# Differential atac gained, lost ,  statics Chondrocyte
# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.05)

sex_res_Shrink <- sex_res_Shrink[!is.na(sex_res_Shrink$padj)]
sexNorm_counts <- counts(sexATACdds, normalized=TRUE)

# Get the condition information
Sex_cond <- colData(sexATACdds)$Sex

# Calculate the baseMean for each condition
baseMean_M <- rowMeans(sexNorm_counts[, Sex_cond == "M"])
baseMean_F <- rowMeans(sexNorm_counts[, Sex_cond == "F"])

sex_res_Shrink_df <- as.data.frame(sex_res_Shrink)

#Getting significant ATAC-signals
diff_chonSex_atac_sig <- sex_res_Shrink_df %>%
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)

#finiding up sig and down sig
gained <- diff_chonSex_atac_sig[diff_chonSex_atac_sig$log2FoldChange > 0 & diff_chonSex_atac_sig$padj < 0.05,]
lost <- diff_chonSex_atac_sig[diff_chonSex_atac_sig$log2FoldChange < 0 & diff_chonSex_atac_sig$padj < 0.05,]
static <- sex_res_Shrink_df[!(sex_res_Shrink_df$peakID %in% c(gained$peakID, lost$peakID)), ]


sex_res_Shrink_df$class <- "static"
sex_res_Shrink_df$class[sex_res_Shrink_df$peakID %in% gained$peakID] <- "gained"
sex_res_Shrink_df$class[sex_res_Shrink_df$peakID %in% lost$peakID] <- "lost"


write.table(gained, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/chondrocyte_FNF_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lost, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/chondrocyt_PBS_specific.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(chon_res_Shrink_df, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/chondrocyt_atac_allPeaks.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


#-------------------------------------------------------------------------------
#PCA

rldSex <- rlog(sexATACdds, blind=TRUE)
vsdSex <- vst(sexATACdds, blind=TRUE)
#plotPCA(rldSex, intgroup = "Sex",ntop=256743) + ggplot2::theme(aspect.ratio = 1)
plotPCA(vsdSex, intgroup = "Sex",ntop=256743) + ggplot2::theme(aspect.ratio = 1)

mat_atac <- assay(vsdSex)
pca_comp <- prcomp(t(mat_atac))

# Create PCA data frame
pca_df <- as.data.frame(pca_comp$x)
variance_explained <- summary(pca_comp)$importance[2, ]

pca_df <- as.data.frame(pca_comp$x)
pca_df$Donor <- colData(vsdSex)$Donor
pca_df$Condition <-colData(vsdSex)$Condition
pca_df$Sex <-colData(vsdSex)$Sex


ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Sex)) +
  geom_point(size = 3) +  # Increase point size here
  #geom_text_repel(size = 5, box.padding = unit(1, "lines")) +
  scale_color_manual(values = c("M" = "#9567FE", "F" = 	"#FD6320")) +
  labs(
    title = "Total peaks (n=256,743) with all samples",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))
