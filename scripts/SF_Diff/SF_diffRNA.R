# Synovial RNA_Seq and ATAC-seq for Doug's grant
#step1 for the Differential RNA-seq
setwd("/work/users/s/e/seyoun/CQTL_AI")
library(data.table)
library(tximeta)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(org.Hs.eg.db)
library(readr)
library(rrvgo)
library(gridtext)
library(ggtext)
library(colorspace)
source("scripts/utils/functions_rna.R")
#-------------------------------------------------------------------------------
#Make MA plot for RNA-seq synovial fibroblast


DonorInfo <- fread("/work/users/s/e/seyoun/CQTL_AI/SF_Donor_info.txt")
rna_sampleInfo <- fread("/work/users/s/e/seyoun/CQTL_AI/rna_samplesheet.txt")
rna_sampleInfo_sub <- rna_sampleInfo |> dplyr::select("Donor","Condition","RNAQubit","RIN")
DonorInfo_sub <- DonorInfo |> dplyr::select("Donor","Sex","Age")

dirnames <- dir("/work/users/s/e/seyoun/CQTL_AI/rna_output/quant",full.names=T)
sample_names <- basename(dirnames)

coldata.df <- data.frame(
  Sample = sample_names,
  SampleID = sapply(strsplit(sample_names, "_"), `[`, 2),
  Condition = sapply(strsplit(sample_names, "_"), `[`, 3),
  Tissue = sapply(strsplit(sample_names, "_"), `[`, 4)
)

coldata_SF <- coldata.df |> filter(Tissue == "synovium")
coldata_SF_fi <- coldata_SF %>%
  left_join(rna_sampleInfo_sub, by = c("SampleID" = "Donor", "Condition" ="Condition")) %>%
  left_join(DonorInfo_sub, by = c("SampleID" = "Donor"))

coldata_SF_fi$SampleID <- gsub("-", "_", coldata_SF_fi$SampleID)

## Add quant paths and names
coldata_SF_fi$files <- file.path("/work/users/s/e/seyoun/CQTL_AI/rna_output/quant", coldata_SF_fi$Sample, "quant.sf")
colnames(coldata_SF_fi) <- gsub("Sample", "names", colnames(coldata_SF_fi))
file.exists(coldata_SF_fi$files)

## Import data with tximeta & summarize to gene
se <- tximeta(coldata_SF_fi)
gse <- summarizeToGene(se)


#-------------------------------------------------------------------------------
#  Gene-level Synovial Fibroblast

## Convert to factors (avoids a warning)
colData(gse)[] <- lapply(colData(gse), factor)

## Build DESeq object
ddsRNA_SF <- DESeqDataSet(gse, design = ~ Condition  +namesID)

## Filter out lowly expressed genes
keep <- rowSums(counts(ddsRNA_SF) >= 10) >= ceiling(nrow(colData(gse))*0.5)
ddsRNA_SF <- ddsRNA_SF[keep,]

## Fit model
ddsRNA_SF <- DESeq(ddsRNA_SF)

# Get results
resRNA_SF <- results(ddsRNA_SF)
print(summary(resRNA_SF,alpha=0.05, na.rm=TRUE))
res_Shrink_RNA_SF_gr <- lfcShrink(ddsRNA_SF, coef='Condition_FNF_vs_CTL', type="apeglm", format =  "GRanges") |>
  plyranges::names_to_column("gene_id")

res_Shrink_RNA_SF <- lfcShrink(ddsRNA_SF, coef='Condition_FNF_vs_CTL', type="apeglm")

save(ddsRNA_SF,resRNA_SF,res_Shrink_RNA_SF,res_Shrink_RNA_SF_gr,
     file="/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/synovium_deseq.Rdata")

#------------------------------------------------------------------------------
#Differential gene expression -PCA plot
load("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/synovium_deseq.Rdata")

res_Shrini_RNA_SF <- res_Shrink_RNA_SF[!is.na(res_Shrink_RNA_SF$padj),]

norm_dds_SF_rna <- normTransform(ddsRNA_SF)
rlog_dds_SF_rna <- rlog(ddsRNA_SF)
vsd_dds_SF_rna <- vst(ddsRNA_SF)
plotPCA(vsd_dds_SF_rna,"namesID")



#assay(vsd_dds_SF_rna) <- limma::removeBatchEffect(assay(rlog_dds_SF_rna), rlog_dds_SF_rna$Sex)
#plotPCA(vsd_dds_SF_rna, "Condition",ntop=26446)

plotPCA(norm_dds_SF_rna, intgroup = "namesID",ntop=26446) + ggplot2::theme(aspect.ratio = 1)
#plotPCA(rlog_dds_SF_rna, intgroup ="Condition")

normalized_counts <- counts(ddsRNA_SF, normalized=TRUE)

# For rlog transformation with PCA visualization:
mat_chon_rna <- assay(rlog_dds_SF_rna)
pca_comp_chon <- prcomp(t(mat_chon_rna))

# Create PCA data frame
pca_df <- as.data.frame(pca_comp_chon$x)
variance_explained <- summary(pca_comp_chon)$importance[2, ]

pca_df <- as.data.frame(pca_comp_chon$x)
pca_df$DonorID <- colData(ddsRNA_SF)$namesID
pca_df$Condition <- colData(ddsRNA_SF)$Condition
pca_df$Sex <- colData(ddsRNA_SF)$Sex
pca_df$Age <- colData(ddsRNA_SF)$Age

ggplot(pca_df, aes(x = PC1, y = PC2, label = DonorID, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("CTL" = "#1775AE", "FNF" = "#FFC200")) +
  #scale_shape_manual(values = c("F" = 16, "M" = 17)) +  # Added this line
  labs(
    title = "Total peaks with all samples-CQTL",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent')) +
  coord_fixed(ratio = 1.2)

#-------------------------------------------------------------------------------
#MA PLOT

#DESeq2::plotMA(res_Shrink_RNA_SF)


#-------------------------------------------------------------------------------
# Differential atac gained, lost ,  statics
# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.05)
load("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/synovium_deseq.Rdata")
res_Shrink_RNA_SF_gr <- res_Shrink_RNA_SF_gr[!is.na(res_Shrink_RNA_SF_gr$padj)]
normalized_counts <- counts(ddsRNA_SF, normalized=TRUE)

# Get the condition information
condition <- colData(ddsRNA_SF)$Condition

# Calculate the baseMean for each condition
baseMean_CTL <- rowMeans(normalized_counts[, condition == "CTL"])
baseMean_FNF <- rowMeans(normalized_counts[, condition == "FNF"])
# Calculate baseMean
baseMean <- rowMeans(normalized_counts)
log2FC <- log2(baseMean_FNF/baseMean_CTL)

#Getting significant RNA-signals
diff_rna_sig <- res_Shrink_RNA_SF_gr %>%
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)


#finiding up sig and down sig
up <- diff_rna_sig[diff_rna_sig$log2FoldChange > 0 & diff_rna_sig$padj < 0.05,]
down <- diff_rna_sig[diff_rna_sig$log2FoldChange < 0 & diff_rna_sig$padj < 0.05,]
static <- diff_rna_sig[!(diff_rna_sig$gene_id %in% c(up$gene_id, down$gene_id)), ]

#adding the class
res_Shrink_RNA_SF_df <- as.data.frame(res_Shrink_RNA_SF_gr)
res_Shrink_RNA_SF_df$class <- "static"
res_Shrink_RNA_SF_df$class[res_Shrink_RNA_SF_df$gene_id %in% up$gene_id] <- "gained"
res_Shrink_RNA_SF_df$class[res_Shrink_RNA_SF_df$gene_id %in% down$gene_id] <- "lost"


#Adding the HGNC

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
res_Shrink_RNA_SF_df_noVer <-res_Shrink_RNA_SF_df |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))

annots <- AnnotationDbi::select(org.Hs.eg.db, res_Shrink_RNA_SF_df_noVer$gene_id,
                                columns=c("SYMBOL"), keytype="ENSEMBL")
annots <- dplyr::rename(annots, gene_id = ENSEMBL)
annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
res_Shrink_RNA_SF_df_hgnc <- left_join(res_Shrink_RNA_SF_df_noVer, annots_unique, by = "gene_id")


# Save R data
save(res_Shrink_RNA_SF_df,res_Shrink_RNA_SF_df_hgnc, diff_rna_sig, results,up, down, static, 
     file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/SF_differential_Deseq2_results.RData")
load("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/SF_differential_Deseq2_results.RData")
#------------------------------------------------------------------------------
# plotting for the MA 

ggplot(res_Shrink_RNA_SF_df_hgnc, aes(x = log10(baseMean), y = log2FoldChange)) +
  # Plot grey points first
  geom_point(data = subset(res_Shrink_RNA_SF_df, padj > 0.05), 
             color = "gray", alpha = 0.5) +
  # Plot colored points on top
  geom_point(data = subset(res_Shrink_RNA_SF_df, padj <= 0.05 & class == "gained"), 
             color = "#e07653", alpha = 0.5) +
  geom_point(data = subset(res_Shrink_RNA_SF_df, padj <= 0.05 & class == "lost"), 
             color = "#1e87a5", alpha = 0.5) +
  theme_classic() +
  ggtitle("RNA") +
  labs(x = "Mean of normalized counts (log10)",
       y = "Log2 fold change FNF vs CTL") +
  coord_cartesian(ylim = c(-10,11)) 
  #scale_x_log10(breaks=c(1,2,5,10,20,50,100,500,2000,12000,50000), 
  #              limits = c(10,150000))


up_significantRNA_SF_df <- res_Shrink_RNA_SF_df_hgnc %>%
  filter(class == "gained", 
         log10(baseMean) > 4,
         log2FoldChange > 4)

down_significantRNA_SF_df <- res_Shrink_RNA_SF_df_hgnc %>%
  filter(class == "lost", 
         log10(baseMean) > 2.7,
         log2FoldChange < - 3.25 )



up_significantRNA_SF_df <- res_Shrink_RNA_SF_df_hgnc %>%
  filter(class == "gained", 
         log2FoldChange > 0) |> dplyr::select( gene_id, SYMBOL, padj, log2FoldChange)

write.table(up_significantRNA_SF_df, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_up_SFRNA.txt",
            sep='\t',quote=F,row.names=F,col.names=T)

down_significantRNA_SF_df <- res_Shrink_RNA_SF_df_hgnc %>%
  filter(class == "lost", log2FoldChange < 0
 ) |> dplyr::select( gene_id, SYMBOL, padj, log2FoldChange)
write.table(down_significantRNA_SF_df, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_down_SFRNA.txt",
            sep='\t',quote=F,row.names=F,col.names=T)

RNA_MAplot <- ggplot(res_Shrink_RNA_SF_df_hgnc, aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(data = subset(res_Shrink_RNA_SF_df, padj > 0.05), 
             color = "gray", alpha = 0.5) +
  geom_point(data = subset(res_Shrink_RNA_SF_df, padj <= 0.05 & class == "gained"), 
             color = "#e07653", alpha = 0.5) +
  geom_point(data = subset(res_Shrink_RNA_SF_df, padj <= 0.05 & class == "lost"), 
             color = "#1e87a5", alpha = 0.5) +
  # Points with black outline
  geom_point(data = up_significantRNA_SF_df,
             color = "black", fill = "#e07653", size = 2, shape = 21, stroke = 1) +
  geom_point(data = down_significantRNA_SF_df,
             color = "black", fill = "#1e87a5", size = 2, shape = 21, stroke = 1) +
  # Text labels with black connecting lines
  ggrepel::geom_text_repel(data = up_significantRNA_SF_df,
                           aes(label = SYMBOL),
                           color = "black",
                           size = 2,
                           segment.color = "black",
                           segment.size = .15,
                           box.padding = 0.75,
                           point.padding = 0.75,
                           max.overlaps = Inf)+
  ggrepel::geom_text_repel(data = down_significantRNA_SF_df,
                           aes(label = SYMBOL),
                           color = "black",  # Changed text color to black
                           size = 2,
                           segment.color = "black",
                           segment.size = .15,
                           box.padding = 0.75,
                           point.padding = 0.75,
                           max.overlaps = Inf)+
  theme_classic() +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_markdown(size = 6),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_text(size = 8, color = "black"),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8))+
  labs(x = "Mean of normalized counts (log10)",
       y = "Log2 Fold change FNF vs. CTL") +
  coord_cartesian(ylim = c(-10,11))

output_dir <- "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/plot"
dir.create(output_dir, showWarnings =TRUE)

ggsave("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/plot/rnaMA_withpadding.pdf", plot = last_plot(), width = 5, height = 5)
save(RNA_MAplot, file = "rna_output/Differential_analysis/plot/RNA_MAplot.rda")

#-------------------------------------------------------------------------------
# GO and KEGG pathway from homer
res_Shrink_RNA_SF_df_noVer <-res_Shrink_RNA_SF_df |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))
up_RNAgeneslfc1padj05 <- diff_rna_sig[diff_rna_sig$log2FoldChange > 0 & diff_rna_sig$padj < 0.05,] |> 
  dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
  as.data.frame() |> 
  dplyr::select(gene_id) |> distinct() 
  

down_RNAgeneslfc1padj05 <- diff_rna_sig[diff_rna_sig$log2FoldChange < 0 & diff_rna_sig$padj < 0.05,] |> 
  dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
  as.data.frame() |> 
  dplyr::select(gene_id) |> distinct() 

all_RNAgeneslfc1padj05 <- diff_rna_sig |> 
  dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
  as.data.frame() |> 
  dplyr::select(gene_id) |> distinct() 


background_genes <- res_Shrink_RNA_SF_df |> 
  dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
  as.data.frame() |> 
  dplyr::select(gene_id) |> distinct() 


write.table(up_RNAgeneslfc1padj05, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_up_RNA_lfc1_padj05.txt",
            sep='\t',quote=F,row.names=F,col.names=F)
write.table(down_RNAgeneslfc1padj05, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_down_RNA_lfc1_padj05.txt",
            sep='\t',quote=F,row.names=F,col.names=F)
write.table(all_RNAgeneslfc1padj05, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_all_RNA_lfc1_padj05.txt",
            sep='\t',quote=F,row.names=F,col.names=F)

write.table(background_genes, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/backgroundRNA_genes.txt",
            sep='\t',quote=F,row.names=F,col.names=F)

#background gene set

#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_homer.sh /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_all_RNA_lfc1_padj05.txt /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/backgroundRNA_genes.txt /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/homer/homer_sig_RNA_all_LFC1_padj05")
#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_homer.sh /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_up_RNA_lfc1_padj05.txt /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/backgroundRNA_genes.txt /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/homer/homer_sig_RNA_up_LFC1_padj05")
#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_homer.sh /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/sig_down_RNA_lfc1_padj05.txt /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/backgroundRNA_genes.txt /work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/homer/homer_sig_RNA_down_LFC1_padj05")


#-------------------------------------------------------------------------------
# Get reduced, significant GO terms for each category
#GO

upsig_go_data <- read_delim("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/homer/homer_sig_RNA_up_LFC1_padj05/biological_process.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01)
upsig_go <- reduceGO(go_data,
                   category = "Upregulated")

## Format and write to table
upgo_table <- upsig_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(upgo_table, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/tables/upsigRNA_GO.csv")


downsig_go_data <- 
  read_delim("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/homer/homer_sig_RNA_down_LFC1_padj05/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
downsig_go <- reduceGO(downsig_go_data,
                       category = "Downregulated")

## Format and write to table
downgo_table <- downsig_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(downgo_table, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/tables/downsigRNA_GO.csv")

# Select 5 each for plotting
upsig_go_plotting <- upsig_go |> 
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("immune response", 
                           "defense response", 
                           "response to cytokine", 
                           "cell surface receptor signaling pathway", 
                           "regulation of leukocyte activation")) |> 
  arrange(`-log10pval`)

downsig_go_plotting <- downsig_go |> 
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("cell adhesion", 
                           "actin filament-based process",
                           "tissue development", 
                           "sodium ion transmembrane transport",
                           "bone remodeling")) |> 
  arrange(`-log10pval`)

# Combine into one
go_plotting <- bind_rows(upsig_go_plotting, downsig_go_plotting)
go_plotting$parentTerm <- factor(go_plotting$parentTerm, levels = go_plotting$parentTerm) 
go_plotting$category <- factor(go_plotting$category, levels = c("Upregulated", "Downregulated"))


# Plot all in barplot
GO_barplots <- ggplot(go_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 30, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 40, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 50, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 51), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 50, 10)) +
  scale_fill_manual(values = c(colorcode_RNA[["up"]], colorcode_RNA[["down"]])) +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = parentTerm), hjust = 0, family = "Helvetica",
            size = 2) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_text(size = 8, color = "black"),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("GO Terms")

#ggsave(filename = "plots/conditionDE_Fig1/GO_barplots.pdf",
#       plot = GO_barplots, width = 5, height = 8, units = "in")
#save(GO_barplots, file = "plots/conditionDE_Fig1/GO_barplots.rda")

#All -GO 

allsig_go_data <- 
  read_delim("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/homer/homer_sig_RNA_all_LFC1_padj05/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
allsig_go <- reduceGO(allsig_go_data,
                       category = "All")

## Format and write to table
allgo_table <- allsig_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(allgo_table, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/tables/allsigRNA_GO.csv")

# Select 5 each for plotting
sig_go_plotting <- allgo_table |> 
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("response to external stimulus", 
                           "immune response", 
                           "response to cytokine",
                           "signaling", 
                           "anatomical structure development",
                           "regulation of developmental process")) |> 
  arrange(`-log10pval`)

sig_go_plotting$parentTerm <- factor(sig_go_plotting$parentTerm, levels = sig_go_plotting$parentTerm) 
sig_go_plotting$category <- c("GO Terms")


GO_term_plot <- ggplot(sig_go_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 10, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 30, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 40, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 50 , color = "grey95", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 50, 10)) +
  scale_fill_manual(values = "grey80") +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = parentTerm), hjust = 0, family = "Helvetica",
            size = 2) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_markdown(size = 6),
        panel.spacing = unit(0, "mm"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) 
ggsave(filename = "rna_output/Differential_analysis/plot/allgenes_RNA_GO_barplots.pdf",
       plot = GO_term_plot, width = 5, height = 8, units = "in")
save(GO_term_plot, file = "rna_output/Differential_analysis/plot/allgenes_RNA_GO_barplots.rda")

#-------------------------------------------------------------------------------
#kEGG pathway


# Read in from Homer
allsig_kegg_data <- read_delim("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/homer/homer_sig_RNA_all_LFC1_padj05/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval)) |> 
  mutate(category = "all")

## Format and write to table
allsig_kegg_table <- allsig_kegg_data |> 
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(allsig_kegg_table, file = "/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/tables/KEGG_RNA_all.csv")


# Plot top 5 significant for each category
allsig_kegg_plotting <- allsig_kegg_data |> 
  filter(Term %in% c("Cytokine-cytokine receptor interaction",
                     "TNF signaling pathway", 
                     "IL-17 signaling pathway", 
                     "Jak-STAT signaling pathway", 
                     "Rheumatoid arthritis",
                     "NOD-like receptor signaling pathway")) |> 
  arrange(`-log10pval`)

allsig_kegg_plotting$Term <- factor(allsig_kegg_plotting$Term, levels = allsig_kegg_plotting$Term) 
allsig_kegg_plotting$category <- factor("KEGG Pathways")


KEGG_barplots <- ggplot(allsig_kegg_plotting, aes(x = `-log10pval`, y = Term, fill = category)) +
  geom_vline(xintercept = 5, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 10, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 15, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 25, color = "grey95", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~pval", limits = c(0, 26),
                     breaks = seq(0, 25, 5)) +
  scale_fill_manual(values ="grey80") +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = Term), hjust = 0, family = "Helvetica",
            size = 2) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_markdown(size = 6),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) 

ggsave(filename = "rna_output/Differential_analysis/plot/allgenes_RNA_KEGG_barplots.pdf",
       plot = KEGG_barplots, width = 5, height = 8, units = "in")
save(KEGG_barplots, file = "rna_output/Differential_analysis/plot/allgenes_RNA_KEGG_barplots.rda")




# Making a figure for Doug's grant
#-------------------------------------------------------------------------------
# Signal data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_AI/rna_output/signals/synovium/merged_signal/"
ctl_signal <- paste0(dir_merged,"CTL_sorted.bw")
fnf_signal <- paste0(dir_merged,"FNF_sorted.bw")

pdf(paste0("rna_output/Differential_analysis/plot/RNAplot.pdf"),width=11,height=3.5)
#RNA----------------------------------------------------------------------------
pageCreate(width = 11, height =3.5 , default.units = "inches", showGuides = FALSE)

#A-MA-plot---------------------------------------------------------------------
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load("rna_output/Differential_analysis/plot/RNA_MAplot.rda")
plotGG(plot = RNA_MAplot, x = 0.4, y = 0.5, height = 3.1, width = 3.5)
#GO and kegg-plot---------------------------------------------------------------

plotText("B", x = 4.0, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("rna_output/Differential_analysis/plot/allgenes_RNA_GO_barplots.rda")
load("rna_output/Differential_analysis/plot/allgenes_RNA_KEGG_barplots.rda")

plotGG(plot = KEGG_barplots, x = 4.3, y = 0.4, height = 1.7, width = 3)
plotGG(plot = GO_term_plot, x = 4.3, y = 1.95, height = 1.7, width = 3)

#RNA_GWAS_signal----------------------------------------------------------------

plotText("C", x = 7.5, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

GWAS_manhattan_plot(
  gwas_data = gwas_data_subset,
  gene_highlights = gene_hl,
  label= paste0("GWAS-", "AllOA"),
  rsID = gwas_rsID,
  ctl_signal = ctl_signal,
  fnf_signal = fnf_signal,
  x_start = 8.3,
  y_start = 0.6,
  width = 2.5,
  height = 1,
  signal_height = 0.43,
  padding = 0
)

dev.off()

