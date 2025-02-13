library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(GenomicFeatures)
library(ChIPseeker)
library(ggseqlogo)
library(rsvg)
library(grid)
library(rrvgo)
library(ggnewscale)
library(patchwork)
library(rsvg)
library(grImport2)
library(ggrepel)
library(ComplexHeatmap)
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_CSAW_filter50_windowCounts.Rdata") # for counts and reg.counts_data
load('/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_cqnNormFactors_filter50.RData') # for cqnNormFactors and results
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_regions_data.RData")

#coldata_SF_fi <- results %>%
#  left_join(rna_sampleInfo_sub, by = c("SampleID" = "Donor", "Condition" ="Condition")) %>%
#  left_join(DonorInfo_sub, by = c("SampleID" = "Donor"))

results$SampleID <- gsub("-", "_", results$SampleID)
data.regions$peakID <- paste0("peak_", seq_len(length(data.regions)))

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

load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_deseq.Rdata")


#-------------------------------------------------------------------------------
#Annotation of the gene
txdb <- loadDb("/work/users/s/e/seyoun/crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)
promoters <- promoters(genes, upstream =1500, downstream =500)
# # Join results with gene info
# de_genes_shrink <-
#   inner_join(x = as.data.frame(res_Shrink),
#              y = as.data.frame(rowData(gse)) %>%
#                dplyr::select(c("gene_id", "tx_ids")),
#              by = "gene_id") %>%
#   makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
#   keepStandardChromosomes(pruning.mode = "coarse") %>%
#   as.data.frame()

rld <- rlog(atacDDS, blind=TRUE)
plotPCA(rld, intgroup = "Condition",ntop=214378) + ggplot2::theme(aspect.ratio = 1)

# For rlog transformation with PCA visualization:
rld <- rlog(atac_dds_all)
mat_atac <- assay(rld)
pca_comp <- prcomp(t(mat_atac))

# Create PCA data frame
pca_df <- as.data.frame(pca_comp$x)
variance_explained <- summary(pca_comp)$importance[2, ]

pca_df <- as.data.frame(pca_comp$x)
pca_df$DonorID <- colData(dds)$DonorID
pca_df$Condition <- colData(dds)$Condition

pca_atac_csaw <- ggplot(pca_df, aes(x = PC1, y = PC2, label = DonorID, color = Condition)) +
  geom_point(size = 3) +  # Increase point size here
  #geom_text_repel(size = 5, box.padding = unit(1, "lines")) +
  scale_color_manual(values = c("CTL" = "#1775AE", "FNF" = "#FFC200")) +
  labs(
    title = "Total peaks with all samples",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))


#-------------------------------------------------------------------------------
# Differential atac gained, lost ,  statics
# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.05)
load("/work/users/s/e/seyoun/CQTL_AI/output/old_peaks/csaw/synovium/synovium_deseq.Rdata")
res_Shrink <- res_Shrink[!is.na(res_Shrink$padj)]

normalized_counts <- counts(atacDDS, normalized=TRUE)

# Get the condition information
condition <- colData(atacDDS)$Condition

# Calculate the baseMean for each condition
baseMean_CTL <- rowMeans(normalized_counts[, condition == "CTL"])
baseMean_FNF <- rowMeans(normalized_counts[, condition == "FNF"])

res_Shrink_df <- as.data.frame(res_Shrink)
#Getting significant ATAC-signals
diff_atac_sig <- res_Shrink %>%
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)

#finiding up sig and down sig
gained <- diff_atac_sig[diff_atac_sig$log2FoldChange > 0 & diff_atac_sig$padj < 0.05,]
lost <- diff_atac_sig[diff_atac_sig$log2FoldChange < 0 & diff_atac_sig$padj < 0.05,]
static <- res_Shrink[!(res_Shrink$peakID %in% c(gained$peakID, lost$peakID)), ]


res_Shrink_df$class <- "static"
res_Shrink_df$class[res_Shrink_df$peakID %in% gained$peakID] <- "gained"
res_Shrink_df$class[res_Shrink_df$peakID %in% lost$peakID] <- "lost"

#make a bed format
library(tibble)
gained_bed <- as.data.frame(gained) |> 
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)

lost_bed <- as.data.frame(lost) |> 
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)

all_bed <- as.data.frame(res_Shrink) |> 
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)

diff_all_bed <- as.data.frame(diff_atac_sig) |> 
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)


write.table(gained_bed, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_open_Peaks.bed", 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names=FALSE)
write.table(lost_bed, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_closed_Peaks.bed", 
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)
write.table(all_bed, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_allPeaks.bed", 
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)
write.table(diff_all_bed, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_diffAllPeaks.bed", 
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)

#-------------------------------------------------------------------------------
# MA plot-atacseq---------------------------------------------------------------

atac_MAplot <- ggplot(res_Shrink_df, aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(data = subset(res_Shrink_df, padj > 0.05), 
             color = "grey60", alpha = 0.5 )+
  geom_point(data = subset(res_Shrink_df, padj <= 0.05 & log2FoldChange > 0), 
             color = "#E07653", alpha = 0.5) +
  geom_point(data = subset(res_Shrink_df, padj <= 0.05 &  log2FoldChange < 0), 
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

output_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/plot"
dir.create(output_dir, showWarnings =TRUE)

ggsave("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/plot/atacMA_withpadding.pdf", plot = last_plot(), width = 5, height = 5)
save(atac_MAplot, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/plot/atacMA_withpadding.rda")

#------------------------------------------------------------------------------
# Gene promoter defined from the TSS  (annotatepeak)
txdb <- loadDb("/work/users/s/e/seyoun/crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)
promoters(txdb, upstre)
gainedanno <- annotatePeak(gained, TxDb=txdb, tssRegion=c(-1000, 1000), flankDistance= 1000)
lostanno <- annotatePeak(lost, TxDb=txdb, tssRegion=c(-1000, 1000), flankDistance= 1000)
bg_anno <- annotatePeak(res_Shrink, TxDb=txdb,tssRegion=c(-1000, 1000), flankDistance= 1000)

save(res_Shrink,gained,lost, static,gainedanno,lostanno,bg_anno,
     file=paste0("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/plot/Deseq",'_LFC1padj05_peakAnno.Rdata'))

# Compare genomic distributions
plotAnnoBar(list(Gained=gainedanno, Lost=lostanno, background=bg_anno))
plotDistToTSS(bg_anno)

gained_ATAC_df <- as.data.frame(gainedanno)

gained_ATAC_df[grepl("Promoter", gained_ATAC_df$annotation), ] |> head()


lost_ATAC_df <- as.data.frame(lostanno)
all_ATAC_df <- as.data.frame(bg_anno)

rbind()

#-------------------------------------------------------------------------------
#seq2gene linking genomic regions to genes in many to many

gained_genes <- seq2gene(gained,
                         tssRegion = c(-1000, 1000),
                         flankDistance = 1000,
                         TxDb = txdb)
gained_genesonly <- sub(".*/(ENSG\\d+\\.\\d+).*", "\\1", gained_genes)

# I have tried with the chipseek-anno but finding promoter is better witht the genomic Features
gained_ATAC_df[grepl("Promoter", gained_ATAC_df$annotation),] |> dim()


#pr <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 1000, by = "gene")

# find overlapping promoters
overlap_gainedatac_promoters <- findOverlaps(gained, promoters)
overlap_staticatac_promoters <- findOverlaps(static, promoters)
overlap_lostatac_promoters   <- findOverlaps(lost, promoters)

# Get the genes overlapping those
promoters_at_gained_atac <- promoters[subjectHits(overlap_gainedatac_promoters)]
promoters_at_static_atac <- promoters[subjectHits(overlap_staticatac_promoters)]
promoters_at_lost_atac   <- promoters[subjectHits(overlap_lostatac_promoters)]

promoters_gained_sf_atac_geneID <-promoters_at_gained_atac$gene_id %>%
  gsub("\\..*$", "", .) %>%
  unique()

promoters_lost_sf_atac_geneID <-promoters_at_lost_atac$gene_id %>%
  gsub("\\..*$", "", .) %>%
  unique()

promoters_static_sf_atac_geneID <-promoters_at_static_atac$gene_id %>%
  gsub("\\..*$", "", .) %>%
  unique()

background_sf_genes <- all_ATAC_df$geneId %>%
  gsub("\\..*$", "", .) %>%
  unique()

all_atac_promoter_diff_genes <-c(promoters_gained_sf_atac_geneID,promoters_lost_sf_atac_geneID) |> unique()


write.table(promoters_gained_sf_atac_geneID, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/promoter_gained_sf_atac_genes.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names=FALSE)
write.table(promoters_lost_sf_atac_geneID, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/promoter_lost_sf_atac_genes.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)
write.table(background_sf_genes, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/background_sf_atac_genes.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)
write.table(all_atac_promoter_diff_genes, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/promoter_all_sf_atac_genes.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)



#Get the genes to save for the homer run
#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_homer.sh /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/promoter_gained_sf_atac_genes.txt /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/background_sf_atac_genes.txt /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/homer/open_sf_atac")
#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_homer.sh /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/promoter_lost_sf_atac_genes.txt /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/background_sf_atac_genes.txt /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/homer/close_sf_atac")
#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_homer.sh /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/promoter_all_sf_atac_genes.txt /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/background_sf_atac_genes.txt /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/homer/all_sf_atac")

#-------------------------------------------------------------------------------
#Running GO and pathway only with the promoter gene associated with the gene expression and violin plot? 
#GO and Pathway  visualization open and closed? 

upsig_go_data <- read_delim("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/homer/open_sf_atac/biological_process.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01)
upsig_go <- reduceGO(go_data,
                     category = "Opened")

## Format and write to table
upgo_table <- upsig_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(upgo_table, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/tables/OpensigATAC_promoter_GO.csv")


downsig_go_data <- 
  read_delim("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/homer/close_sf_atac/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
downsig_go <- reduceGO(downsig_go_data,
                       category = "Closed")

## Format and write to table
downgo_table <- downsig_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(downgo_table, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/tables/ClosesigATAC_promoter_GO.csv")

# Select 5 each for plotting
upsig_go_plotting <- upsig_go |> 
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("immune response", 
                           "defense response", 
                           "response to cytokine", 
                           "cell communication", 
                           "signaling")) |> 
  arrange(`-log10pval`)

downsig_go_plotting <- downsig_go |> 
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("negative regulation of cell population proliferation", 
                           "cell adhesion",
                           "negative regulation of collagen biosynthetic process", 
                           "Leydig cell differentiation",
                           "very-low-density lipoprotein particle clearance")) |> 
  arrange(`-log10pval`)

# Combine into one
go_plotting <- bind_rows(upsig_go_plotting, downsig_go_plotting)
go_plotting$parentTerm <- factor(go_plotting$parentTerm, levels = go_plotting$parentTerm) 
go_plotting$category <- factor(go_plotting$category, levels = c("Opened", "Closed"))

colorcode_atac <- c("Close"="#287C6F", "Open"="#E07653")
"#FFB6A6"
# Plot all in barplot
ATAC_GO_barplots <- ggplot(go_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 30, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 40, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 50, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 60, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 70, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 71), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 70, 10)) +
  scale_fill_manual(values = c(lighten(colorcode_atac[["Open"]],amount= 0.3), lighten(colorcode_atac[["Close"]],amount=0.3))) +
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
        strip.text = element_text(size = 8, color =  c(colorcode_atac[["Open"]],colorcode_atac[["Close"]])),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("GO Terms")

ggsave(filename = "output/Differential_analysis/plot/atac_SF_GO_barplots.pdf",
       plot = ATAC_GO_barplots, width = 5, height = 8, units = "in")
save(ATAC_GO_barplots, file = "output/Differential_analysis/plot/atac_SF_GO_barplots.rda")

#pathway------------------------------------------------------------------------

# Read in from Homer
open_sig_kegg_data <- read_delim("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/homer/open_sf_atac/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval)) |> 
  mutate(category = "Opened")

## Format and write to table
open_sig_kegg_table <- open_sig_kegg_data |> 
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(open_sig_kegg_table, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/tables/KEGG_atac_SF_open.csv")

#close
close_sig_kegg_data <- read_delim("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/homer/close_sf_atac/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval)) |> 
  mutate(category = "Closed")

## Format and write to table
close_sig_kegg_table <- close_sig_kegg_data |> 
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(close_sig_kegg_table, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/tables/KEGG_atac_SF_close.csv")

# Plot top 5 significant for each category
open_sig_kegg_plotting <- open_sig_kegg_data |> 
  filter(Term %in% c("Cytokine-cytokine receptor interaction",
                     "TNF signaling pathway", 
                     "IL-17 signaling pathway", 
                     "NOD-like receptor signaling pathway", 
                     "Rheumatoid arthritis")) |> 
  arrange(`-log10pval`)

close_sig_kegg_plotting <- close_sig_kegg_data |> 
  filter(Term %in% c("Tyrosine metabolism",
                     "Catecholamine biosynthesis, tyrosine => dopamine => noradrenaline => adrenaline", 
                     "Focal adhesion", 
                     "Hypertrophic cardiomyopathy (HCM)", 
                     "Antifolate resistance")) |> 
  arrange(`-log10pval`)

# Combine into one
kegg_plotting <- bind_rows(open_sig_kegg_plotting, close_sig_kegg_plotting)
kegg_plotting$Term <- factor(kegg_plotting$Term, levels = kegg_plotting$Term) 
kegg_plotting$category <- factor(kegg_plotting$category , levels = c("Opened","Closed"))


atac_SF_KEGG_barplots <- ggplot(kegg_plotting, aes(x = `-log10pval`, y = Term, fill = category)) +
  geom_vline(xintercept = 2, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 4, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 6, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 8, color = "grey95", alpha = 0.4) +
  geom_vline(xintercept = 10, color = "grey95", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~pval", limits = c(0, 12),
                     breaks = seq(0, 10, 2)) +
  scale_fill_manual(values = c(lighten(colorcode_atac[["Open"]],amount= 0.3), lighten(colorcode_atac[["Close"]],amount=0.3))) +
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
        strip.text = element_text(size = 8, color =  c(colorcode_atac[["Open"]],colorcode_atac[["Close"]])),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("KEGG Pathways")

ggsave(filename = "output/Differential_analysis/plot/SF_ATAC_KEGG_barplots.pdf",
       plot = atac_SF_KEGG_barplots, width = 5, height = 8, units = "in")
save(atac_SF_KEGG_barplots, file = "output/Differential_analysis/plot/SF_ATAC_KEGG_barplots.rda")


# Filter the gene from the RNA data---------------------------------------------
load("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/SF_differential_Deseq2_results.RData")
gained_rna_sf_filtered_promoter_genes <- res_Shrink_RNA_SF_df |> dplyr::filter(gene_id %in% promoters_at_gained_atac$gene_id)
lost_rna_sf_filtered_promoter_genes <- res_Shrink_RNA_SF_df |> dplyr::filter(gene_id %in% promoters_at_lost_atac$gene_id)
static_rna_sf_filtered_promoter_genes <- res_Shrink_RNA_SF_df |> dplyr::filter(gene_id %in% promoters_at_static_atac$gene_id)


gained_rna_sf_filtered_promoter_genes$atac_class = "Opened"
static_rna_sf_filtered_promoter_genes$atac_class = "Static"
lost_rna_sf_filtered_promoter_genes$atac_class = "Closed"

# combine into one df
atac_rna_df = rbind(gained_rna_sf_filtered_promoter_genes ,static_rna_sf_filtered_promoter_genes,lost_rna_sf_filtered_promoter_genes)

atac_rna_df$atac_class <- factor(atac_rna_df$atac_class,levels=c("Closed","Static","Opened"))



violoinPlot_atac <- ggplot(
  atac_rna_df, aes(x = atac_class, y = log2FoldChange, fill=atac_class,color=atac_class)) +
  geom_violin() +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.75,color="black") +
  #stat_summary(fun.data=mean_sdl,
  #             geom="linerange",color="black",size=0.25) + 
  stat_summary(fun = mean, geom = "crossbar",color="black",width = 0.1, linewidth = 0.5) +
  coord_cartesian(ylim = c(-10, 10)) +
  #ylim(-8, 8) +
  scale_fill_manual(values  = c("Static" = "grey", "Closed" = "#287C6F", "Opened" = "#E07653")) +
  scale_color_manual(values = c("Static" = "grey", "Closed" = "#287C6F", "Opened" = "#E07653")) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(color = "black", size = 6),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        axis.line.y = element_line(linewidth = 0.25),
        strip.text = element_text(size = 8, color = "black"),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
labs(x = "Chromatin accessibility at promoter",
     y =  expression("RNA Log"[2]*" Fold change FN-f / PBS"))
ggsave("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/plot/atac_violin_RNApromoter.pdf", plot = last_plot(), width = 5, height = 5)
save(violoinPlot_atac, file = "/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/plot/atac_violin_RNApromoter.rda")


#-------------------------------------------------------------------------------
#Running Motif

#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_motif_homer.sh /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_open_Peaks.bed /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_allPeaks.bed /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/open_sf_atac")
#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_motif_homer.sh /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_closed_Peaks.bed /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/SF_atac_allPeaks.bed /work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis/motif/closed_sf_atac")

motifs_up = read_tsv("output/Differential_analysis/motif/open_sf_atac/knownResults.txt")
motifs_do = read_tsv("output/Differential_analysis/motif/closed_sf_atac/knownResults.txt")

# Read in and subset
upsig_knownmotifs <- read_delim("output/Differential_analysis/motif/open_sf_atac/knownResults.txt") |> 
  # Convert percentages to numbers
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`),
                ~ gsub("%", "", .))) |> 
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`), as.numeric)) |> 
  # Calculate log2 enrichment
  mutate(log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`)) |> 
  # Calculate -log10pval
  mutate(log10pval = -log10(exp(`Log P-value`))) |> 
  slice_max(order_by = log10pval, n = 4) |> 
  mutate(motifLogo = paste0("output/Differential_analysis/motif/open_sf_atac/knownResults/known", row_number(), ".logo.svg"))





# Pull out first part of motif name
upsig_knownmotifs$Name <- unlist(lapply(str_split(upsig_knownmotifs$`Motif Name`, 
                                                  '[(]'), `[[`, 1))


motifImages <- list()
# Make motif images with name and pvalues
for (motif in unique(upsig_knownmotifs$Name)){
  
  # Get image
  motifImg <- pictureGrob(
    readPicture(rawToChar(rsvg::rsvg_svg(upsig_knownmotifs |> 
                                           filter(Name == motif) |> 
                                           pull(motifLogo) |> 
                                           unique()))),
    width = unit(1.25, "npc"),    # Set small width
    height = unit(1.25, "npc")    # Set small height
  )
  
  motifGrab <- grid.grabExpr(expr = {
    grid.newpage()
    grid.draw(motifImg)
  })
  
  motifImages[[motif]] <- motifGrab
  
}

save(motifImages, file = "output/Differential_analysis/motif/motifImages.rda")


pageCreate(width = 5,height = 4, showGuides = TRUE)
plotgardener::plotGG(motifImages$`NFkB-p65-Rel`, x = 2, y = 0.5, 
       width = 1.25, height = 1.25, just = c("right", "top"))









# Function to process motif results
process_motifs <- function(file_path, n_top = 4) {
  read_delim(file_path) |> 
    mutate(
      # Convert percentages to numbers
      across(c(`% of Target Sequences with Motif`, 
               `% of Background Sequences with Motif`),
             ~ as.numeric(gsub("%", "", .)))) |> 
    # Calculate enrichment and p-values
    mutate(
      log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`),
      log10pval = -log10(exp(`Log P-value`))
    ) |> 
    # Get top motifs
    slice_max(order_by = log10pval, n = n_top) |> 
    # Add path to logo files
    mutate(
      motifLogo = paste0(dirname(file_path), "/knownResults/known", row_number(), ".logo.svg"),
      # Extract motif family name
      Name = unlist(lapply(str_split(`Motif Name`, '[(]'), `[[`, 1))
    )
}



up_motifs <- process_motifs("output/Differential_analysis/motif/open_sf_atac/knownResults.txt")
down_motifs <- process_motifs("output/Differential_analysis/motif/closed_sf_atac/knownResults.txt")


# Function to create motif images
create_motif_images <- function(motif_data) {
  motif_images <- list()
  
  for (motif in unique(motif_data$Name)) {
    # Get current motif data
    current_motif <- motif_data |> 
      filter(Name == motif)
    
    # Create image
    motif_img <- pictureGrob(readPicture(
      rawToChar(rsvg::rsvg_svg(current_motif$motifLogo[1]))
    ),
    width = unit(1, "npc"),    # Set small width
    height = unit(1, "npc")    # Set small height
    )
    
   
    
    # Grab the image
    motif_grab <- grid.grabExpr(expr = {
      grid.newpage()
      grid.draw(motif_img)
    })
    
    
    motif_images[[motif]] <- list(
      image = motif_grab,
      pvalue = current_motif$`P-value`[1],
      log2enrichement = current_motif$log2enrichment[1],
      name = motif
    )
  }
  
  return(motif_images)
}

#Select that I like 

up_set <- c("NFkB-p65","IRF1","HLF","c-Jun-CRE")
up_motifs_sub <- up_motifs |> dplyr::filter(Name %in% up_set)
# Create images for both up and down regulated motifs
up_images <- create_motif_images(up_motifs_sub)
down_images <- create_motif_images(down_motifs)

plot_motif_set <- function(motif_images, x_start = 0.5, y_start = 0.5) {
  # Calculate offsets based on input coordinates
  x_offset <- x_start - 0.5  # Difference from original x
  y_offset <- y_start - 0.5  # Difference from original y
  
  # Add title with adjusted positions
  plotText("Motif", 
           x = 0.5 + x_offset, 
           y = 0.5 + y_offset,
           fontfamily = "Helvetica",
           fontsize = 8, 
           fontface = "bold",
           just = "left")
  
  plotText("Name", 
           x = 1.75 + x_offset, 
           y = 0.5 + y_offset,
           fontfamily = "Helvetica",
           fontsize = 8, 
           fontface = "bold",
           just = "left")
  
  plotText("p-value", 
           x = 2.5 + x_offset, 
           y = 0.5 + y_offset,
           fontfamily = "Helvetica",
           fontsize = 8, 
           fontface = "bold",
           just = "left")
  
  #plotText("log2 Enrichment", 
  #         x = 3 + x_offset, 
  #         y = 0.5 + y_offset,
  #         fontfamily = "Helvetica",
  #         fontsize = 8, 
  #         fontface = "bold",
  #         just = "left")
  
  # Plot motifs in a grid
  for(i in seq_along(motif_images)) {
    y_pos <- (0.1 + (i-1) * 0.4) + y_offset
    
    # Plot motif logo with adjusted position
    plotGG(motif_images[[i]]$image, 
           x = 0.45 + x_offset, 
           y = y_pos,
           width = 1.25, 
           height = 1.25, 
           just = c("left", "top"))
    
    # Add motif information with adjusted positions
    plotText(motif_images[[i]]$name,
             x = 1.75 + x_offset, 
             y = y_pos + 0.62,
             just = "left",
             fontfamily = "Helvetica",
             fontsize = 8)
    
    plotText(motif_images[[i]]$pvalue,
             x = 2.5 + x_offset, 
             y = y_pos + 0.62,
             just = "left",
             fontfamily = "Helvetica",
             fontsize = 8)
    
    #plotText(sprintf("%.3f", motif_images[[i]]$log2enrichement),
    #         x = 3 + x_offset, 
    #         y = y_pos + 0.62,
    #         just = "left",
    #         fontfamily = "Helvetica",
    #         fontsize = 8)
  }
}

# Usage example:
pageCreate(width = 5, height = 4, showGuides = FALSE)
plot_motif_set(up_images, x_start = 0.5, y_start = 0.5)  # Move to x=1, y=2


#------------------------------------------------------------------------------


# Making a figure for Doug's grant -AIM2
#-------------------------------------------------------------------------------
# Signal data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_AI/rna_output/signals/synovium/merged_signal/"
ctl_signal <- paste0(dir_merged,"CTL_sorted.bw")
fnf_signal <- paste0(dir_merged,"FNF_sorted.bw")

pdf(paste0("output/Differential_analysis/plot/ATACplot_v2.pdf"),width=7.25,height=7)
#ATAC---------------------------------------------------------------------------
pageCreate(width = 7.25, height =7 , default.units = "inches", showGuides =FALSE)

#A-MA-plot----------------------------------------------------------------------
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load("output/Differential_analysis/plot/atacMA_withpadding.rda")
plotGG(plot = atac_MAplot, x = 0.4, y = 0.5, height = 3.1, width = 3.5)
#A-MA-plot----------------------------------------------------------------------

plotText("B", x = 4.0, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plot_motif_set(up_images, x_start = 4.65, y_start = 0.5)
plot_motif_set(down_images, x_start = 4.65, y_start = 2.3)
plotRect(
  x=4.55,
  y=1.3,
  width=0.1,
  height=1.45,
  just = c("center"),
  default.units = "inches",
  linecolor =  "#E07653",
  lwd = 1,
  lty = 1,
  fill = "#E07653",
  alpha = 0.5)

plotText(label = "Opened", 
         x = 4.5, 
         y = 1.25, 
         rot = 90, 
         fontsize = 8, 
         just = "top",
         default.units = "inches", 
         fontfamily = "Helvetica",
         fontcolor = "black" ,
         fontface="bold")




plotRect(
  x=4.55,
  y=3.1,
  width=0.1,
  height=1.45,
  just = c("center"),
  default.units = "inches",
  linecolor =  "#287C6F",
  lwd = 1,
  lty = 1,
  fill = "#287C6F",
  alpha = 0.5)
plotText(label = "Closed", 
         x = 4.5, 
         y = 3.1, 
         rot = 90, 
         fontsize = 8, 
         just = "top",
         default.units = "inches", 
         fontfamily = "Helvetica",
         fontcolor = "black" ,
         fontface="bold")
#Signal plot--------------------------------------------------------------------

plotText("C", x = 0.1, y = 4.0, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")


# Signal indv data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_AI/output/signals/synovium"
indv_atac_ctl_signal <- list.files(dir_merged,full.names = TRUE,pattern = "CTL_synovium_1.bw")
indv_atac_fnf_signal <- list.files(dir_merged,full.names = TRUE,pattern = "FNF_synovium_1.bw")

GWAS_manhattan_plot_atac(
  gwas_data = gwas_data_subset,
  #gene_highlights = gene_hl,
  label= paste0("GWAS-", i),
  start_anno=start_atac,
  end_anno=end_atac,
  rsID = gwas_rsID,
  atac_ctl_signal = atac_ctl_signal,
  atac_fnf_signal = atac_fnf_signal,
 #ctl_signal = ctl_signal,
 #fnf_signal = fnf_signal,
  x_start = 0.8,
  y_start = 4.55,
  width = 2.5,
  height = 1,
  signal_height = 0.43,
  padding = -5e+05,
  atac_ctl_indv = indv_atac_ctl_signal,
  atac_fnf_indv = indv_atac_fnf_signal
  
)

dev.off()
