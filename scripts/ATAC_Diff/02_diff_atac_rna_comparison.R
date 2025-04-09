library(DESeq2)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
#-----------------------------------------------------------------------
# Load the necessary data
load("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/differential_gene_expression_dds.rda") # dds_gene
#atacdds,atac_res,atac_res_Shrink
load("atac_cqnnorm_deseq2_re.RData")  

# Finnd the prmoter regions. 
txdb <- loadDb("/work/users/s/e/seyoun/crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)
promoters <- promoters(txdb_genes, upstream =1500, downstream =500)

atac_res_Shrink_noNA <- atac_res_Shrink[!is.na(atac_res_Shrink$padj)]
#Getting significant ATAC-signals
diff_atac_sig <- atac_res_Shrink_noNA %>%
  dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.05)
#finiding up sig and down sig
gained <- diff_atac_sig[diff_atac_sig$log2FoldChange > 0 & diff_atac_sig$padj < 0.05,]
lost <- diff_atac_sig[diff_atac_sig$log2FoldChange < 0 & diff_atac_sig$padj < 0.05,]
static <- atac_res_Shrink_noNA[!(atac_res_Shrink_noNA$peakID %in% c(gained$peakID, lost$peakID)), ]


atac_res_Shrink_noNA$class <- "static"
atac_res_Shrink_noNA$class[atac_res_Shrink_noNA$peakID %in% gained$peakID] <- "gained"
atac_res_Shrink_noNA$class[atac_res_Shrink_noNA$peakID %in% lost$peakID] <- "lost"


# find overlapping promoters
overlap_gainedatac_promoters <- findOverlaps(gained, promoters)
overlap_staticatac_promoters <- findOverlaps(static, promoters)
overlap_lostatac_promoters   <- findOverlaps(lost, promoters)

# Get the genes overlapping those
promoters_at_gained_atac <- promoters[subjectHits(overlap_gainedatac_promoters)]
promoters_at_static_atac <- promoters[subjectHits(overlap_staticatac_promoters)]
promoters_at_lost_atac   <- promoters[subjectHits(overlap_lostatac_promoters)]


#------------------------------------------------------------
# Get the genes from the RNA-seq differential expression results

de_genes_shrink <- lfcShrink(dds_gene,
                             coef = "Condition_FNF_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")
gained_rna_filtered_promoter_genes <- de_genes_shrink |> dplyr::filter(gene_id %in% promoters_at_gained_atac$gene_id)
lost_rna_filtered_promoter_genes <- de_genes_shrink |> dplyr::filter(gene_id %in% promoters_at_lost_atac$gene_id)
static_rna_filtered_promoter_genes <- de_genes_shrink |> dplyr::filter(gene_id %in% promoters_at_static_atac$gene_id)

gained_rna_filtered_promoter_genes$atac_class = "Gained"
static_rna_filtered_promoter_genes$atac_class = "Static"
lost_rna_filtered_promoter_genes$atac_class = "Lost"

# Convert each GRanges object to a data frame
gained_df <- as.data.frame(gained_rna_filtered_promoter_genes)
static_df <- as.data.frame(static_rna_filtered_promoter_genes)
lost_df   <- as.data.frame(lost_rna_filtered_promoter_genes)

# Combine the data frames
atac_rna_df_df <- rbind(gained_df, static_df, lost_df)

# Ensure atac_class is a factor with the proper levels
atac_rna_df_df$atac_class <- factor(atac_rna_df_df$atac_class, levels = c("Lost", "Static", "Gained"))


save(atac_rna_df_df, 
gained_rna_filtered_promoter_genes, 
static_rna_filtered_promoter_genes, 
lost_rna_filtered_promoter_genes, de_genes_shrink,
     file = "atac_rna_promoter_comparison.RData")
# Violin plot of RNA log2 fold changes by ATAC class

atac_RNApromoter_violinPlot <- ggplot(
  atac_rna_df_df, 
  aes(x = atac_class, y = log2FoldChange, fill = atac_class, color = atac_class)) +
  geom_violin() +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.75, color = "black") +
  stat_summary(fun = mean, geom = "crossbar", color = "black", width = 0.1, linewidth = 0.5) +
  coord_cartesian(ylim = c(-10, 10)) +
  scale_fill_manual(values = c("Static" = "#999999", "Lost" = "#2057A7", "Gained" = "#F2BC40")) +
  scale_color_manual(values = c("Static" = "#999999", "Lost" = "#2057A7", "Gained" = "#F2BC40")) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background  = element_rect(fill = 'transparent', color = "transparent"),
        text             = element_text(family = "Helvetica"),
        legend.position  = "None",
        axis.text.y      = element_text(color = "black", size = 6),
        axis.title.y     = element_text(size = 8),
        axis.title.x     = element_text(size = 8),
        axis.text.x      = element_text(color = "black", size = 6),
        strip.background = element_blank(),
        axis.ticks       = element_blank(),
        axis.line.x      = element_line(linewidth = 0.25),
        axis.line.y      = element_line(linewidth = 0.25),
        strip.text       = element_text(size = 8, color = "black"),
        panel.spacing    = unit(0, "mm"),
        plot.title       = element_text(hjust = 0.5, face = "bold", size = 8)) +
  labs(x = "Chromatin accessibility at promoter",
       y =  expression("RNA Log"[2] * " Fold change FN-f / PBS"))

ggsave("plots/violinplot_promoter_RNA.pdf", 
plot = atac_RNApromoter_violinPlot, width= 5, height = 5)
save(atac_RNApromoter_violinPlot, 
file = "plots/atac_RNApromoter_violinPlot.rda")
