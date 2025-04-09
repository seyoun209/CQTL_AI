#Visualize the differentially accessible peaks from the ATAC-seq analysis
library(ggplot2)
library(ggrepel)
library(data.table)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(DESeq2)
library(httpgd)
library(ChIPseeker)
library(plotgardener)

setwd("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition")
#atacdds,atac_res,atac_res_Shrink
load("atac_cqnnorm_deseq2_re.RData")  
#vsd, rld,normCounts
load("QC/Chon_macs2_normalized.RData") 
#se, macs2_gr, counts_macs2, meta_final
load("QC/Chon_macs2_se.RData")

#--------------------------------------------------------------------
## Finding the significant peaks
atac_res_Shrink_noNA <- atac_res_Shrink[!is.na(atac_res_Shrink$padj)]
atac_res_Shrink_df <- as.data.frame(atac_res_Shrink_noNA)
#Getting significant ATAC-signals
diff_atac_sig <- atac_res_Shrink_df %>%
  dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.05)

#finiding up sig and down sig
gained <- diff_atac_sig[diff_atac_sig$log2FoldChange > 0 & diff_atac_sig$padj < 0.05,]
lost <- diff_atac_sig[diff_atac_sig$log2FoldChange < 0 & diff_atac_sig$padj < 0.05,]
static <- atac_res_Shrink_df[!(atac_res_Shrink_df$peakID %in% c(gained$peakID, lost$peakID)), ]


atac_res_Shrink_df$class <- "static"
atac_res_Shrink_df$class[atac_res_Shrink_df$peakID %in% gained$peakID] <- "gained"
atac_res_Shrink_df$class[atac_res_Shrink_df$peakID %in% lost$peakID] <- "lost"

# Save the significant peaks
save(diff_atac_sig, gained, lost, static,atac_res_Shrink_df,
     file = "atac_diff_significant_peaks.RData")

#--------------------------------------------------------------------
# 1. Heatmap of significant peaks

load("atac_diff_significant_peaks.RData") # diff_atac_sig, gained, lost, static, atac_res_Shrink_df
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Extract normalized counts for significant peaks
rownames(normCounts) <- rowData(se)$peakID  # Ensure rownames match peak IDs
# Add sequence range info from the SummarizedExperiment 'se' to normCounts
range_info <- data.frame(
  seqnames = as.character(seqnames(rowRanges(se))),
  start = start(rowRanges(se)),
  end = end(rowRanges(se))
)
rownames(range_info) <- rowData(se)$peakID
normCounts <- cbind(range_info[rownames(normCounts), ], normCounts)

norm_counts_sig <- normCounts[rownames(normCounts) %in% diff_atac_sig$peakID, , drop = FALSE]
# Extract only the numeric part and calculate z-scores by row
num_matrix <- as.matrix(norm_counts_sig[, -(1:3)])
zscore_matrix <- t(apply(num_matrix, 1, cal_z_score))

# Convert annotation columns to factors and set levels explicitly
sample_meta <- meta_final %>%
  filter(sampleID %in% colnames(zscore_matrix)) %>%
  dplyr::select("sampleID", "Condition", "Sex", "AgeGroup")

# Set the rownames to match the colnames of zscore_matrix
rownames(sample_meta) <- sample_meta$sampleID

# Ensure factors have correct levels
sample_meta$Condition <- factor(sample_meta$Condition, levels = c("CTL", "FNF"))
sample_meta$Sex       <- factor(sample_meta$Sex, levels = c("M", "F"))
sample_meta$AgeGroup  <- factor(sample_meta$AgeGroup, levels = c("25-44", "45-64", "65-84"))

# Define the color mapping for all factor levels
my_colour <- list(
  Condition = c("CTL" = "#2057A7", "FNF" = "#F2BC40"),
  Sex       = c("M" = "#c864da", "F" = "#55e0bd"),
  AgeGroup  = c("25-44" = "#b5d1ae", "45-64" = "#80ae9a", "65-84" = "#568b87")
)

#age_palette <- c("#b5d1ae","#80ae9a","#568b87","#326b77","#1b485e","#122740")
#age_palette <- c("#80ae9a","#326b77","#122740")

# Prepare colors for the ComplexHeatmap
heatmapColors <- colorRampPalette(c("#097EA4", "black", "#BFA527"))(5)

ID <- colnames(zscore_matrix)
meta_ctl_fnf <- sample_meta %>% 
  arrange(factor(sampleID, levels = ID)) %>%
  arrange(Condition) 
  
sample_INFO <- meta_ctl_fnf |>dplyr::select("AgeGroup", "Sex", "Condition")

colAnn <- HeatmapAnnotation(
  df = sample_INFO,   
  col = my_colour,
  gap = unit(0.5, 'mm'),
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    AgeGroup = list(
      title = "AgeGroup",
      title_position = "leftcenter",
      title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(1, "mm"),
      at = c("25-44", "45-64", "65-84"),
      labels = c("25-44", "45-64", "65-84"),
      ncol = 1
    ),
    Condition = list(
      title = "Condition",
      title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 6),
      title_position = "leftcenter",
      grid_width = unit(2, "mm"),
      grid_height = unit(1, "mm"),
      at = c("CTL", "FNF"),
      labels = c("PBS", "FNF"),
      ncol = 1
    ),
    Sex = list(
      title = "Sex",
      title_position = "leftcenter",
      title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(1, "mm"),
      at = c("M", "F"),
      labels = c("Male", "Female"),
      ncol = 1
    )
  )
)

column_order <- meta_ctl_fnf$sampleID  # meta_ctl_fnf must have an "ID" column
zscore_matrix_ordered <- zscore_matrix[, column_order, drop = FALSE]

# Create a ComplexHeatmap
hmap <- Heatmap(
  zscore_matrix_ordered,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # Disable column clustering
  column_order = column_order,  # Specify column order
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE,
  top_annotation = colAnn,
  col = colorRamp2(seq(-2, 2), heatmapColors)
)

# Create a legend for the heatmap
heatmapLegend <- Legend(at = c(-2, 2),
                        col_fun = colorRamp2(breaks = seq(-2, 2),
                                             colors = heatmapColors),
                        border = NA,
                        title_gp = gpar(fontsize = 0),
                        labels_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
                        legend_width = unit(4.325, "in"),
                        grid_height = unit(0.11, "in"),
                        direction = "horizontal")

# Grab plot grobs for further use or saving
heatmapGrob <- grid.grabExpr(draw(hmap,
                                  show_annotation_legend = FALSE,
                                  show_heatmap_legend = FALSE,
                                  background = "transparent"))
heatmapLegendGrob <- grid.grabExpr(draw(heatmapLegend))

# Save the resulting grobs to file
save(heatmapGrob, file = "plots/heatmapGrob.rda")
save(heatmapLegendGrob, file = "plots/heatmapLegendGrob.rda")

#---------------------------------------------------------------------
# Finding the 
#---------------------------------------------------------------------
pdf("plots/fig1_atac.pdf",width=13,height=6,bg="transparent")  # Set background to transparent
pageCreate(width = 13, height = 6, showGuides = FALSE)  
plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

  load("plots/heatmapGrob.rda")  # Load the heatmap grob
  load("plots/heatmapLegendGrob.rda")  # Load the heatmap legend grob
plotGG(plot = heatmapGrob, x = 0.5, y = 0.5, height = 4.5, width = 5)
plotGG(plot = heatmapLegendGrob, x = 0.6, y = 5,
       width = 4.325, height = 0.11)
# Colorbar title
plotText(label = "Relative chromatin accessibility (FN-f/PBS)", fontfamily = "Helvetica",
         fontsize = 8, x = 2.75, y = 5.2, just = "top")

# Use the same colors from my_colour:
age_colors <- my_colour$AgeGroup
sex_colors <- my_colour$Sex
cond_colors <- my_colour$Condition

x_start <- 5.6
y_start <- 0.65
box_height <- 0.1    
box_width  <- 0.05 
# --- AgeGroup Legend with no vertical gap ---
ageText <- names(age_colors)  # "25-44", "45-64", "65-84"
for(i in seq_along(ageText)){
  # Place the rectangle: stack boxes with no gap
  plotRect(x = unit(x_start, "in"), 
           y = unit(y_start + (i - 1) * box_height, "in"), 
           width = unit(box_width, "in"), 
           height = unit(box_height, "in"), 
           linecolor = NA, 
           fill = age_colors[ageText[i]])
  # Place the text shifted upward by 0.05 inches
  plotText(label = ageText[i], 
           x = unit(x_start + 0.25, "in"), 
           y = unit(y_start + (i - 1) * box_height + box_height/2 - 0.05, "in"),
           fontsize = 8, fontfamily = "Helvetica")
}

x_start <- 5.6
y_start <- 1.1
box_height <- 0.1    
box_width  <- 0.05 

# --- Sex Legend ---
sexText <- names(sex_colors)  
for(i in seq_along(sexText)){
  # Place the rectangle: stack boxes with no gap
  plotRect(x = unit(x_start, "in"), 
           y = unit(y_start + (i - 1) * box_height, "in"), 
           width = unit(box_width, "in"), 
           height = unit(box_height, "in"), 
           linecolor = NA, 
           fill = sex_colors[sexText[i]])
  # Place the text shifted upward by 0.05 inches
  plotText(label = sexText[i], 
           x = unit(x_start + 0.25, "in"), 
           y = unit(y_start + (i - 1) * box_height + box_height/2 - 0.05, "in"),
           fontsize = 8, fontfamily = "Helvetica")
}

  x_start <- 5.6
  y_start <- 1.35
  box_height <- 0.1    
  box_width  <- 0.05 

# --- Condition Legend ---
condText <- c( "PBS","FN-f")  
for(i in seq_along(condText)){
  # Place the rectangle: stack boxes with no gap
  plotRect(x = unit(x_start, "in"), 
           y = unit(y_start + (i - 1) * box_height, "in"), 
           width = unit(box_width, "in"), 
           height = unit(box_height, "in"), 
           linecolor = NA, 
           fill = cond_colors[i])
  # Place the text shifted upward by 0.05 inches
  plotText(label = condText[i], 
           x = unit(x_start + 0.25, "in"), 
           y = unit(y_start + (i - 1) * box_height + box_height/2 - 0.05, "in"),
           fontsize = 8, fontfamily = "Helvetica")
}

#----------------------------------------------------------------------
# Figure 1b  Motif analysis
plotText("b", x = 6.35, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plot_motif_set(up_images, x_start = 7, y_start = 0.5)
plot_motif_set(down_images, x_start = 7, y_start = 2.1)

plotRect(
  x = 6.9,  
  y = 1.25,
  width = 0.15,
  height = 1.3,
  just = c("center"),
  default.units = "inches",
  linecolor = "#F2BC40",
  lwd = 1,
  lty = 1,
  fill = "#F2BC40",
  alpha = 0.5
)

plotText(label = "Gained",
         x = 6.85,  
         y = 1.25,
         rot = 90,
         fontsize = 8,
         just = "top",
         default.units = "inches",
         fontfamily = "Helvetica",
         fontcolor = "black",
         fontface = "bold")

plotRect(
  x = 6.9,
  y = 2.85,
  width = 0.15,
  height = 1.3,
  just = c("center"),
  default.units = "inches",
  linecolor = "#2057A7",
  lwd = 1,
  lty = 1,
  fill = "#2057A7",
  alpha = 0.5
)

plotText(label = "Lost",
         x = 6.85, 
         y = 2.9,
         rot = 90,
         fontsize = 8,
         just = "top",
         default.units = "inches",
         fontfamily = "Helvetica",
         fontcolor = "black",
         fontface = "bold")

#---------------------------------------------------------------------
# viologin plot for the ATACseq promoter and RAN gene  
#---------------------------------------------------------------------
plotText("c", x = 9.7, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load("plots/atac_RNApromoter_violinPlot.rda")

plotGG(plot = atac_RNApromoter_violinPlot,
 x = 10.2, y = 0.5, height = 3, width = 2.75)

dev.off()
