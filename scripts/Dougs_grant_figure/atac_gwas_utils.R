#functions
library(plotgardener)
library(grid)
library(RColorBrewer)

pageCreate(width = 7, height =4 , default.units = "inches", showGuides = TRUE)

# manhattan plot 
yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")

prepare_locus_plot <- function(data) {
  data |> 
    mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) |>
    mutate(LDgrp = addNA(LDgrp)) |>
    dplyr::rename(chrom = "var_chr",
                  pos = "var_from",
                  p = "nom_pval",
                  LD = "R2",
                  snp = "rsID") |>
    dplyr::select("chrom", "pos", "p", "snp", "LD", "LDgrp") |> 
    mutate(LDgrp = factor(LDgrp, levels = c(NA,"(0.8,1]","(0.6,0.8]", "(0.4,0.6]","(0.2,0.4]","(0,0.2]"), ordered = TRUE)) |>
    arrange(desc(LDgrp)) |>  # Sort by LDgrp so that higher LD values are at the end
    data.frame()
}

GWAS_manhattan_plot_atac <- function(gwas_data, gene_highlights, rsID,label,start_anno,end_anno,atac_ctl_signal,atac_fnf_signal, atac_ctl_indv, atac_fnf_indv, 
                                x_start = 0.7, y_start = 1, width = 2.5, height = 1,
                                signal_height = 0.3, padding = 10000) {
  
  # Calculate region parameters
  chrom <- gwas_data$var_chr[1]
  minregion <- min(gwas_data$var_from) - (padding)
  maxregion <- max(gwas_data$var_from) + (padding) -350000
  # Set up plot parameters
  region_pg <- pgParams(assembly = "hg38", 
                        chrom = chrom,
                        chromstart = minregion,
                        chromend = maxregion,
                        x = x_start, 
                        y = y_start, 
                        width = width, 
                        height = height)
  
  gwas_locus_plot <- prepare_locus_plot(gwas_data)
  # Calculate y-limit for GWAS plot
  gwas_ylim <- ceiling(max(log10(gwas_locus_plot$p) * -1.1))
  
  
  # GWAS Manhattan Plot
  locus_plot <- plotManhattan(
    data = gwas_locus_plot,
    params = region_pg,
    range = c(0, gwas_ylim),
    fill = colorby("LDgrp",
                   palette = colorRampPalette(c("#DD3931", "#EEA741", 
                                                "#499A53", "#98CDED", "#262C74"))),
    y = y_start, 
    height = height,
    snpHighlights = data.frame(
      snp = c(rsID),
      pch = c(24),
      cex = c(0.75),
      col = c("black")
    )
  )
  
  # Add y-axis
  annoYaxis(plot = locus_plot, 
            at = seq(0, gwas_ylim, 2),
            axisLine = TRUE, 
            fontsize = 6)
  
  # Add y-axis label
  plotText(label = "-log10(p-value)", 
           x = x_start - 0.3, 
           y = y_start + height/2, 
           rot = 90, 
           fontsize = 6, 
           just = "center",
           default.units = "inches", 
           fontfamily = "Helvetica")
  
  
  
  # Add rsID marker and label
  plotText(label =  rsID, 
           x =x_start + width/2 , y =y_start-0.05, 
           just = c("center","top"),
           fontfamily = "Helvetica", fontsize = 7,
           lineheight = 0.8,
           fontcolor = "#DD3931")
  # grid.points(x = x_start + width - 0.6, 
  #             y = y_start + (height - 0.3) / 3 + 0.15, 
  #             default.units = "native", 
  #             pch = 24,
  #             size = unit(0.5, "char"))
  # Add label
  plotText(label =  label, 
           x = x_start + 0.05, y = y_start , 
           just = c("left", "top"),
           fontfamily = "Helvetica", fontsize = 7)
  
  
  # Add ATAC signal track
  atac_signals <- plotMultiSignal(
    data = list(atac_ctl_signal, atac_fnf_signal),
    params = region_pg,
    y = y_start + height + 0.1, 
    height = signal_height, 
    linecolor = c(yl_gn_bu[8], yl_gn_bu[8]), 
    fill = c(yl_gn_bu[8], yl_gn_bu[8]),
    default.units = "inches",
    gapdistance = 0.03
  )
  
  # Add signal track labels
  plotText(label = "PBS", 
           x = x_start, 
           y = y_start + height + 0.1,
           just = c("left", "top"),
           fontfamily = "Helvetica", 
           fontsize = 7,
           lineheight = 0.8,
           fontcolor = yl_gn_bu[8])
  
  plotText(label = "FN-f", 
           x = x_start, 
           y = y_start + height + 0.13 + (signal_height - 0.03)/2,
           just = c("left", "top"),
           fontfamily = "Helvetica", 
           fontsize = 7,
           lineheight = 0.8,
           fontcolor = yl_gn_bu[8])
  
  plotText(label = "ATAC", 
           x = x_start - 0.1, 
           y = y_start + height + 0.1 + signal_height/2, 
           rot = 90, 
           fontsize = 6, 
           just = "center",
           default.units = "inches", 
           fontfamily = "Helvetica",
           fontcolor = yl_gn_bu[8])
  
  
  #RNA_signals
  
  # RNA_signals <- plotMultiSignal(
  #   data = list(ctl_signal, fnf_signal),
  #   params = region_pg,
  #   y = y_start + height + 0.2 +signal_height, 
  #   height = signal_height, 
  #   linecolor = c(yl_gn_bu[6], yl_gn_bu[6]), 
  #   fill = c(yl_gn_bu[6], yl_gn_bu[6]),
  #   default.units = "inches",
  #   gapdistance = 0.03
  # )
  # 
  # # Add signal track labels
  # plotText(label = "PBS", 
  #          x = x_start, 
  #          y = y_start + height + signal_height + 0.2,
  #          just = c("left", "top"),
  #          fontfamily = "Helvetica", 
  #          fontsize = 7,
  #          lineheight = 0.8,
  #          fontcolor = yl_gn_bu[6])
  # 
  # plotText(label = "FN-f", 
  #          x = x_start, 
  #          y = y_start + height + 0.23 + signal_height+ (signal_height - 0.03)/2,
  #          just = c("left", "top"),
  #          fontfamily = "Helvetica", 
  #          fontsize = 7,
  #          lineheight = 0.8,
  #          fontcolor = yl_gn_bu[6])
  # 
  # plotText(label = "RNA", 
  #          x = x_start - 0.1, 
  #          y = y_start + height + 0.2 +signal_height+ signal_height/2, 
  #          rot = 90, 
  #          fontsize = 6, 
  #          just = "center",
  #          default.units = "inches", 
  #          fontfamily = "Helvetica",
  #          fontcolor = yl_gn_bu[6])
  # 
  
  # Add gene track
  plotgenes <- plotGenes(
    params = region_pg, 
    y = y_start + height + signal_height +  + 0.2,
    height = 0.4,
    #geneHighlights = data.frame(
    #  "gene" = gene_highlights,
    #  "color" = "#37a7db"
    #), 
    fontsize = 6
    #geneOrder = gene_highlights
  )
  
  # Add genome label
  annoGenomeLabel(
    plot = plotgenes, 
    params = region_pg, 
    fontsize = 6, 
    y = y_start + height + signal_height + 0.65
  )
  
  #Add highlight anno
  
  annoHighlight(
    plot = plotgenes,
    chrom = gwas_data$var_chr[1],
    chromstart = start_anno, chromend = end_anno+7000,
    y = y_start, height = 2.05, just = c("center", "top"),
    default.units = "inches",
    fill = "grey90",
    linecolor = "grey90",
    alpha= 0.4
  )
  
  atacindv_x_start <- x_start + width +0.5
  atacindv_y_start <- y_start+0.2
  # Set up plot parameters
  atac_indv_pg <- pgParams(assembly = "hg38", 
                        chrom = chrom,
                        chromstart = start_anno-700,
                        chromend = end_anno+700,
                        x = atacindv_x_start, 
                        y = atacindv_y_start, 
                        width = width/2)
  
  #ATAC_indv_signals
  
  ctl_atac_signals <- plotMultiSignal(
    data = atac_ctl_indv,
    params = atac_indv_pg,
    height = signal_height*5-0.5,
    linecolor = c(yl_gn_bu[8], yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8]),
    fill = c(yl_gn_bu[8], yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8]),
    default.units = "inches",
    gapdistance = 0.03
  )
  
  ATAC_plotgenes <- plotGenes(
    params = atac_indv_pg, 
    y = y_start + signal_height*5 -0.5 + 0.25,
    height = 0.4,
    #geneHighlights = data.frame(
    #  "gene" = gene_highlights,
    #  "color" = "#37a7db"
    #), 
    fontsize = 6,draw = FALSE
    #geneOrder = gene_highlights
  )
  
  # Add genome label
  annoGenomeLabel(
    plot = ATAC_plotgenes, 
    params = atac_indv_pg, 
    fontsize = 6, 
    y = y_start + signal_height*5 -0.5 + 0.25
  )
  
  ctl_atac_signals <- plotMultiSignal(
    data = atac_fnf_indv,
    params = atac_indv_pg,
    x= atacindv_x_start +width/2 +0.5,
    height = signal_height*5-0.5,
    linecolor = c(yl_gn_bu[8], yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8]),
    fill = c(yl_gn_bu[8], yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8],yl_gn_bu[8]),
    default.units = "inches",
    gapdistance = 0.03,alpha=0.5
  )
  
  ATAC_plotgenes <- plotGenes(
    params = atac_indv_pg, 
    x= atacindv_x_start +width/2 +0.5,
    y = y_start + signal_height*5 -0.5 + 0.25,
    height = 0.4,
    #geneHighlights = data.frame(
    #  "gene" = gene_highlights,
    #  "color" = "#37a7db"
    #), 
    fontsize = 6,draw = FALSE
    #geneOrder = gene_highlights
  )
  
  # Add genome label
  annoGenomeLabel(
    plot = ATAC_plotgenes, 
    params = atac_indv_pg, 
    fontsize = 6, 
    y = y_start + signal_height*5 -0.5 + 0.25,
    x= atacindv_x_start +width/2 +0.5
  )
  

  # Add signal track labels
  plotText(label = "PBS",
           x = atacindv_x_start+width/4,
           y = y_start,
           just = c("left", "top"),
           fontfamily = "Helvetica",
           fontsize = 8,
           lineheight = 0.8,
           fontcolor = yl_gn_bu[8])

  plotText(label = "FN-f",
           x = atacindv_x_start +width/2 +0.5 +width/4,
           y = y_start,
           just = c("center", "top"),
           fontfamily = "Helvetica",
           fontsize = 8,
           lineheight = 0.8,
           fontcolor = yl_gn_bu[8])

  
  
  
  
  # Add legend
  legend_x <- x_start + width - 0.8
  legend_y <- y_start
  plotLegend(
    legend = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
    fill = c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"),
    x = legend_x,
    y = legend_y,
    width = 0.05,
    height = 0.3,
    border = FALSE,
    fontsize = 6
  )
  
  # Return the plots for potential further modification
  # return(list(manhattan = locus_plot, 
  #             signal = RNA_signals, 
  #             genes = plotgenes))
}
