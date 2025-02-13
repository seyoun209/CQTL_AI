# siganl check
## Author: Seyoun Byun
## Date: 12.19.2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_AI/output/signals/synovium")
library(dplyr)
library(plotgardener)
library(RColorBrewer)
library(rtracklayer)
yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")

dir_merged <- "/work/users/s/e/seyoun/CQTL_AI/output/signals/synovium"
# Create vectors of file paths
ctl_files <- c(
  "FACT_23-0025_CTL_synovium_1.bw",
  "FACT_23-0016_CTL_synovium_1.bw",
  "FACT_23-0023_CTL_synovium_1.bw",
  "FACT_23-0022_CTL_synovium_1.bw",
  "FACT_23-0013_CTL_synovium_1.bw"
)

fnf_files <- c(
  "FACT_23-0025_FNF_synovium_1.bw",
  "FACT_23-0023_FNF_synovium_1.bw",
  "FACT_23-0013_FNF_synovium_1.bw",
  "FACT_23-0016_FNF_synovium_1.bw",
  "FACT_23-0022_FNF_synovium_1.bw"
)

# Create full paths
all_files <- c(ctl_files, fnf_files)
all_paths <- file.path(dir_merged, all_files)


create_gene_params <- function(chr, start, end) {
  pgParams(
    chrom = chr,
    chromstart = start - 20000,
    chromend = end + 20000,
    assembly = "hg38",
    resolution = 100,
    width = 4.0,
    x = 0.5,
    y = 0.35
  )
}

# Create a page for 2x2 layout
pdf(file = "/work/users/s/e/seyoun/CQTL_AI/output/signals/synovium/plot/signal_plot_v1.pdf",   # The directory you want to save the file in
    width = 9.75, # The width of the plot in inches
    height = 11.5)
pageCreate(width = 9.75, height = 11.5, default.units = "inches", showGuides = FALSE)

# Modified panel labeling function
get_panel_label <- function(pos) {
  switch(paste(pos$x, pos$y),
         "0.1 0.5" = "A",
         "5 0.5" = "B",
         "0.1 5.75" = "C",
         "5 5.75" = "D",
         "Unknown")
}

# Function to plot a single panel
plot_gene_panel <- function(gene_params, position, tracks_paths) {
  track_height <- 0.3
  track_spacing <- 0.2
  
  # Plot CTL tracks
  for(i in 1:5) {
    plotSignal(
      data = tracks_paths[i],
      params = gene_params,
      x = position$x,
      y = position$y + (i-1) * (track_height + track_spacing),
      height = track_height,
      width = 4.0,
      linecolor = yl_gn_bu[3],
      fill = yl_gn_bu[3],
      default.units = "inches",
      label = gsub("_synovium_1.bw", "", basename(tracks_paths[i])),
      fontsize = 6
    )
  }
  
  # Plot FNF tracks
  for(i in 6:10) {
    plotSignal(
      data = tracks_paths[i],
      params = gene_params,
      x = position$x,
      y = position$y + (i-1) * (track_height + track_spacing),
      height = track_height,
      width = 4.0,
      linecolor = yl_gn_bu[6],
      fill = yl_gn_bu[6],
      default.units = "inches",
      label = gsub("_synovium_1.bw", "", basename(tracks_paths[i])),
      fontsize = 7
    )
  }
  
  # Add gene annotation
  plotGenes(
    params = gene_params,
    x = position$x,
    y = position$y + 10 * (track_height + track_spacing),
    height = 0.3,
    width = 4.0,
    default.units = "inches",
    fontsize = 6
  )
  
  # Add panel label using the new function
  plotText(
    label = get_panel_label(position),
    x = position$x + 4.2,
    y = position$y,
    default.units = "inches",
    fontsize = 6,
    fontface = "bold"
  )
}

# Create parameters for all genes
gene_params <- list(
  cxcl8 = create_gene_params("chr4", 73740506, 73743716),
  mmp3 = create_gene_params("chr11", 102706690, 102714615),
  mmp13 = create_gene_params("chr11", 102813721, 102826461),
  il6 = create_gene_params("chr7", 22725889, 22732002)
)

# Define panel positions
panel_positions <- list(
  A = list(x = 0.1, y = 0.35),
  B = list(x = 5.0, y = 0.35),
  C = list(x = 0.1, y = 5.75),
  D = list(x = 5.0, y = 5.75)
)

# Plot all panels
plot_gene_panel(gene_params$cxcl8, panel_positions$A, all_paths)
plot_gene_panel(gene_params$mmp3, panel_positions$B, all_paths)
plot_gene_panel(gene_params$mmp13, panel_positions$C, all_paths)
plot_gene_panel(gene_params$il6, panel_positions$D, all_paths)

dev.off()
