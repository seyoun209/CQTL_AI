# Making a figure for Doug's grant
setwd("/work/users/s/e/seyoun/CQTL_AI")
source("scripts/utils/functions_rna.R")
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

