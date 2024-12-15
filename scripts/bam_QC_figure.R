#QC check
#Author Jess Byun
#date_start: December 6th 2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_AI/")
library(ATACseqQC)
library(data.table)
library(dplyr)
library(tidyr)
library(tximeta)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(GenomicAlignments)
#-------------------------------------------------------------------------------


qc_li <- list.files("output/bamQC", pattern="bamqQC", full.names = TRUE)
bamFileLabels <- gsub("_bamqQC", "", gsub("output/bamQC/", "", basename(qc_li)))

metadata <- data.frame(file_name = bamFileLabels) %>%
  separate(file_name, 
           into = c("project", "donor", "condition", "tissue", "replicate"),
           sep = "_",
           remove = FALSE) %>%
  mutate(
    is_replicate = ifelse(replicate == "replicate", TRUE, 
                          ifelse(replicate == "1", FALSE, TRUE)),
    condition = ifelse(condition == "CTL", "PBS", condition),
    tissue_type = case_when(
      tissue == "Ankle" ~ "Ankle",
      tissue == "Femur" ~ "Femur",
      tissue == "synovium" ~ "Synovium",
      TRUE ~ NA_character_
    ),
    sample_id = ifelse(is_replicate,
                       paste(donor, condition, tissue, "rep2", sep="_"),
                       paste(donor, condition, tissue, "rep1", sep="_"))
  )
#load(qc_li[1])


read_qc_data <- function(file_path) {
     qc_data <- load(file_path)
     data.frame(
         file_name = gsub("_bamqQC", "", basename(file_path)),
         total_reads = bamqc$totalQNAMEs *2,
         duplication_rate = bamqc$duplicateRate,
       mito_rate = bamqc$mitochondriaRate,
         proper_pair_rate = bamqc$properPairRate,
         nrf = bamqc$nonRedundantFraction,
         pbc1 = bamqc$PCRbottleneckCoefficient_1, 
        pbc2= bamqc$PCRbottleneckCoefficient_2
       )
   }
# Combine QC data with metadata
qc_df <- do.call(rbind, lapply(qc_li, read_qc_data)) %>%
     left_join(metadata, by = c("file_name" = "file_name"))


create_qc_plot <- function(data, metric, ylabel, tissue_filter = NULL, replicate_filter = NULL) {
  # Filter data if specified
  plot_data <- data
  if (!is.null(tissue_filter)) {
    plot_data <- plot_data[plot_data$tissue_type == tissue_filter, ]
  }
  if (!is.null(replicate_filter)) {
    plot_data <- plot_data[plot_data$is_replicate == replicate_filter, ]
  }
  
  # Check if data is empty
  if (nrow(plot_data) == 0) {
    stop("No data available for the specified filters.")
  }
  
  # Create facet formula based on filters
  facet_formula <- if (is.null(tissue_filter) && is.null(replicate_filter)) {
    tissue_type ~ is_replicate
  } else if (!is.null(tissue_filter)) {
    ~ is_replicate
  } else if (!is.null(replicate_filter)) {
    tissue_type ~ .
  } else {
    NULL
  }
  
  # Create base plot
  base_plot <- ggplot(plot_data, aes(x = sample_id, y = .data[[metric]], fill = condition)) +
    geom_bar(stat = "identity", width = 0.7) +  # Fix bar width
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Sample",
      y = ylabel,
      fill = "Condition"
    ) +
    scale_fill_manual(values = c("PBS" = "#2057A7", "FNF" = "#F2BC40")) +
    scale_y_continuous(
      breaks = seq(0, 160000000, by = 10000000),  # Breaks every 10M
      labels = scales::label_number(scale = 1e-6, suffix = "M")  # Show in millions
    )
  
  # Add faceting if needed
  if (!is.null(facet_formula)) {
    base_plot <- base_plot + facet_grid(facet_formula, scales = "free_x")
  }
  
  # Add 50M threshold line for total reads plot
  if (metric == "total_reads") {
    base_plot <- base_plot +
      geom_hline(yintercept = 50000000, linetype = "dashed", color = "black")
  }
  
  return(base_plot)
}


# Create separate plots for Femur Synovium
create_qc_plot(qc_df, "total_reads", "Total Reads (Millions)", replicate_filter = FALSE, tissue_filter = "Ankle" )
create_qc_plot(qc_df, "total_reads", "Total Reads (Millions)", replicate_filter = FALSE, tissue_filter = "Synovium")
duplication = create_qc_plot(qc_df, "duplication_rate", "Duplication Rate (%)", replicate_filter = FALSE, tissue_filter = "Ankle")
mito = create_qc_plot(qc_df, "mito_rate", "Mitochondrial Rate (%)", replicate_filter = FALSE, tissue_filter = "Ankle")
proper_pair = create_qc_plot(qc_df, "proper_pair_rate", "Proper Pair Rate (%)", replicate_filter = FALSE, tissue_filter = "Ankle")
nrf = create_qc_plot(qc_df, "nrf", "Non-Redundant Fraction (%)", replicate_filter = FALSE, tissue_filter = "Ankle")
  pbc1 = create_qc_plot(qc_df, "pbc1", "PBC1 (%)", tissue_filter = "Femur_Synovium")

