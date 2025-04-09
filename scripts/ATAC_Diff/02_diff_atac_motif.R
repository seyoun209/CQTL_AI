setwd("/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition")
library(tibble)
library(dplyr)
library(readr)
library(stringr)
library(grid)
library(gridExtra)
library(rsvg)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)
library(plotgardener)
library(grImport2)

# Load data

#diff_atac_sig, gained, lost, static,atac_res_Shrink_df
load("atac_diff_significant_peaks.RData") 
#atac_rna_df_df, gained_rna_filtered_promoter_genes, static_rna_filtered_promoter_genes, lost_rna_filtered_promoter_genes, de_genes_shrink
load("atac_rna_promoter_comparison.RData")

#--------------------------------------------------------------------------------
#make a bed format

gained_bed <- as.data.frame(gained) |>
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)

lost_bed <- as.data.frame(lost) |>
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)

all_bed <- as.data.frame(atac_res_Shrink_df) |>
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)

diff_all_bed <- as.data.frame(diff_atac_sig) |>
  dplyr::select(seqnames, start, end, peakID, strand) |>
  tibble::add_column(score = 0, .after = 4)


write.table(gained_bed, file = "motif/atac_gain_Peaks.bed",
            sep = "\t", row.names = FALSE, quote = FALSE, col.names=FALSE)
write.table(lost_bed, file = "motif/atac_lost_Peaks.bed",
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)
write.table(all_bed, file = "motif/atac_allPeaks.bed",
            sep = "\t", row.names = FALSE, quote = FALSE,col.names=FALSE)


#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_motif_homer.sh /work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/motif/atac_gain_Peaks.bed /work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/motif/atac_allPeaks.bed /work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/motif/gained_atac")
#system("/work/users/s/e/seyoun/CQTL_AI/scripts/run_motif_homer.sh /work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/motif/atac_lost_Peaks.bed /work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/motif/atac_allPeaks.bed /work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/motif/lost_atac")

# Function to process motif results
process_motifs <- function(file_path, n_top = 20) {
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

# Create motif images for gained and lost motifs

up_motifs <- process_motifs("motif/gained_atac/knownResults.txt")
down_motifs <- process_motifs("motif/lost_atac/knownResults.txt")

up_set <- c("NFkB-p65","Atf3","AP-1","BATF","Fos")
up_motifs_sub <- up_motifs |> dplyr::filter(Name %in% up_set)
down_set <- c("TEAD","FOXA1","NF1","Mef2d","ETV4")
down_motifs_sub <- down_motifs |> dplyr::filter(Name %in% down_set)
# Create images for both up and down regulated motifs
up_images <- create_motif_images(up_motifs_sub)
down_images <- create_motif_images(down_motifs_sub)

plot_motif_set(up_images, x_start = 4.65, y_start = 0.5)
plot_motif_set(down_images, x_start = 4.65, y_start = 2.3)

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
  # Plot motifs in a grid
  for(i in seq_along(motif_images)) {
    y_pos <- (0.1 + (i-1) * 0.27) + y_offset

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
             y = y_pos + 0.57,
             just = "left",
             fontfamily = "Helvetica",
             fontsize = 8)
    # Add motif information with adjusted positions
    plotText(motif_images[[i]]$name,
             x = 1.75 + x_offset,
             y = y_pos + 0.57,
             just = "left",
             fontfamily = "Helvetica",
             fontsize = 8)

    plotText(motif_images[[i]]$pvalue,
             x = 2.5 + x_offset,
             y = y_pos + 0.57,
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
pageCreate(width = 5, height = 4, showGuides = TRUE)
plot_motif_set(up_images, x_start = 0.5, y_start = 0.5)  # Move to x=1, y=2

