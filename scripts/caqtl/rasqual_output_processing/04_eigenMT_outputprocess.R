#!/usr/bin/env Rscript
# filepath: /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/04_eigenMT_outputprocess.R

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(grid)  # for unit()
library(httpgd)
options(stringsAsFactors = FALSE)

# Define parameters
base_dir <- "/work/users/s/e/seyoun/CQTL_AI"
window_type <- "1kb"  # or "100kb"
conditions <- c("pbs", "fnf")
pcs <- c("pc0", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10",
         "pc11", "pc12", "pc13", "pc14", "pc15", "pc16", "pc17", "pc18", "pc19", "pc20")
chromosomes <- 1:22

# Loop over conditions and PCs
for (cond in conditions) {
  for (pc in pcs) {
    pc_cov <- pc  # using same naming as in input files
    Allchr_results <- NULL
    
    # Process each chromosome
    for (chr in chromosomes) {
      File <- file.path(base_dir, "output/rasqual_output/eigenMT/results", 
                        window_type, cond, pc_cov,
                        paste0("eigenMT_chr", chr, "_", cond, "_", pc_cov, ".txt"))
      
      if (!file.exists(File)) {
        cat("Warning: File not found -", File, "\n")
        next
      }
      
      Eigen_results <- fread(File, header = TRUE)
      Allchr_results <- rbind(Allchr_results, Eigen_results)
    }
    
    if (is.null(Allchr_results)) {
      cat("No results found for condition:", cond, "and PC:", pc_cov, "\n")
      next
    }
    
    # Ensure numeric conversion
    Allchr_results$`p-value` <- as.numeric(Allchr_results$`p-value`)
    Allchr_results$BF <- as.numeric(Allchr_results$BF)
    
    # Sort by nominal p-value
    Allchr_results <- Allchr_results %>% arrange(`p-value`)
    
    # Calculate global FDR adjusted values (using BF column) 
    Allchr_results$fdr <- p.adjust(Allchr_results$BF, method = "fdr", n = nrow(Allchr_results))
    
    # Find significant SNPs (FDR < 0.05)
    sig.snp <- Allchr_results %>% filter(fdr < 0.05)
    
    if (nrow(sig.snp) > 0) {
      Pval_MAX <- max(sig.snp$`p-value`)
    } else {
      Pval_MAX <- NA
    }
    
    # Create output directory for combined results (e.g., /work/users/s/e/seyoun/CQTL_AI/output/rasqual_output/eigenMT/combined_results/1kb/pbs/pc0)
    output_dir <- file.path(base_dir, "output/rasqual_output/eigenMT/combined_results", 
                            window_type, cond, pc_cov)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write full combined results
    full_results_path <- file.path(output_dir, paste0("eigenMT_combined_", cond, "_", pc_cov, ".txt"))
    fwrite(Allchr_results, full_results_path, sep = "\t", quote = FALSE)
    
    # Write p-value cutoff file
    cutoff_file <- file.path(output_dir, paste0(cond, "_", pc_cov, "_MaxPvalueCutoff_05.txt"))
    write.table(Pval_MAX, cutoff_file, col.names = TRUE, row.names = FALSE, quote = FALSE)
    
    cat("Processed:", cond, pc_cov, "- Total results:", nrow(Allchr_results), 
        "- Significant SNPs:", nrow(sig.snp), "\n")
  }
}

cat("All eigenMT results processed and combined.\n")

#---------------------------------------------------------------------------
#Make a dot plot


# Generate a summary data frame of significant SNP counts automatically
summary_list <- list()
for (cond in conditions) {
  for (pc in pcs) {
    
    # Construct path to the combined results file created earlier
    combined_file <- file.path(base_dir, "output/rasqual_output/eigenMT/combined_results", 
                               window_type, cond, pc,
                               paste0("eigenMT_combined_", cond, "_", pc, ".txt"))
    
    # Check if the combined file exists
    if (file.exists(combined_file)) {
      df <- fread(combined_file, header = TRUE)
      # Ensure numeric conversion (in case it hasn't been done)
      df$fdr <- as.numeric(df$fdr)
      # Count significant SNPs (FDR < 0.05)
      count_sig <- nrow(dplyr::filter(df, fdr < 0.05))
    } else {
      count_sig <- NA
    }
    
    # Store summary info; converting condition to uppercase (e.g., "pbs" to "PBS")
    summary_list[[paste0(cond, "_", pc)]] <- data.frame(
      type = toupper(cond),
      pc = pc,
      counts = count_sig,
      stringsAsFactors = FALSE
    )
  }
}
# Combine all summary rows into one data frame
df_long <- do.call(rbind, summary_list)
df_long$pc <- factor(df_long$pc, levels = pcs)

# Determine the maximum counts per type for annotation
max_counts <- df_long %>% 
  group_by(type) %>% 
  filter(counts == max(counts)) %>% 
  mutate(pc_max = pc, max_count = counts,
         color = ifelse(type == "PBS", "#2057A7", "#F2BC40")) %>% 
  ungroup()

# Create the dot plot
counts_dotPlot <- ggplot(df_long, aes(x = pc, y = counts, color = type)) +
  geom_point(size = 2) +
  labs(x = "PC", y = "Number of significant introns", color = "Type") +
  scale_color_manual(values = c("PBS" = "#2057A7", "FNF" = "#F2BC40")) +
  scale_y_continuous(limits = c(500, 5000)) +
  theme_classic() +
  geom_text(data = max_counts, aes(x = pc_max, y = max_count + 100,
                                   label = paste0(pc_max, ": ", max_count)),
            color = max_counts$color, vjust = 0, size = 3.5) +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.9, 0.9),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8),
        axis.ticks.x = element_line(color = "black", linewidth = 0.1),
        axis.text.x = element_text(color = "black", size = 8),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

# Save the plot to a file and display it
save(counts_dotPlot, file = "output/plots/01_rasqual_sigSNPs/counts_dotplot.rda")
ggsave(filename = "output/plots/01_rasqual_sigSNPs/counts_dotplot.pdf",
       plot = counts_dotPlot, width = 7.5, height = 4, units = "in")
print(counts_dotPlot)