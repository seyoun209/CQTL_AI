# EigenMT output processing
setwd("/work/users/s/e/seyoun/CQTL_AI/")
library(data.table)
library(dplyr)
library(stringr)

# Define parameters
Conditions <- c("fnf", "pbs")
Chromosomes <- c(1:22)
PCs <- c(1:10)


##############
# Process each condition and PC #
##############

for (condition in Conditions) {
  for (pc in PCs) {
    pc_cov <- paste0("PC", pc)
    
    # Initialize results dataframe
    Allchr_results <- NULL
    
    # Process each chromosome
    for (chr in Chromosomes) {
      # Construct file path 
      File <- file.path("/work/users/s/e/seyoun/CQTL_AI/output/eigenMT/results", 
                        condition, 
                        pc_cov, 
                        paste0("eigenMT_chr", chr, "_", condition, "_", pc_cov, ".txt"))
      
      # Check if file exists
      if (!file.exists(File)) {
        cat("Warning: File not found -", File, "\n")
        next
      }
      
      # Read eigenMT results
      Eigen_results <- fread(File, header = TRUE)
      
      # Combine results
      Allchr_results <- rbind(Allchr_results, Eigen_results)
    }

    
    # Ensure numeric columns
    Allchr_results$`p-value` <- as.numeric(Allchr_results$`p-value`)
    Allchr_results$BF <- as.numeric(Allchr_results$BF)
    
    # Sort by p-value
    Allchr_results <- Allchr_results[order(Allchr_results$`p-value`),]
    
    # Create output directory if it doesn't exist
    output_dir <- file.path("/work/users/s/e/seyoun/CQTL_AI/output/eigenMT/combined_results", 
                            condition, 
                            pc_cov)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write full results
    full_results_path <- file.path(output_dir, 
                                   paste0("Fullresults_", condition, "_", pc_cov, ".txt"))
    write.table(Allchr_results, 
                full_results_path, 
                row.names = FALSE, 
                col.names = TRUE)
    
    # Calculate global FDR
    Allchr_results$global <- p.adjust(Allchr_results$BF, 
                                      method = "fdr", 
                                      n = nrow(Allchr_results))
    
    # Find significant SNPs
    sig.snp <- Allchr_results[which(Allchr_results$global < 0.05),]
    
    # Get maximum p-value for significant SNPs
    if (nrow(sig.snp) > 0) {
      Pval_MAX <- max(sig.snp$`p-value`)
    } else {
      Pval_MAX <- NA
    }
    
    # Write p-value threshold
    pval_threshold_path <- file.path(output_dir, 
                                     paste0(condition, "_", pc_cov, "_MaxPvalueCutoff_05.txt"))
    write.table(Pval_MAX, 
                pval_threshold_path, 
                col.names = TRUE, 
                row.names = FALSE)
    
    # Print progress
    cat("Processed:", condition, pc_cov, 
        "- Total results:", nrow(Allchr_results), 
        "- Significant SNPs:", nrow(sig.snp), "\n")
  }
}

cat("All eigenMT results processed and combined.\n")
