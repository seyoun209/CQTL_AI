# script for filtering of RASQUAL output for eigenMT input
library(dplyr)
library(stringr)
library(stats)

#-------------------------------------------------------------------------------

conditions <- c("fnf","pbs")
chromosomes <- c(1:22, "X")
pcs <- paste0("PC", 1:10) 

dir.create("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results/combined", 
           recursive = TRUE, showWarnings = FALSE)

# Loop through conditions
for (condition in conditions) {
  print(paste("Processing condition:", condition))
  
  # Loop through PCs
  for (pc in pcs) {
    print(paste("Processing", condition, pc))
    
    # Initialize dataframe for this PC
    full_pc <- as.data.frame(matrix(0, nrow=1, ncol=26))
    colnames(full_pc) <- c("Feature", "rs ID", "Chromosome", "SNP Position", "Ref Allele", 
                          "Alt allele", "Allele Frequency", "HWE Chi_Square Statistic", 
                          "Imputation Quality (IA)", "Log_10 Benjamini-Hochberg Q-value", 
                          "Chi square statistic (2 x log Likelihood ratio)", "Effect Size", 
                          "Sequencing/mapping error rate (Delta)", "Reference allele mapping bias (Phi)", 
                          "Overdispersion", "SNP ID within the region", "No. of feature SNPs", 
                          "No. of tested SNPs", "No. of iterations for null hypothesis", 
                          "No. of iterations for alternative hypothesis", "Random location of ties", 
                          "Log likelihood of the null hypothesis", "Convergence status", 
                          "Squared correlation between prior and posterior genotypes (fSNPs)", 
                          "Squared correlation between prior and posterior genotypes (rSNP)", 
                          "PValue")
    
    # Loop through chromosomes for this PC
    for (chr in chromosomes) {
      print(paste("Processing chromosome:", chr))
      
      # Construct input filename
      input_file <- file.path("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results", 
                             condition,
                             paste0("chr", chr, "_", condition, "_", pc, "_RASQUALResults_25kb.txt"))
      
      # Check if file exists
      if (!file.exists(input_file)) {
        warning(paste("File not found:", input_file))
        next
      }
      
      # Read and process RASQUAL output
     rasqual <- try(read.delim(input_file, header = FALSE))
      
 
      # Set column names
      colnames(rasqual) <- colnames(full_pc)[1:25]
      
      # Calculate p-value
      rasqual$PValue <- pchisq(rasqual$`Chi square statistic (2 x log Likelihood ratio)`, 
                              1, lower.tail = FALSE)
      
      # Filter results
      rasqual <- subset(rasqual, `SNP Position` != "-1")
      rasqual <- subset(rasqual, 
                       `Squared correlation between prior and posterior genotypes (rSNP)` >= 0.8)
      
      # Add to full PC dataframe
      full_pc <- rbind(full_pc, rasqual)
    }
    
    # Remove the initial empty row
    full_pc <- full_pc[-1,]
    
    # Write filtered results for this PC
    output_file <- file.path("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results/combined",
                            paste0(condition, "_", pc, "_combined.txt"))
    
    write.table(full_pc, output_file, 
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
}