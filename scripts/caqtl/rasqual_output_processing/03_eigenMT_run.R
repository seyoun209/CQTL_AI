#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)

# Define parameters
base_dir <- "/work/users/s/e/seyoun/CQTL_AI"
window_type <- "100kb"  # or "100kb"
conditions <- c("pbs","fnf")
#conditions <- c("pbs")
#pcs <- c("pc0", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10",
#         "pc11", "pc12", "pc13", "pc14", "pc15", "pc16", "pc17", "pc18", "pc19", "pc20")
pcs <- c("pc0", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10","pc11", "pc12", "pc13")
chromosomes <- c(1:22)


# Create the logs directory
logs_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_output/eigenMT/logs"
if (!dir.exists(logs_dir)) {
  dir.create(logs_dir, recursive = TRUE)
}

for (cond in conditions) {
  for (pc in pcs) {
    for (chr in chromosomes) {
      QTL_input <- file.path(base_dir, "output/rasqual_output/eigenMT/inputs", window_type, cond, pc,
                             paste0("qtl_chr", chr, "_", cond, "_", pc, ".txt"))
      Gen_input <- file.path(base_dir, "output/rasqual_output/eigenMT/inputs", window_type, cond, pc,
                             paste0("genotypes_chr", chr, "_", cond, "_", pc, ".txt"))
      Genopos_input <- file.path(base_dir, "output/rasqual_output/eigenMT/inputs", window_type, cond, pc,
                                 paste0("gen.positions_chr", chr, "_", cond, "_", pc, ".txt"))
      Phepos_input <- file.path(base_dir, "output/rasqual_output/eigenMT/inputs", window_type, cond, pc,
                                paste0("phe.positions_chr", chr, "_", cond, "_", pc, ".txt"))
      out_dir <- file.path(base_dir, "output/rasqual_output/eigenMT", "results", window_type, cond, pc)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      Outputpath <- file.path(out_dir, paste0("eigenMT_chr", chr, "_", cond, "_", pc, ".txt"))

      
      # Construct the sbatch command with logs written into logs_dir
      commandline <- paste0("sbatch -p general --mem 10G -t 1-00:00:00 ",
                            "-o ", logs_dir, "/Eigen_Submit_chr", chr, "_", cond, "_", pc, ".out ",
                            "-e ", logs_dir, "/Eigen_Submit_chr", chr, "_", cond, "_", pc, ".err ",
                            "-J Eigen_Submit_chr", chr, "_", cond, "_", pc, " ",
                            "/work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/03_eigenMT_run.sh ",
                            QTL_input, " ", Gen_input, " ", Genopos_input, " ", Phepos_input, " ", Outputpath, " ", chr)
      
      cat("Submitting job:\n", commandline, "\n")
      system(commandline)
    }
  }
}
