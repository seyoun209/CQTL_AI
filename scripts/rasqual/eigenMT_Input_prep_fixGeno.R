#Preparing EigenMT input - FIXED VERSION

library(data.table)
library(dplyr)
library(stringr)

# Define parameters
#conditions <- c("fnf", "pbs")
conditions <- c("pbs")
chromosomes <- c(1:22, "X")
PCs <- c(1:10)

# Set directories
base_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/results/combined"
output_base_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/eigenMT"
genotype_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf/recode_012"

# Create output directories
for (condition in conditions) {
  dir.create(file.path(output_base_dir, condition), showWarnings = FALSE, recursive = TRUE)
}

# Process each condition
for (j in 1:length(conditions)) {
  condition <- conditions[j]
  
  # Process each chromosome
  for (c in chromosomes) {
    cat(sprintf("Processing chromosome %s for condition %s\n", c, condition))
    
    # Process each PC covariate with the RASQUAL results
    for (pc in PCs) {
      pc_cov <- paste0("PC", pc)
      cat(sprintf("Processing condition: %s, chromosome: %s, with covariate: %s\n", condition, c, pc_cov))
      
      # Read RASQUAL results
      rasqual_file <- file.path(base_dir, paste0(condition, "_", pc_cov, "_combined.txt"))
      
      if (!file.exists(rasqual_file)) {
        cat(sprintf("Warning: RASQUAL file not found: %s\n", rasqual_file))
        next
      }
      
      cat(sprintf("Reading RASQUAL file: %s\n", rasqual_file))
      rasqual_re <- fread(rasqual_file, header = TRUE)
      
      # Filter for current chromosome
      chr_pattern <- paste0("chr", c)
      rasqual_re_chr <- rasqual_re %>% 
        filter(Chromosome == chr_pattern)
      
      if (nrow(rasqual_re_chr) == 0) {
        cat(sprintf("Warning: No data for chromosome %s in file %s\n", c, rasqual_file))
        next
      }
      
      # Construction of Gen.position.txt
      Gen.pos <- rasqual_re_chr %>% 
        mutate(Chromosome = gsub("chr", "", Chromosome)) %>% 
        select(`rs ID`, Chromosome, `SNP Position`) %>%
        rename(snp = `rs ID`, pos = `SNP Position`, chr_snp = Chromosome) %>%
        filter(snp != "SKIPPED")  # Remove SKIPPED entries
      
      # Write gen.positions.txt
      gen.position.output <- file.path(output_base_dir, condition, 
                                       paste0("gen.positions_chr", c, "_", condition, "_", pc_cov, ".txt"))
      write.table(Gen.pos, gen.position.output, quote = FALSE, row.names = FALSE)
      cat(sprintf("Saved gen.positions file: %s\n", gen.position.output))
      
      # Construction of Phe.position.txt
      Phe.pos <- rasqual_re_chr %>% 
        select(Feature) %>%
        # Extract components using regex
        mutate(
          chrom_probe = sub("chr([^_]+)_.*", "\\1", Feature),
          coords = sub("chr[^_]+_(.*)", "\\1", Feature)
        ) %>%
        # Split coordinates - FIXED to handle all formats
        mutate(
          s1 = gsub("_.*$|[-].*$", "", coords),  # Get first number before _ or -
          s2 = gsub("^.*[-_]", "", coords)       # Get last number after _ or -
        ) %>%
        # Clean up and rename
        select(-coords) %>%
        select(Feature, chrom_probe, s1, s2) %>%
        rename(peak_id = Feature)
      
      # Remove duplicates
      Phe.pos <- distinct(Phe.pos, peak_id, .keep_all = TRUE)
      
      # Write phe.positions.txt
      phe.position.output <- file.path(output_base_dir, condition, 
                                       paste0("phe.positions_chr", c, "_", condition, "_", pc_cov, ".txt"))
      write.table(Phe.pos, phe.position.output, quote = FALSE, row.names = FALSE)
      cat(sprintf("Saved phe.positions file: %s\n", phe.position.output))
      
      # Construction of qtl.txt
      qtl <- rasqual_re_chr %>% 
        select(`rs ID`, Feature, PValue) %>%
        rename(snp = `rs ID`, peak = Feature, `p-value` = PValue) %>%
        filter(snp != "SKIPPED")  # Remove SKIPPED entries
      
      # Write qtl.txt
      qtl_output <- file.path(output_base_dir, condition, 
                              paste0("qtl_chr", c, "_", condition, "_", pc_cov, ".txt"))
      write.table(qtl, qtl_output, quote = FALSE, row.names = FALSE)
      cat(sprintf("Saved qtl file: %s\n", qtl_output))
      
      # Now process genotypes AFTER we have the QTL file
      # This ensures we only include SNPs that are actually in the QTL file
      geno_file_path <- file.path(genotype_dir, condition, paste0("recodeA_chr", c, "_", condition, ".traw"))
      
      if (!file.exists(geno_file_path)) {
        cat(sprintf("Warning: Genotype file not found: %s\n", geno_file_path))
        next
      }
      
      cat(sprintf("Reading genotype file: %s\n", geno_file_path))
      
      # Get the list of SNPs from the QTL file
      qtl_snps <- unique(qtl$snp)
      cat(sprintf("Found %d unique SNPs in QTL file\n", length(qtl_snps)))
      
      # Read genotype file in chunks to avoid memory issues
      chunk_size <- 10000
      genotypes_matched <- NULL
      header_read <- FALSE
      header <- NULL
      
      # Open the file for reading
      con <- file(geno_file_path, "r")
      
      # First read the header
      header <- readLines(con, n = 1)
      header_parts <- unlist(strsplit(header, "\t"))
      id_col_idx <- 2  # CHROM column in traw format
      
      # Prepare column indices for the data we want
      cols_to_keep <- c(id_col_idx, 7:length(header_parts))
      
      # Initialize counter
      total_lines <- 0
      matched_lines <- 0
      
      cat("Starting to process genotype file in chunks...\n")
      
      # Process the file in chunks
      repeat {
        chunk <- readLines(con, n = chunk_size)
        if (length(chunk) == 0) break
        
        # Process chunk
        chunk_data <- lapply(chunk, function(line) {
          parts <- unlist(strsplit(line, "\t"))
          if (length(parts) >= max(cols_to_keep)) {
            snp_id <- parts[id_col_idx]
            if (snp_id %in% qtl_snps) {
              return(parts[cols_to_keep])
            }
          }
          return(NULL)
        })
        
        # Remove NULL entries
        chunk_data <- chunk_data[!sapply(chunk_data, is.null)]
        
        # Count
        total_lines <- total_lines + length(chunk)
        matched_lines <- matched_lines + length(chunk_data)
        
        # Add matched data
        if (length(chunk_data) > 0) {
          if (is.null(genotypes_matched)) {
            genotypes_matched <- do.call(rbind, lapply(chunk_data, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
          } else {
            genotypes_matched <- rbind(genotypes_matched, 
                                       do.call(rbind, lapply(chunk_data, function(x) as.data.frame(t(x), stringsAsFactors = FALSE))))
          }
        }
        
        # Progress update
        cat(sprintf("Processed %d lines, found %d matches so far...\n", total_lines, matched_lines))
      }
      
      # Close the connection
      close(con)
      
      cat(sprintf("Completed genotype processing. Found %d matching SNPs out of %d total.\n", 
                  matched_lines, total_lines))
      
      # If we found matches, process and save them
      if (!is.null(genotypes_matched) && nrow(genotypes_matched) > 0) {
        # Fix column names
        colnames(genotypes_matched) <- c("ID", header_parts[7:length(header_parts)])
        
        # Clean up duplicated sample names if they exist
        if (any(grepl("_CQTL_", names(genotypes_matched)[-1]))) {
          new_colnames <- c("ID", sapply(strsplit(names(genotypes_matched)[-1], "_CQTL_"), `[`, 1))
          names(genotypes_matched) <- new_colnames
        }
        
        # Save chromosome-specific genotype file
        genotypes_output <- file.path(output_base_dir, condition, paste0("genotypes_chr", c, "_", condition, "_", pc_cov, ".txt"))
        write.table(genotypes_matched, genotypes_output, quote = FALSE, row.names = FALSE)
        cat(sprintf("Saved matched genotype file: %s\n", genotypes_output))
      } else {
        cat("WARNING: No matching SNPs found in genotype file! eigenMT will fail.\n")
        cat("You need to ensure your genotype file contains the SNPs referenced in your QTL file.\n")
      }
      
      cat(sprintf("Completed processing chr%s with PC%s for %s\n", c, pc, condition))
    }
  }
  
  cat(sprintf("Completed processing condition: %s\n", condition))
}

cat("All processing complete. Please check the output files.\n")
