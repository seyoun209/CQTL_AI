# Two way interaction analysis 
# Load libraries
library(DESeq2)
library(data.table)
library(dplyr)
library(stringr)

setwd("/work/users/s/e/seyoun/CQTL_AI/output/")

#--------------------------------------------------------------
# Load the data counts

#se, macs2_gr, counts_macs2, meta_final,
load("diff_deseq2/condition/QC/Chon_macs2_se.RData")
#cqn.counts, GCcontent, peakwidths,cqnNormFactors,
load("diff_deseq2/condition/QC/CQN_results.RData")




#-------------------------------------------------------
# Load RASQUAL DATA for both condtions
#-------------------------------------------------------
#eigenMT results for pbs and fnf
pbs_qtl <- fread("rasqual_output/eigenMT/combined_results/1kb/pbs/pc1/eigenMT_combined_pbs_pc1.txt")
fnf_qtl <- fread("rasqual_output/eigenMT/combined_results/1kb/fnf/pc2/eigenMT_combined_fnf_pc2.txt")

#Combining 

pbs_qtl_sig <- pbs_qtl |> filter(fdr < 0.05) |> 
mutate(Condition = "pbs") 

fnf_qtl_sig <- fnf_qtl |> filter(fdr < 0.05) |> 
mutate(Condition = "fnf")

#Get it from the rasqual

pbs_rasqual <- fread("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_output/window_1kb/combined_1kb/pbs_pc1_1kb_combined.txt")
fnf_rasqual <- fread("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_output/window_1kb/combined_1kb/fnf_pc2_1kb_combined.txt")

#Get the significant rasqual
pbs_rasqual_sig <- inner_join(pbs_qtl_sig, pbs_rasqual, by = c("peak" = "Feature", "snp" = "rs_ID")) |>
mutate(Condition = "pbs")
fnf_rasqual_sig <- inner_join(fnf_qtl_sig, fnf_rasqual, by = c("peak" = "Feature", "snp" = "rs_ID")) |>
mutate(Condition = "fnf")

all_rasqual_qtl <- rbind(pbs_rasqual_sig, fnf_rasqual_sig)
all_rasqual_qtl$Peak_SNP <- paste(all_rasqual_qtl$peak, all_rasqual_qtl$snp, sep = "_")

#Making a duplicated SNPs object, can use to reference caQTLs significant in BOTH conditions
all_rasqual_qtl_duplicated <- all_rasqual_qtl[duplicated(all_rasqual_qtl[, "Peak_SNP"]),]
#Unique Peak-SNP pairings list significant in EITHER condition; those with duplicates noted above
all_rasqual_qtl <- all_rasqual_qtl[!duplicated(all_rasqual_qtl[, "Peak_SNP"]),]

#-------------------------------------------------------
# Filter counts objects 
#-------------------------------------------------------

condition <- c('pbs','fnf')
Chromosomes <- c(1:22)

response_dds <- DESeqDataSetFromMatrix(countData = counts_macs2,
                              colData = meta_final,
                              rowRanges=macs2_gr,
                              design = ~Donor)

# Apply CQN normalization factors
normalizationFactors(response_dds) <- cqnNormFactors

# Run DESeq2
response_dds <- DESeq(response_dds)
response_vsd <- vst(response_dds)
vsd_response_counts = assay(response_vsd) 
peak_regions <- as.data.frame(rowRanges(response_dds)) |> 
dplyr::select(seqnames, start, end, peakID) |>
mutate(peak= paste0(seqnames, "_", start, "_", end)) 

head(peak_regions)
counts_vsd_regions <- cbind( peak_regions, vsd_response_counts)

#Filter the counts from the significant rasqual
  counts_vsd_region_final <- counts_vsd_regions |>
                filter(peak %in% all_rasqual_qtl$peak)
#-------------------------------------------------------
# Load Genotype Data Frame over all Chromosomes
#-------------------------------------------------------
pbs_geno_matrix <- NULL
for (i in Chromosomes) {
    print(i)
    chr_file <- paste0("caQTL_rasqual_input/vcf/recode_012/pbs/recodeA_chr", i, "_pbs.traw")
    chr_file_df <- fread(chr_file)
    pbs_geno_matrix <- rbind(pbs_geno_matrix, chr_file_df)
}

pbs_geno_matrix_final <- pbs_geno_matrix |> filter(SNP %in% all_rasqual_qtl$snp) 

fnf_geno_matrix <- NULL
for (i in Chromosomes) {
    print(i)
    chr_file <- paste0("caQTL_rasqual_input/vcf/recode_012/fnf/recodeA_chr", i, "_fnf.traw")
    chr_file_df <- fread(chr_file)
    fnf_geno_matrix <- rbind(fnf_geno_matrix, chr_file_df)
}

fnf_geno_matrix_final <- fnf_geno_matrix |> filter(SNP %in% all_rasqual_qtl$snp)


# Define your key columns
key_cols <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")

# Merge the PBS and FNF genotype matrices by the key columns.
combined_geno_matrix <- full_join(pbs_geno_matrix_final, fnf_geno_matrix_final,
                                  by = key_cols)


#-------------------------------------------------------
# Covariate data
#-------------------------------------------------------
  covar_pbs <- fread("caQTL_rasqual_input/covar/pbs_covariates_pc2.txt")
  covar_fnf <- fread("caQTL_rasqual_input/covar/fnf_covariates_pc2.txt")

 meta_enhanced <- meta_final %>%
  mutate(
    Sex_numeric = as.numeric(Sex == "M"),  # 1=Male, 0=Female
    Protocol_batch = as.numeric(factor(ATACProtocolDate))  # Protocol batch as numeric factor
  )

ctl_indices <- which(meta_enhanced$Condition == "CTL")
ctl_samples <- meta_enhanced[ctl_indices, ]
fnf_indices <- which(meta_enhanced$Condition == "FNF")
fnf_samples <- meta_enhanced[fnf_indices, ]
  
covar_pbs_final <- bind_cols(covar_pbs, sampleID = ctl_samples$sampleID)
covar_fnf_final <- bind_cols(covar_fnf, sampleID = fnf_samples$sampleID)

full_covar <- rbind(covar_pbs_final, covar_fnf_final)


#-------------------------------------------------------
# Save the data
#-------------------------------------------------------
save(counts_vsd_region_final, all_rasqual_qtl, all_rasqual_qtl_duplicated,
     combined_geno_matrix, full_covar,meta_enhanced,
     file = "rasqual_output/response_qtl/caQTL_prepdata.RData")



