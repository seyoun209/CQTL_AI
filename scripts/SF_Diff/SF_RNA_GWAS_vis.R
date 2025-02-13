#SF- GWAS data finding anything high LD with it. 
library(GenomicRanges)
library(data.table)
library(AnnotationDbi)
library(tidyr)
library(dplyr)
#-------------------------------------------------------------------------------
#Make a range for the RNA gene
load("/work/users/s/e/seyoun/CQTL_AI/rna_output/Differential_analysis/SF_differential_Deseq2_results.RData") #res_Shrink_RNA_SF_df_hgnc, diff_rna_sig, results,up, down, static, 

SF_diffRNA_sig_df <- res_Shrink_RNA_SF_df_hgnc |> filter(class != "static")
SF_diffRNA_gr_coord <- makeGRangesFromDataFrame(SF_diffRNA_sig_df)
#Make it grange
SF_diffRNA_sig_gr <- GRanges(
  seqnames = SF_diffRNA_sig_df$seqnames,
  ranges = IRanges(start = SF_diffRNA_sig_df$start, end = SF_diffRNA_sig_df$end),
  strand = SF_diffRNA_sig_df$strand,
  mcols = SF_diffRNA_sig_df[, !colnames(SF_diffRNA_sig_df) %in% c("seqnames", "start", "end", "strand")]
)

SF_diffRNA_sig_promoter_gr <- SF_diffRNA_sig_gr %>%
  mutate(start = ifelse(strand == "+", start - 2000, start),
         end = ifelse(strand == "-", end + 2000, end))



#-------------------------------------------------------------------------------
# GWAS data 


OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA",
                "SpineOA", "THR", "ThumbOA", "TJR", "TKR")

OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA",
                "THR", "ThumbOA", "TJR", "TKR")

for (i in OAsubtypes) {
  print(i)
  pdf(paste0("rna_output/Differential_analysis/plot/signal_plot_",i,".pdf"),width=4,height=4)
  subtype_leads <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/",
                                i,
                                "/leads/EUR_",i,"_leads_ld_final.csv")) |>
    dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
    dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "ldbuddy_rsID","rsID","CHR:hg38POS", "p")|>
    dplyr::mutate(
      chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
      pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2])
    ) |> dplyr::rename(gwas_rsID = "rsID") |>
    dplyr::rename(R2 = "ldbuddy_R2", var_chr = "chr_ldbuddy", var_from ="pos_ldbuddy",rsID="ldbuddy_rsID",nom_pval ="p")
  
  
  GWAS_subtype_leads <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/",
                                     i,
                                     "/leads/EUR_",i,"_leads_ld_final.csv"))
  
  LDs <- separate(subtype_leads, `ldbuddy_CHR:hg38POS`, c("LD_chr", "LD_pos"), ":")
  LDs$LD_pos <- as.numeric(LDs$LD_pos)
  LDs$LD_chr <- paste0("chr",LDs$LD_chr)
  LDs_gr <- with(LDs[LDs$R2 > 0.7,], GRanges(LD_chr, IRanges(LD_pos, width=1)))
  
  overlaps <- findOverlaps(SF_diffRNA_sig_promoter_gr, LDs_gr)
  if(length(overlaps) == 0){next}
  first_hits <- data.frame(
    queryHits = queryHits(overlaps),
    subjectHits = subjectHits(overlaps)
  ) %>%
    group_by(queryHits) %>%
    slice(1) %>%
    ungroup()
  overlapping_genes_promoter <- SF_diffRNA_sig_promoter_gr[queryHits(overlaps) |> unique()]
 
  for (j in 1:length(first_hits$queryHits)){
    print(j)
  overlapping_genes_gwasLD <- LDs_gr[first_hits$subjectHits[j]]
  matching_GWASrows <- subtype_leads[subtype_leads$`ldbuddy_CHR:hg38POS` %in% gsub("chr","",as.character(overlapping_genes_gwasLD)),]
  overlapping_genes_promoter_Row <- SF_diffRNA_sig_promoter_gr[first_hits$queryHits[j]]
  gwas_rsID <- unique(matching_GWASrows$gwas_rsID)
  gene_hl <- overlapping_genes_promoter_Row$mcols.SYMBOL

  #reading all the related GWAS
  gwas_data_subset <- GWAS_subtype_leads |> dplyr::filter(rsID  == gwas_rsID) |>
    dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
    dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "ldbuddy_rsID","rsID","CHR:hg38POS", "p")|>
    dplyr::mutate(
      chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
      pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2])
    ) |> dplyr::rename(gwas_rsID = "rsID") |>
    dplyr::rename(R2 = "ldbuddy_R2", var_chr = "chr_ldbuddy", var_from ="pos_ldbuddy",rsID="ldbuddy_rsID",nom_pval ="p")
 
  pageCreate(width = 4, height =4 , default.units = "inches", showGuides = FALSE)
  
  GWAS_manhattan_plot(
    gwas_data = gwas_data_subset,
    gene_highlights = gene_hl,
    label= paste0("GWAS-", i),
    rsID = gwas_rsID,
    ctl_signal = ctl_signal,
    fnf_signal = fnf_signal,
    x_start = 0.7,
    y_start = 1,
    width = 2.5,
    height = 1,
    signal_height = 0.43,
    padding = 10000
  )
  
  
  
  }
  dev.off()
}






#------------------------------------------------------------------------------

GWAS_subtype_leads <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/",
                              i,
                              "/leads/EUR_",i,"_leads_ld_final.csv"))

gwas_data_subset <- GWAS_subtype_leads |> dplyr::filter(rsID  == gwas_rsID) |>
  dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
  dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "ldbuddy_rsID","rsID","CHR:hg38POS", "p")|>
  dplyr::mutate(
    chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
    pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2])
  ) |> dplyr::rename(gwas_rsID = "rsID") |>
  dplyr::rename(R2 = "ldbuddy_R2", var_chr = "chr_ldbuddy", var_from ="pos_ldbuddy",rsID="ldbuddy_rsID",nom_pval ="p")

gwas_locus_plot <- prepare_locus_plot(gwas_data_subset)
chrom <- gwas_data_subset$var_chr |> unique()
minregion <-min(gwas_data_subset$var_from) - 10000
maxregion <- max(gwas_data_subset$var_from) + 10000
x_start =0.7
y_start = 1
width = 2.5
height = 1
# Set up plot parameters
region_pg <- pgParams(assembly = "hg38", chrom = chrom,
                      chromstart = minregion,
                      chromend = maxregion,
                      x = x_start, y = y_start, width = width, height = height)

gwas_ylim <- ceiling(max(log10(gwas_locus_plot$p)*-1.1)) 

pageCreate(width = 4, height =4 , default.units = "inches", showGuides = TRUE)

GWAS_manhattan_plot(
  gwas_data = gwas_data_subset,
  gene_highlights = gene_hl,
  label= paste0("GWAS-", "AllOA"),
  rsID = gwas_rsID,
  ctl_signal = ctl_signal,
  fnf_signal = fnf_signal,
  x_start = 0.7,
  y_start = 1,
  width = 2.5,
  height = 1,
  signal_height = 0.43,
  padding = 75000
)



# Signal data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_AI/rna_output/signals/synovium/merged_signal/"
ctl_signal <- paste0(dir_merged,"CTL_sorted.bw")
fnf_signal <- paste0(dir_merged,"FNF_sorted.bw")

signal_height= 0.43



