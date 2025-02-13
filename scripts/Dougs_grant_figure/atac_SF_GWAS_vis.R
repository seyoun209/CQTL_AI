# Signal plot 
#functions
library(plotgardener)
library(grid)
library(RColorBrewer)
library(biomaRt)

# Signal data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_AI/output/signals/synovium/merged_signal/"
atac_ctl_signal <- paste0(dir_merged,"CTL_sorted.bw")
atac_fnf_signal <- paste0(dir_merged,"FNF_sorted.bw")

load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/synovium/synovium_deseq.Rdata") # atacDDS,res,res_Shrink 

diff_atac_sig_gr <- res_Shrink %>%
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) 

diffATAC_anno <- annotatePeak(diff_atac_sig_gr, TxDb=txdb,tssRegion=c(-1000, 500), flankDistance= 0)
diffATAC_anno_df <-as.data.frame(diffATAC_anno)

#diffATAC_promoter_df  <- diffATAC_anno_df[grepl("Promoter",diffATAC_anno_df$annotation ),]
ensembl_version = "https://dec2016.archive.ensembl.org"
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host=ensembl_version)
diffATAC_promoter_df_noV <-diffATAC_anno_df |> dplyr::mutate(geneId = gsub("\\..*$", "", geneId))

annots <- AnnotationDbi::select(org.Hs.eg.db, diffATAC_promoter_df_noV$geneId,
                                columns=c("SYMBOL"), keytype="ENSEMBL")
annots <- dplyr::rename(annots, geneId = ENSEMBL)
annots_unique <- annots %>% distinct(geneId, .keep_all = TRUE)
diffATAC_promoter_df_hgnc <- left_join(diffATAC_promoter_df_noV, annots_unique, by = "geneId")


#Make it grange
diffATAC_anno_gr <- GRanges(
  seqnames = diffATAC_promoter_df_hgnc$seqnames,
  ranges = IRanges(start = diffATAC_anno_df$start, end = diffATAC_promoter_df_hgnc$end),
  strand = diffATAC_promoter_df_hgnc$strand,
  mcols = diffATAC_promoter_df_hgnc[, !colnames(diffATAC_promoter_df_hgnc) %in% c("seqnames", "start", "end", "strand")]
)

#diffATAC_anno_gr <- diffATAC_anno_gr |> filter(mcols.SYMBOL %in%  rna_promotergenes_filtered_hgnc$SYMBOL)
#Make it grange

#diff_atac_sig_gr <- diff_atac_sig %>%
#  mutate(start = ifelse(strand == "+", start - 100, start),
#         end = ifelse(strand == "-", end + 100, end))



OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA",
                "THR", "ThumbOA", "TJR", "TKR")
for (i in OAsubtypes) {
  print(i)
  pdf(paste0("output/Differential_analysis/plot/sigATAC_signal_plot_",i,".pdf"),width=4,height=4)
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
  LDs_gr <- with(LDs[LDs$R2 > 0.6,], GRanges(LD_chr, IRanges(LD_pos, width=1)))
  
  overlaps <- findOverlaps(diff_atac_sig_gr, LDs_gr)
  if(length(overlaps) == 0){next}
  first_hits <- data.frame(
    queryHits = queryHits(overlaps),
    subjectHits = subjectHits(overlaps)
  ) %>%
    group_by(queryHits) %>%
    slice(1) %>%
    ungroup()
  overlapping_genes_promoter <- diff_atac_sig_gr[queryHits(overlaps) |> unique()]
  
  for (j in 1:length(first_hits$queryHits)){
    print(j)
    overlapping_genes_gwasLD <- LDs_gr[first_hits$subjectHits[j]]
    matching_GWASrows <- subtype_leads[subtype_leads$`ldbuddy_CHR:hg38POS` %in% gsub("chr","",as.character(overlapping_genes_gwasLD)),]
    overlapping_genes_promoter_Row <- diff_atac_sig_gr[first_hits$queryHits[j]]
    gwas_rsID <- unique(matching_GWASrows$gwas_rsID)
    gene_hl <- overlapping_genes_promoter_Row$mcols.SYMBOL
    
    start_atac <- start(overlapping_genes_promoter_Row)
    end_atac <- end(overlapping_genes_promoter_Row)
    
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
    
    GWAS_manhattan_plot_atac(
      gwas_data = gwas_data_subset,
      gene_highlights = gene_hl,
      label= paste0("GWAS-", i),
      start_anno=start_atac,
      end_anno=end_atac,
      rsID = gwas_rsID,
      atac_ctl_signal = atac_ctl_signal,
      atac_fnf_signal = atac_fnf_signal,
      ctl_signal = ctl_signal,
      fnf_signal = fnf_signal,
      x_start = 0.7,
      y_start = 1,
      width = 2.5,
      height = 1,
      signal_height = 0.43,
      padding = -10000
    )
    
    
    
  }
  dev.off()
}

