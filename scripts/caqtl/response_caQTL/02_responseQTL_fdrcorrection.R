library(dplyr)
library(data.table)
library(stringr)

# Load the data
setwd("/work/users/s/e/seyoun/CQTL_AI/output/")
#response_qtl_results,
load("rasqual_output/response_qtl/caQTL_response_qtl_results.RData")
head(response_qtl_results)

#--------------------------------------------------------
# FDR correction for response QTL results

response_qtl_results_order <- response_qtl_results[order(response_qtl_results$summary_pvalue),] 

# Fixing class 

response_qtl_results_order$summary_pvalue <- as.numeric(response_qtl_results_order$summary_pvalue)

#response_qtl_results_order$annov_fdr <- p.adjust(response_qtl_results_order$anova_interaction_pval, 
#method = "fdr", n=nrow(response_qtl_results_order))

response_qtl_results_order$lmer_fdr <- p.adjust(response_qtl_results_order$summary_pvalue, 
method = "fdr", n=nrow(response_qtl_results_order))

#sig.snp <- response_qtl_results_order |> filter(annov_fdr < 0.05)
sig_response_qtl <- response_qtl_results_order |> filter(lmer_fdr < 0.05) |> 
arrange(lmer_fdr)


#--------------------------------------------------------
# effect size and filter the rasqual data
#--------------------------------------------------------

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
mutate(Condition = "pbs",
         Peak_SNP = paste(peak, snp, sep = "_"))
fnf_rasqual_sig <- inner_join(fnf_qtl_sig, fnf_rasqual, by = c("peak" = "Feature", "snp" = "rs_ID")) |>
mutate(Condition = "fnf",
         Peak_SNP = paste(peak, snp, sep = "_"))

pbs_rasqual_update <- pbs_rasqual |> mutate(Condition = "pbs",
         Peak_SNP = paste(Feature, rs_ID, sep = "_"))
fnf_rasqual_update <- fnf_rasqual |> mutate(Condition = "fnf",
         Peak_SNP = paste(Feature, rs_ID, sep = "_"))

response_caQTL_final <- NULL
for (i in 1:nrow(sig_response_qtl)){
  print(i)
  Sample <- sig_response_qtl[i,]
  
  pbs <- pbs_rasqual_sig |> filter(Peak_SNP == Sample$Peak_SNP)
  fnf <- fnf_rasqual_sig |> filter(Peak_SNP == Sample$Peak_SNP)
  pbs_bulk <- pbs_rasqual_update |> filter(Peak_SNP == Sample$Peak_SNP)
  fnf_bulk <- fnf_rasqual_update |> filter(Peak_SNP == Sample$Peak_SNP)
  
  if (nrow(pbs) > 0){
    
    Sample$Ref <- pbs$Ref_Allele
    Sample$Alt <- pbs$Alt_allele 
    Sample$RASQUAL_PValue_pbs <- pbs$PValue
    Sample$Pi_pbs <- pbs$Effect_Size
    Sample$EffectSize_pbs <- Sample$Pi_pbs-0.5
    # adding FNF
    Sample$RASQUAL_PValue_fnf <- fnf_bulk$PValue
    Sample$Pi_fnf <- fnf_bulk$Effect_Size
    Sample$EffectSize_fnf <- Sample$Pi_fnf-0.5


  } else {
    
    Sample$Ref <- fnf$Ref_Allele
    Sample$Alt <- fnf$Alt_allele 
    Sample$RASQUAL_PValue_fnf <- fnf$PValue
    Sample$Pi_fnf <- fnf$Effect_Size
    Sample$EffectSize_fnf <- Sample$Pi_fnf-0.5
    #adding pbs
    Sample$RASQUAL_PValue_pbs <- pbs_bulk$PValue
    Sample$Pi_pbs <- pbs_bulk$Effect_Size
    Sample$EffectSize_pbs <- Sample$Pi_pbs-0.5

    
  }
  

  response_caQTL_final <- rbind(response_caQTL_final, Sample)
  
}


#--------------------------------------------------------
# Save the final response caQTL results

save(response_caQTL_final, sig_response_qtl, 
     file = "rasqual_output/response_qtl/response_caQTL_final.RData")




