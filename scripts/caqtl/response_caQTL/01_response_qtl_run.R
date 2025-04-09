# Load the packages
library(data.table)
library(dplyr)
library(stringr)
library(lme4)
library(lmerTest)
library(tidyverse)

#Load the data
setwd("/work/users/s/e/seyoun/CQTL_AI/output/")
#counts_vsd_region_final, all_rasqual_qtl, all_rasqual_qtl_duplicated,combined_geno_matrix, full_covar,meta_enhanced
load("rasqual_output/response_qtl/caQTL_prepdata.RData")

#------------------------------------------------
#interaction testing
#-------------------------------------------------
analyze_response_qtl <- function(qtl_results, norm_counts, geno.df, meta_df, covariates_df) {

  for (i in 1:nrow(qtl_results)) {
    print(i)
    FeatureID <- qtl_results$peak[i]

    norm_counts_matrix <- norm_counts %>%
        filter(peak %in% FeatureID) %>%
        column_to_rownames(var = "peak") %>%
        dplyr::select(-seqnames, -start, -end, -peakID)

    norm_counts_data <- data.frame(sampleID = names(norm_counts_matrix),
                           vsd = as.double(norm_counts_matrix))

    variantID <- qtl_results$snp[i]
    
    geno_data <- geno.df %>%
      filter(SNP == variantID) %>%
    column_to_rownames(var = "SNP") %>%
    dplyr::select(-c("CHR", "(C)M", "POS", "COUNTED", "ALT")) 
 

     geno_data_filtered <- data.frame(sampleID = names(geno_data),
                           genotype = as.double(geno_data))

    meta_df_filtered <- meta_df |> dplyr::select(sampleID, Donor, Condition) 
    meta_data_fi <- geno_data_filtered |>
    left_join(meta_df_filtered, by=c("sampleID"= "sampleID")) |>
      left_join(norm_counts_data,by = c("sampleID"="sampleID")) |>
      left_join(covariates_df, by=c("sampleID"= "sampleID")) 
        

    summarize_count <- meta_data_fi |>
      group_by(genotype) |>
      summarise(n = n(), .groups = "drop")

    meta_combined_all <- meta_data_fi %>%
      mutate(
        Donor = as.factor(Donor),
        Condition = ifelse(Condition == "CTL", 0, ifelse(Condition == "FNF", 1, NA))
      )

    # Define Formulas
    pc_columns <- paste(colnames(meta_combined_all)[6:length(colnames(meta_combined_all))], collapse=" + ")
    reduced_equation <- as.formula(paste("vsd ~ genotype + Condition + (1|Donor)", pc_columns, sep = "+"))
    full_equation <- as.formula(paste("vsd ~ genotype + Condition + Condition:genotype + (1|Donor)", pc_columns, sep = "+"))


    # Compare two models
    reduced_model <- lme4::lmer(reduced_equation, meta_combined_all, REML=FALSE)
    full_model <- lme4::lmer(full_equation, meta_combined_all, REML=FALSE)
    interaction_model <- lmerTest::lmer(full_equation, meta_combined_all, REML=FALSE)
    result <- anova(reduced_model, full_model)
    qtl_results$anova_interaction_pval[i] <-  result$`Pr(>Chisq)`[2]
    qtl_results$lrt_chisq[i] <- result$Chisq[2]  # Chi-square statistic
    qtl_results$lrt_df[i] <- result$Df[2]  # Degrees of freedom
    qtl_results$minor_alle_count[i] <-   summarize_count$n[3]/2
    interaction_sum <- summary(interaction_model)
    qtl_results$summary_pvalue[i] <- coef(summary(interaction_sum))["genotype:Condition", "Pr(>|t|)"] #Pvalue
    qtl_results$summary_beta[i] <- coef(summary(interaction_sum))["genotype:Condition", "Estimate"] #Beta

  }
  return(qtl_results)
}

response_qtl_results <- analyze_response_qtl(qtl_results=all_rasqual_qtl,
                                             norm_counts=counts_vsd_region_final,
                                             geno.df=combined_geno_matrix,
                                             meta_df=meta_enhanced,
                                             covariates_df=full_covar)
#------------------------------------------------

# Save the results
save(response_qtl_results,
     file = "rasqual_output/response_qtl/caQTL_response_qtl_results.RData")