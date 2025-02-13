# Gene expression with the 101 donors? maybe stick to the 101 donors
# Not sure I should look into the first 15 donors?  since the AM7769-CTL sample have high CHIPMIX 
# So maybe 14 Donors instead of 16?

setwd("/work/users/s/e/seyoun/CQTL_AI")

library(DESeq2)
library(dplyr)
library(tximeta)
#-------------------------------------------------------------------------------

# Load configuration from YAML file
config <- yaml::read_yaml("config/RNAconfig.yaml")
aligned_samples <- fread(config$samplesheet)
rna_samples <- aligned_samples |> dplyr::filter(Tissue != "synovium")
donorsheet <- fread("/work/users/s/e/seyoun/CQTL_AI/Chon_DonorInfo.txt") 

#Check
donorsheet[grep("\\(ankle\\)|\\(femur\\)", donorsheet$Donor), ]
donorsheet <- donorsheet[!grepl("\\(femur\\)", Donor)]
donorsheet[, Donor := gsub(" \\(ankle\\)", "", Donor)]

#verify
donorsheet[grep("AM7778", donorsheet$Donor), ]

RNAmerged_meta <- left_join(rna_samples, donorsheet, by = "Donor") |> 
  filter(Tissue == "Ankle") |>
  filter(Donor != "AM7769") |>
  dplyr::select("Proj","Donor","Condition","Tissue","Time","Tech_Rep","Seq_Rep","RNAextractDate",
                "RIN","RNAQubit","Sex","Age","Race","RNAlibrariesseqDate",
                "OAGradeAvg","CauseOfDeath","FragmentBatch")



#There is no replicate from the sequencing

#Add the age group
RNAmerged_meta$AgeGroup <- cut(RNAmerged_meta$Age,
                                      breaks=c(24, 44, 64, 84),
                                      labels=c("25_44", "45_64", "65_84"))



RNAmerged_meta_fi <- RNAmerged_meta |> 
  mutate(sampleID = paste(Proj, Donor, Condition, Tissue,Tech_Rep, sep="_"))

RNAmerged_meta_fi$Donor <- factor(RNAmerged_meta_fi$Donor)
RNAmerged_meta_fi$Condition <- factor(RNAmerged_meta_fi$Condition, levels = c("CTL", "FNF"))
RNAmerged_meta_fi$Sex <- factor(RNAmerged_meta_fi$Sex)

## Add quant paths and names
RNAmerged_meta_fi$files <- file.path("rna_output", "quant", RNAmerged_meta_fi$sampleID, "quant.sf")
colnames(RNAmerged_meta_fi) <- gsub("sampleID", "names", colnames(RNAmerged_meta_fi))
file.exists(RNAmerged_meta_fi$files)
coldata <- RNAmerged_meta_fi

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)

#Gene level---------------------------------------------------------------------
## Convert to factors (avoids a warning)
colData(gse)[] <- lapply(colData(gse), factor)

## Build DESeq object
dds_gene <- DESeqDataSet(gse, design = ~ Condition+Donor)

## Filter out lowly expressed genes
## (at least 10 counts in at least 2 samples)
keep <- rowSums(counts(dds_gene) >= 10) >= ceiling(nrow(colData(gse))*0.10)
dds_gene <- dds_gene[keep,]

## Fit model
dds_gene <- DESeq(dds_gene)


## Gene results
rna_res <- results(dds_gene)
print(summary(rna_res,alpha=0.05, na.rm=TRUE))

# Shrink l2fc
de_RNAgenes_shrink <- lfcShrink(dds_gene,
                             coef = "Condition_FNF_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

# Join results with gene info
de_RNAgenes_shrink <-
  inner_join(x = as.data.frame(de_RNAgenes_shrink),
             y = as.data.frame(rowData(gse)) %>%
               dplyr::select(c("gene_id", "tx_ids")),
             by = "gene_id") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  as.data.frame()

save(dds_gene,rna_res,de_RNAgenes_shrink,
     file="/work/users/s/e/seyoun/CQTL_AI/rna_output/Diff_analysis_chon/RNA_deseq.RData")
