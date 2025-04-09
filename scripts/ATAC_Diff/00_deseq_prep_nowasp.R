#This is updated to use the no-wasp
setwd("/work/users/s/e/seyoun/CQTL_AI/output/")
# Differential atac-seq -Chondrocyte
#-------------------------------------------------------------------------------
#library
library(DESeq2)
library(Rsubread,lib.loc = "/nas/longleaf/home/seyoun/R/x86_64-pc-linux-gnu-library/4.2")
library(GenomicRanges)
library(base)
library(data.table)
library(tidyverse)
library(httpgd)
#--------------------------------------------------------------------------------
#Functions 
# Read and clean data
read_and_clean_peaks <- function(file) {
  df <- read.table(file, header=TRUE, skip=1, sep="\t")
  colnames(df) <- gsub("output.filtered.blk_filter.", "", colnames(df))
  colnames(df) <- gsub(".sorted_final.bam", "", colnames(df))
  return(df)
}
#--------------------------------------------------------------------------------
#Read in the data macs2

macs2 <- read_and_clean_peaks("peaks/Ankle/merged/allsamples_macs2_merged_counts.txt")

#Make the meta for samples
# sample config file generate (Metasamplesheet)
config <- yaml::read_yaml("../config/ATACconfig.yaml")
samplesheet <- fread(paste0("/work/users/s/e/seyoun/CQTL_AI/",config$samplesheet))
donorsheet <- fread("/work/users/s/e/seyoun/CQTL_AI/Chon_DonorInfo.txt") 

#Check
donorsheet[grep("\\(ankle\\)|\\(femur\\)", donorsheet$Donor), ]
donorsheet <- donorsheet[!grepl("\\(femur\\)", Donor)]
donorsheet[, Donor := gsub(" \\(ankle\\)", "", Donor)]

#verify
donorsheet[grep("AM7778", donorsheet$Donor), ]

merged_meta <- left_join(samplesheet, donorsheet, by = "Donor") |> 
  filter(Tissue == "Ankle") |>
  dplyr::select("Proj","Donor","Condition","Tissue","Protocol_notes","Time","Tech_Rep","Seq_Rep","TreatmentDate",
                "ATACProtocolDate","ATACLibrarySubmissionDate","Sex","Age","Race",
                "OAGradeAvg","CauseOfDeath","FragmentBatch")



#Take the replicate out
pairs_count <- merged_meta[, .N, by=.(Donor, Condition)]
duplicated_pairs <- pairs_count[N > 1, .(Donor, Condition)]

meta_data_rep_cleaned <- merged_meta[!(Donor %in% unique(duplicated_pairs$Donor)) | 
                                       (Donor %in% unique(duplicated_pairs$Donor) & Tech_Rep == 1)]

#Add the age group
meta_data_rep_cleaned$AgeGroup <- cut(meta_data_rep_cleaned$Age,
                                      breaks=c(24, 44, 64, 84),
                                      labels=c("25-44", "45-64", "65-84"))



meta_final <- meta_data_rep_cleaned |> 
  mutate(sampleID = paste(Proj, Donor, Condition, Tissue, Protocol_notes, sep="_"))

meta_final$Donor <- factor(meta_final$Donor)
meta_final$Condition <- factor(meta_final$Condition, levels = c("CTL", "FNF"))

#Change to macs2_gr
macs2_gr <- makeGRangesFromDataFrame(macs2)
macs2_gr$peakID <- macs2$Geneid
counts_macs2 <- as.matrix(macs2[,7:ncol(macs2)])
se <- SummarizedExperiment(assays=list(counts = counts_macs2),
                               rowRanges=macs2_gr, colData=meta_final)

save(se, macs2_gr, counts_macs2, meta_final,
file = "/work/users/s/e/seyoun/CQTL_AI/output/diff_deseq2/condition/QC/Chon_macs2_se.RData")
