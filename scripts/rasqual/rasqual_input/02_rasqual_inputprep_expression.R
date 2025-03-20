#rasqual count and the offset
setwd("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs")
library(rasqualTools)
library(dplyr)
library(stringr)

#Expression data is ready

#Adding placeholder and make the condition specific

#CQN Normalized factors (cqnNormFactors, mapped_reads_ordered, counts_macs2, macs2_gr,macs2_filtered)
load("/work/users/s/e/seyoun/CQTL_AI/output/Differential_analysis_CQTL/CQTL_cqnNormFactors.RData")

# Create ID format
macs2_filtered$ID <- paste(macs2_filtered$Chr, macs2_filtered$Start, macs2_filtered$End, sep="_")

# Identify CTL and FNF columns
ctl_cols <- grep("CTL", colnames(macs2_filtered), value=TRUE)
fnf_cols <- grep("FNF", colnames(macs2_filtered), value=TRUE)

# Create CTL file
ctl_count <- macs2_filtered %>%
  select(ID, Chr, Start, End, all_of(ctl_cols))

# Create FNF file
fnf_count <- macs2_filtered %>%
  select(ID, Chr, Start, End, all_of(fnf_cols))

# Write the two output files
write.table(ctl_count, "peaks_CTL.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(fnf_count, "peaks_FNF.txt", row.names=FALSE, quote=FALSE, sep="\t")




#Create the count data file for Rasqual

Chromosomes <- paste0("chr",c(1:22, "X"))

# For CTL counts
for (i in Chromosomes){
  Chr_specific <- subset(ctl_count, ctl_count$Chr == i)
  counts_data <- Chr_specific[,5:ncol(Chr_specific)]
  row.names(counts_data) <- Chr_specific[,1]
  counts_data[] <- lapply(counts_data, function(x) as.numeric(x))
  ID <- paste0(i, "_CTL") 
  countslist <- list(counts_data)
  names(countslist) <- ID
  saveRasqualMatrices(countslist, "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs", file_suffix = "expression")
}
countslist <- NULL
# For FNF counts

for (i in Chromosomes){
  Chr_specific <- subset(fnf_count, fnf_count$Chr == i)
  counts_data <- Chr_specific[,5:ncol(Chr_specific)]
  row.names(counts_data) <- Chr_specific[,1]
  counts_data[] <- lapply(counts_data, function(x) as.numeric(x))
  ID <- paste0(i, "_FNF")
  countslist <- list(counts_data)
  names(countslist) <- ID
  saveRasqualMatrices(countslist, "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs", file_suffix = "expression")
}


#create size factor file--------------------------------------------------------

# Split for CTL
ctl_samples <- grep("CTL", colnames(cqnNormFactors), value=TRUE)
CQN_CTL <- cqnNormFactors[,ctl_samples]

# Split for FNF
fnf_samples <- grep("FNF", colnames(cqnNormFactors), value=TRUE)
CQN_FNF <- cqnNormFactors[,fnf_samples]

# Process CTL normalization factors
CQN_CTL <- as.data.frame(cbind(CQN_CTL, ctl_count$Chr, ctl_count$ID))
Colcount_ctl <- ncol(CQN_CTL)
colnames(CQN_CTL)[Colcount_ctl-1] <- "Chr"
colnames(CQN_CTL)[Colcount_ctl] <- "PeakID"

# Process FNF normalization factors
CQN_FNF <- as.data.frame(cbind(CQN_FNF, fnf_count$Chr, fnf_count$ID))
Colcount_fnf <- ncol(CQN_FNF)
colnames(CQN_FNF)[Colcount_fnf-1] <- "Chr"
colnames(CQN_FNF)[Colcount_fnf] <- "PeakID"

# Process by chromosome for CTL
Chromosomes <- paste0("chr",c(1:22, "X"))
for (i in Chromosomes){
  Chr_specific <- subset(CQN_CTL, CQN_CTL$Chr == i)
  CQNstuff <- ncol(CQN_CTL)-2
  normfactors <- Chr_specific[,1:CQNstuff]
  row.names(normfactors) <- Chr_specific$PeakID
  normfactors[] <- lapply(normfactors, function(x) as.numeric(x))
  ID <- paste0(i, "_CTL")
  normlist <- list(normfactors)
  names(normlist) <- ID
  saveRasqualMatrices(normlist, "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs", file_suffix = "size_factors")
}

# Process by chromosome for FNF
for (i in Chromosomes){
  Chr_specific <- subset(CQN_FNF, CQN_FNF$Chr == i)
  CQNstuff <- ncol(CQN_FNF)-2
  normfactors <- Chr_specific[,1:CQNstuff]
  row.names(normfactors) <- Chr_specific$PeakID
  normfactors[] <- lapply(normfactors, function(x) as.numeric(x))
  ID <- paste0(i, "_FNF")
  normlist <- list(normfactors)
  names(normlist) <- ID
  saveRasqualMatrices(normlist, "/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs", file_suffix = "size_factors")
}

#Create peak info file----------------------------------------------------------
Chromosomes <- paste0("chr",c(1:22, "X"))

for (i in Chromosomes){
  Chr <- subset(macs2_filtered, macs2_filtered$Chr == i)
  
  peak_file <- data.frame(
    PeakID = paste(Chr$Chr, Chr$Start, Chr$End, sep="_"),
    Chr = Chr$Chr,
    Start = Chr$Start,
    End = Chr$End
  )
  file.name <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/inputs/peak_info_", i, ".txt")
  write.table(peak_file, file = file.name, quote = FALSE, row.names = FALSE, col.names = TRUE)
}

#Create peak info file_This is with GeneID ---------------------------------------------------------

Chromosomes <- paste0("chr",c(1:22, "X"))

for (i in Chromosomes){
  Chr <- subset(macs2_filtered, macs2_filtered$Chr == i)

  peak_file <- data.frame(
    GeneID =Chr$Geneid,
    PeakID = paste(Chr$Chr, Chr$Start, Chr$End, sep="_"),
    Chr = Chr$Chr,
    Start = Chr$Start,
    End = Chr$End
  )
  file.name <- paste0("/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/peak_info/peak_info_", i, ".txt")
  write.table(peak_file, file = file.name, quote = FALSE, row.names = FALSE, col.names = TRUE)
}


