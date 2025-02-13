#ATACqc_another
## Started: Nov 27th 2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_AI")
library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(Rsamtools)

args <- commandArgs(TRUE)

bamfiles_li <- args[1]  # Input BAM path
output_dir <- args[2]
name_subset <- args[3]

print(bamfiles_li)
print(output_dir)
print(name_subset)

# Verify input file exists
if (!file.exists(bamfiles_li)) {
  stop(paste("Input BAM file does not exist:", bamfiles_li))
}

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create label from filename
bamFileLabels <- gsub(name_subset, "", basename(bamfiles_li))

print(paste("Processing BAM file:", bamfiles_li))
print(paste("Output label:", bamFileLabels))

bamqc <- bamQC(bamfiles_li, outPath = NULL)
save(bamqc, file=paste0(output_dir,bamFileLabels,"_bamQC"))
print("bamqc qc done")




# TSS NEIFHTMENT
# tss <- read.table("/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.tss.bed.gz", 
#                   col.names=c("chr", "start", "end", "name", "gene", "strand"))
# tss <- makeGRangesFromDataFrame(tss, keep.extra.columns=TRUE)
# possibleTag <- combn(LETTERS, 2)
# possibleTag <- c(paste0(possibleTag[1,], possibleTag[2,]),
#                  paste0(possibleTag[2,], possibleTag[1,]))
# 
# bamTop100 <- scanBam(BamFile(bamfiles_li, yieldSize=100),
#                      param=ScanBamParam(tag=possibleTag))[[1]]$tag
# tags <- names(bamTop100)[lengths(bamTop100)>0]
# ojs <- splitGAlignmentsByCut(gal, txs=tss, genome=Hsapiens)
# gal <- readBamFile(bamfiles_li, tag=tags)
# save(gal, file= paste0)
# 
# # Calculate library size
# librarySize <- estLibSize(bamfiles_li)
# 
# # Calculate enrichment
# TSSE <- TSSEscore(gal,tss,librarySize=librarySize,upstream=1000,downstream=1000)
# 
# # Plot enrichment profile
# plotTSSEnrichment(TSSE, tss)
