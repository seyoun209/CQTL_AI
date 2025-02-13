library(csaw)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(cqn)
library(DESeq2)
library(edgeR)
library(dplyr)
library(data.table)

#-------------------------------------------------------------------------------
# peak calling- fragment length

cqtl_bam_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/wasp/rmdup/"
cqtl_bam_files <- dir(cqtl_bam_dir, pattern="Ankle_1.sorted_wasp.bam$", full.names=TRUE) 

output_dir <- "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only"
dir.create(output_dir, showWarnings =TRUE)

# Start timing
message("Analysis started at: ", Sys.time())
ptm <- proc.time()


# Create PDF for fragment length histograms
pdf(file.path(output_dir, 'CQTL_fragment_lengths_histogram.pdf'))

# Initialize vectors to store results
frag_lens <- c()
sample_names <- c()
mapped_reads <- c()
frag_sizes_li <- NULL


# Process each BAM file
for (i in 1:length(cqtl_bam_files)) {
  current_file <- cqtl_bam_files[i]
  sample_name <- gsub(".sorted_wasp.bam","",basename(current_file))
  
  message("Processing file ", i, " of ", length(cqtl_bam_files), ": ", sample_name)
  
  # Get fragment sizes
  out <- getPESizes(current_file)
  
  # Extract and store the number of mapped reads
  num_mapped <- out$diagnostics["mapped.reads"]  # Extract mapped reads
  mapped_reads <- c(mapped_reads, num_mapped)

  # Store median fragment length
  frag_lens <- c(frag_lens, median(out$sizes))
  sample_names <- c(sample_names, sample_name)
  frag_sizes_li[[i]] <- out$sizes
  names(frag_sizes_li)[i] <- sample_name
  
  # Extract sample info from filename
  sample_info <- strsplit(sample_name, "_")[[1]]
  sample_id <- sample_info[2]
  condition <- sample_info[3]  # CTL or FNF
  
  # Create histogram
  frag_sizes <- out$sizes
  hist(frag_sizes, 
       breaks = 50, 
       xlab = "Fragment sizes (bp)",
       ylab = "Frequency", 
       main = paste0(sample_id, " - ", condition), 
       col = if(condition == "CTL") "lightblue" else "lightcoral")
  abline(v = 1500, col = "red")
  abline(v = median(out$sizes), col = "blue", lty = 2)
  legend("topright", 
         legend = c("1500bp cutoff", "Median length"), 
         col = c("red", "blue"), 
         lty = c(1, 2))
  
}
dev.off()

# Create results dataframe with sample information
results <- data.frame(
  Sample = sample_names,
  SampleID = sapply(strsplit(sample_names, "_"), `[`, 2),
  Condition = sapply(strsplit(sample_names, "_"), `[`, 3),
  Tissue = sapply(strsplit(sample_names, "_"), `[`, 4),
  MedianFragmentLength = frag_lens
)

write.csv(results, 
          file.path(output_dir, "CQTL_fragment_analysis_summary.csv"), 
          row.names = FALSE)



# Save R data
save(sample_names, frag_lens, results,mapped_reads, frag_sizes_li, 
     file = file.path(output_dir, "CQTL_fragment_analysis.RData"))

# Print execution time
print(proc.time() - ptm)
message("Analysis completed at: ", Sys.time())




#------------------------------------------------------------------------------
#peak calling-counts

message('Loading fragment length data')
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_fragment_analysis.RData")

# Set parameters for counting
pe.param = readParam(max.frag=1000,    # Maximum fragment length
                     pe="both",         # Use both reads of pair
                     minq=20)           # Minimum mapping quality
multi.frag.len = list(frag_lens, NA)

# Set window parameters
win.width = 10                        # Binding site length
spacing = win.width                   # Space between windows


mean_frag_len <- mean(frag_lens)
median_frag_len <- median(frag_lens)
message("Mean fragment length: ", round(mean_frag_len, 2))
message("Median fragment length: ", round(median_frag_len, 2))


#spacing should not be larger than ext/2 for analyses with small windows. If ext is also very small, spacing should be set to width to avoid loading too many small windows.
message('Starting window counts')
# Generate window counts
data = windowCounts(cqtl_bam_files, 
                    ext=multi.frag.len, 
                    filter=(5*length(cqtl_bam_files)), # Filter low-abundance windows
                    spacing=spacing, 
                    param=pe.param, 
                    width=win.width)
# Merge nearby windows
message('Merging windows')
merged <- mergeWindows(rowRanges(data), tol=100L)
my.regions = merged$region

save(data, file="/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_CSAW_windows_count")
save(my.regions, 
     file="/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_regions.RData")

reg.counts <- regionCounts(cqtl_bam_files, my.regions, ext=multi.frag.len, param=pe.param)
counts=assay(reg.counts)
save(counts,my.regions, results, file=paste0("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_CSAW_windowCounts.Rdata"))
#------------------------------------------------------------------------------
#peak cqn normalization


# Load the window counts and regions
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_fragment_analysis.RData")
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_CSAW_windows_count")
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_regions.RData")
load("/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_CSAW_windowCounts.Rdata")

# flagstat_data <- fread("/work/users/s/e/seyoun/CQTL_AI/output/wasp/rmdup/qc_cqtl/multiqc_data/multiqc_samtools_flagstat.txt", 
#                        header = TRUE, 
#                        stringsAsFactors = FALSE)
# mapped_reads <- data.frame(
#   File = gsub("_stats$", "", flagstat_data$Sample),  # Remove "_stats" suffix
#   MappedReads = flagstat_data$mapped_passed          # Use mapped_passed column
# )


# write.table(mapped_reads,
#             file = "/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/mapped_read_counts.txt",
#             sep = "\t",
#             quote = FALSE,
#             row.names = FALSE)
# 
# mapped_reads <- mapped_reads[match(sample_names, mapped_reads$File), ]

mappedsReads_df <- cbind(sample_names, as.data.frame(mapped_reads))
colnames(mappedsReads_df) <- c("sample_names","MappedReads")
#------------------------------------------------------------------------------
# CQN Offset
# Get GC content for peaks

windowsGR = my.regions
peakwidths = width(windowsGR)
message('calculating GC content')
peakseqs = getSeq(Hsapiens, 
                  seqnames(windowsGR),
                  start(windowsGR),
                  end(windowsGR))
GCcontent = sapply(peakseqs, 
                   function(x) (length(gregexpr("[GC]",x)[[1]]))/nchar(x))

# Save GC content
save(GCcontent, sample_names, 
     file = '/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_GCContent.RData')


#CQN----------------------------------------------------------------------------

# Run CQN normalization
message('running cqn')
cqn.counts.sqn = cqn(counts,
                     lengths = peakwidths,
                     x = GCcontent,
                     sizeFactors = mappedsReads_df$MappedReads,
                     sqn = TRUE)



# Calculate normalization factors
cqnOffset = cqn.counts.sqn$glm.offset
cqnNormFactors = exp(cqnOffset) #exponential function to the offsets to get the raw normalization factors.
cqnNormFactors = cqnNormFactors / exp(rowMeans(log(cqnNormFactors))) #normalization factors around 1.

# Save normalization factors
save(cqnNormFactors, results,
     file = '/work/users/s/e/seyoun/CQTL_AI/output/peaks/csaw/ankle_only/CQTL_cqnNormFactors.RData')
