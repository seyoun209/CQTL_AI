#Insert size plot
#Author Jess Byun
#date_start: December 9th 2024

setwd("/work/users/s/e/seyoun/CQTL_AI/")
library(ggplot2)
library(tidyverse)
library(scales)
library(RColorBrewer)

#-------------------------------------------------------------------------------

# Get all metric files
files <- list.files("./output/metrics",pattern = "*.txt",full.names = TRUE)

# Function to extract metrics from each file
read_metrics <- function(file) {
  data <- read.delim(file, skip = 10, header = TRUE)
  # Extract metadata from the filename
  parts <- strsplit(basename(file), "_")[[1]]
  data$proj <- parts[1]
  data$donor <- parts[2]
  data$condition <- parts[3] 
  data$tissue <- parts[4]   # Assuming all are "Ankle" in this example
  data$replicate <- parts[5]
  return(data)
}

# Read all files and combine into a single data frame
metrics_data <- files %>%
  map_dfr(read_metrics)


# Calculate mean All_Reads.fr_count for each group
metrics_mean <- metrics_data %>%
  group_by(condition, tissue, replicate, insert_size) %>%
  summarise(Mean_Reads = mean(All_Reads.fr_count), .groups = "drop")

metrics_mean$Condition <- factor(metrics_mean$condition, levels = c("CTL", "FNF"))

#individual 
# Get the unique donors in each condition and tissue
metrics_data$condition <- ifelse(metrics_data$condition == "CTL", "PBS", metrics_data$condition)

# Get the unique donors for each condition and tissue
pbs_donors <- unique(metrics_data$donor[metrics_data$condition == "PBS" & metrics_data$tissue == "synovium"])
fnf_donors <- unique(metrics_data$donor[metrics_data$condition == "FNF" & metrics_data$tissue == "synovium"])

# Number of unique donors for each condition
num_pbs_samples <- length(pbs_donors)
num_fnf_samples <- length(fnf_donors)

# Choose color palettes for both conditions (adjust the palette size accordingly)
pbs_colors <- brewer.pal(min(num_pbs_samples, 9), "OrRd")  # Max of 9 colors for YlOrRd
fnf_colors <- brewer.pal(min(num_fnf_samples, 9), "Blues")  # Max of 9 colors for Blues

# Combine the colors for both conditions
color_palette <- c(
  setNames(pbs_colors, paste("PBS", pbs_donors, sep = ".")),
  setNames(fnf_colors, paste("FNF", fnf_donors, sep = "."))
)
ggplot(metrics_data %>% filter(tissue == "synovium"),  # Filter for synovium in metrics_data
       aes(x = insert_size, y = All_Reads.fr_count, color = interaction(condition, donor))) +
  geom_line(aes(group = interaction(condition, donor, replicate)), size = 0.5)+
  scale_color_manual(values = color_palette) +  # Use the custom color palette for both conditions
  theme_classic() +
  labs(
    title = "Insert Size vs. Read Counts for Synovium Tissue",
    x = "Insert Size (bp)",
    y = "Frequency Read Counts",
    color = "Sample ID"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#Group

ggplot(metrics_mean %>% filter(tissue == "synovium"),  # Filter for synovium
       aes(x = insert_size, y = Mean_Reads, color = condition)) +
  geom_line(size = 1, aes(group = condition)) +  # Line for mean values
  scale_color_manual(values = c("CTL" = "#2057A7", "FNF" = "#F2BC40")) +  # Custom colors
  theme_classic() +
  labs(
    title = "Synovium Tissue",
    x = "Insert Size (bp)",
    y = "Mean Frequency Read Counts",
    color = "Condition"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )


#------------------------------------------------------------------------------
# Finding the insert_size metris for the mono

metrics_class <- function(file) {
  # Skip header lines and read the metrics
  data <- read.delim(file, skip=6, nrows=1)
  data$sample <- gsub("_insert_size_metrics.txt", "", basename(file))
  return(data)
}

# Read all files into a list and combine
metrics_class_calc <- do.call(rbind, lapply(files, metrics_class))

# Filter samples to exclude femur, replicate, and synovium
ankle_metrics <- metrics_class_calc %>%
  filter(!grepl("femur|replicate|synovium", tolower(sample)))


# Summary statistics
ankle_summary_stats <- ankle_metrics %>%
  summarise(
    mean_median = mean(MEDIAN_INSERT_SIZE),
    sd_median = sd(MEDIAN_INSERT_SIZE),
    mean_mode = mean(MODE_INSERT_SIZE),
    sd_mode = sd(MODE_INSERT_SIZE)
  )

print(ankle_summary_stats)


synovium_metrics <- metrics_class_calc %>%
  filter(grepl("synovium", tolower(sample)))


# Summary statistics
synovium_summary_stats <- synovium_metrics %>%
  summarise(
    mean_median = mean(MEDIAN_INSERT_SIZE),
    sd_median = sd(MEDIAN_INSERT_SIZE),
    mean_mode = mean(MODE_INSERT_SIZE),
    sd_mode = sd(MODE_INSERT_SIZE)
  )

print(synovium_summary_stats)


