create_correlation_plot <- function(pca_df, condition = NULL, title = "") {
  # Filter data by condition if specified
  if (!is.null(condition)) {
    pca_df <- pca_df[pca_df$Condition == condition, ]
  }
  
  # Select relevant columns for correlation
  # Exclude certain columns that shouldn't be in correlation
  cols_to_exclude <- c("Row.names", "sampleID", "Condition", "Proj", "Time","Protocol_notes",
                       "Tissue", "TreatmentDate", "sizeFactor", "Seq_Rep","FragmentBatch")
  
  # Separate PCs and predictors
  pc_cols <- grep("^PC[0-9]+$", colnames(pca_df), value = TRUE)
  pred_cols <- setdiff(colnames(pca_df), c(pc_cols, cols_to_exclude))
  
  # Select only numeric/factor columns from predictors
  pred_df <- pca_df[, pred_cols]
  num_df <- pred_df[, sapply(pred_df, function(x) is.numeric(x) )]
  
  pc_matrix <- as.matrix(pca_df[, pc_cols[1:10]])
  num_matrix <- as.matrix(num_df)
  
  cat_df <- pred_df[, sapply(pred_df, function(x) is.character(x) | is.factor(x))]
  
  # Convert categorical to numeric
  cat_df_numeric <- DataFrame(lapply(cat_df@listData, function(x) {
    if(is.factor(x) | is.character(x)) {
      as.numeric(factor(x))
    } else {
      x
    }
  }))
  cat_matrix <- as.matrix(cat_df_numeric)
  
  # Numeric variables - using Pearson
  numeric_cors <- cor(pc_matrix, 
                      num_matrix, 
                      method = "pearson", 
                      use = "pairwise.complete.obs")
  
  cat_cors <- cor(pc_matrix, 
                  cat_matrix, 
                  method = "spearman",
                  use = "pairwise.complete.obs")
  
  # Combine correlations
  correlations <- cbind(numeric_cors, cat_cors)

  # Calculate p-values for combined correlations
  # Calculate p-values matching correlation structure
  p_values <- matrix(NA, nrow=nrow(correlations), ncol=ncol(correlations))
  colnames(p_values) <- colnames(correlations)
  rownames(p_values) <- rownames(correlations)
  
  for(i in 1:nrow(correlations)) {
    for(j in 1:ncol(correlations)) {
      if(colnames(correlations)[j] %in% colnames(num_matrix)) {
        test <- cor.test(pc_matrix[,i], num_matrix[,colnames(correlations)[j]], 
                         method="pearson")
      } else {
        test <- cor.test(pc_matrix[,i], cat_matrix[,colnames(correlations)[j]], 
                         method="spearman")
      }
      p_values[i,j] <- test$p.value
    }
  }
  
  # Prepare for plotting
  melted_corr <- reshape2::melt(correlations)
  melted_p <- reshape2::melt(p_values)
  names(melted_corr) <- c("x", "y", "cor")
  names(melted_p) <- c("x", "y", "pval")
  
  all_data <- merge(melted_corr, melted_p, by = c("x", "y"))
  
  # Add significance labels
  all_data$pval_sig <- ifelse(all_data$pval < 0.001, "***",
                              ifelse(all_data$pval < 0.01, "**",
                                     ifelse(all_data$pval < 0.05, "*", "")))
  
  # Format correlation values
  all_data$cor_text <- sprintf("%.2f", all_data$cor)
  
  # Ensure correct PC ordering
  all_data$x <- factor(all_data$x, levels = rev(pc_cols))
  
  # Create plot
  plot <- ggplot(all_data, aes(x = y, y = x, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#418C82", mid = "white", high = "#C55E2D") +
    geom_text(aes(label = cor_text), family = "Helvetica", size = 2) +
    geom_text(aes(label = pval_sig), family = "Helvetica", size = 2, fontface = "bold",
              vjust = -0.75) +
    guides(fill = guide_colorbar(title = "",
                                 title.position = "top",
                                 title.hjust = 1, 
                                 direction = "vertical")) +
    theme(axis.text.x = element_text(angle = 45, size = 7, vjust = 1, hjust = 0.8,
                                     color = "black"),
          axis.text.y = element_text(size = 7, color = "black"),
          panel.background = element_rect(fill = 'transparent', color = "transparent"),
          plot.background = element_rect(fill = 'transparent', color = "transparent"),
          strip.background = element_blank(),
          text = element_text(family = "Helvetica"),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = 4, color = "black"),
          legend.text = element_text(size = 4, color = "black"),
          legend.position = "right",
          legend.key.width = unit(0.6, 'mm'),
          legend.margin = margin(t = -35),
          strip.text = element_text(size = 10),
          title = element_blank())
  
  return(plot)
}
