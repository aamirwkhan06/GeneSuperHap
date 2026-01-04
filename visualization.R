#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap: Visualization Module (Scientifically Validated)
# Author: Aamir Khan
# Description: Publication-quality plots for all 5 unique features
# Methods based on: Breiman (2001), Smith & Cullis (2018), and standard practices
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(pheatmap)
  library(igraph)
  library(ggrepel)
  library(gridExtra)
})

# ==============================================================================
# FEATURE 1 PLOTS: Machine Learning Results
# ==============================================================================

#' Plot ML model performance and feature importance
#' Based on: Breiman (2001) Random Forests, Variable Importance
#' @param ml_result Result from train_haplotype_ranker()
#' @param output_prefix Prefix for output files
plot_ml_results <- function(ml_result, output_prefix = "ml_analysis") {
  
  cat("\n========================================\n")
  cat("PLOTTING ML RESULTS\n")
  cat("========================================\n\n")
  
  # 1. Predicted vs Observed plot
  pred_data <- ml_result$predictions
  
  # Calculate regression line
  lm_fit <- lm(Predicted ~ Actual, data = pred_data)
  r_squared <- ml_result$performance$r_squared
  rmse <- ml_result$performance$rmse
  
  p1 <- ggplot(pred_data, aes(x = Actual, y = Predicted)) +
    geom_point(aes(color = Haplotype), size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = "gray40", size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
    annotate("text", x = Inf, y = -Inf, 
             label = sprintf("R² = %.3f\nRMSE = %.3f", r_squared, rmse),
             hjust = 1.1, vjust = -0.5, size = 5) +
    labs(title = "Machine Learning Model Performance",
         subtitle = "Random Forest Predictions on Test Set",
         x = "Observed Values",
         y = "Predicted Values") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 13, face = "bold")
    )
  
  ggsave(paste0(output_prefix, "_predictions.pdf"), p1, width = 10, height = 8)
  ggsave(paste0(output_prefix, "_predictions.png"), p1, width = 10, height = 8, dpi = 300)
  
  # 2. Haplotype Rankings Bar Plot
  rankings <- ml_result$haplotype_rankings %>%
    arrange(ML_Rank) %>%
    head(15)  # Top 15
  
  p2 <- ggplot(rankings, aes(x = reorder(Haplotype, -Mean_Predicted), 
                              y = Mean_Predicted)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean_Predicted - SD/sqrt(N), 
                      ymax = Mean_Predicted + SD/sqrt(N)),
                  width = 0.3, size = 0.8) +
    geom_text(aes(label = paste0("n=", N)), 
              vjust = -0.5, hjust = 1.2, size = 3.5) +
    labs(title = "Haplotype Rankings by ML Prediction",
         subtitle = "Top 15 Haplotypes (Error bars = SE)",
         x = "Haplotype",
         y = "Predicted Mean Performance") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      panel.grid.major.y = element_line(color = "gray90")
    )
  
  ggsave(paste0(output_prefix, "_rankings.pdf"), p2, width = 12, height = 8)
  ggsave(paste0(output_prefix, "_rankings.png"), p2, width = 12, height = 8, dpi = 300)
  
  cat("  Saved:", paste0(output_prefix, "_predictions.*"), "\n")
  cat("  Saved:", paste0(output_prefix, "_rankings.*"), "\n\n")
  
  return(list(p1 = p1, p2 = p2))
}


# ==============================================================================
# FEATURE 2 PLOTS: Multi-Trait Selection Index
# ==============================================================================

#' Plot multi-trait selection index results
#' Based on: Smith (1936) hazel (1943) selection index theory
#' @param selection_result Result from calculate_selection_index()
#' @param output_prefix Prefix for output files
plot_selection_index <- function(selection_result, output_prefix = "selection_index") {
  
  cat("\n========================================\n")
  cat("PLOTTING SELECTION INDEX\n")
  cat("========================================\n\n")
  
  # 1. Selection Index Rankings
  summary_data <- selection_result$summary %>%
    arrange(Rank) %>%
    head(15)
  
  p1 <- ggplot(summary_data, aes(x = reorder(Haplotype, -Mean_Index), 
                                  y = Mean_Index)) +
    geom_bar(stat = "identity", aes(fill = Mean_Index), color = "black") +
    scale_fill_gradient2(low = "red", mid = "white", high = "darkgreen",
                        midpoint = 0, name = "Index\nValue") +
    geom_text(aes(label = sprintf("%.2f", Mean_Index)), 
              vjust = ifelse(summary_data$Mean_Index >= 0, -0.5, 1.5),
              size = 3.5, fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    labs(title = "Multi-Trait Selection Index Rankings",
         subtitle = paste("Weights:", paste(names(selection_result$weights), 
                         selection_result$weights, sep = "=", collapse = ", ")),
         x = "Haplotype",
         y = "Selection Index (standardized)") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      panel.grid.major.y = element_line(color = "gray90")
    )
  
  ggsave(paste0(output_prefix, "_rankings.pdf"), p1, width = 12, height = 8)
  ggsave(paste0(output_prefix, "_rankings.png"), p1, width = 12, height = 8, dpi = 300)
  
  # 2. Distribution of Index Values
  individual_data <- selection_result$individual_values
  
  p2 <- ggplot(individual_data, aes(x = Selection_Index, fill = Haplotype)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
    labs(title = "Distribution of Selection Index Values",
         subtitle = "Density plot by haplotype",
         x = "Selection Index",
         y = "Density") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "right",
      axis.title = element_text(size = 13, face = "bold")
    )
  
  ggsave(paste0(output_prefix, "_distribution.pdf"), p2, width = 10, height = 8)
  ggsave(paste0(output_prefix, "_distribution.png"), p2, width = 10, height = 8, dpi = 300)
  
  cat("  Saved:", paste0(output_prefix, "_rankings.*"), "\n")
  cat("  Saved:", paste0(output_prefix, "_distribution.*"), "\n\n")
  
  return(list(p1 = p1, p2 = p2))
}


# ==============================================================================
# FEATURE 3 PLOTS: Haplotype Stability
# ==============================================================================

#' Plot haplotype stability analysis
#' Based on: Finlay & Wilkinson (1963), Eberhart & Russell (1966)
#' @param stability_result Result from calculate_haplotype_stability()
#' @param output_prefix Prefix for output files
plot_stability_analysis <- function(stability_result, output_prefix = "stability") {
  
  cat("\n========================================\n")
  cat("PLOTTING STABILITY ANALYSIS\n")
  cat("========================================\n\n")
  
  # 1. Mean vs CV plot (classic stability plot)
  p1 <- ggplot(stability_result, aes(x = Mean_trait, y = CV)) +
    geom_point(aes(size = N_samples, color = Stability_Rank), alpha = 0.7) +
    geom_text_repel(aes(label = Haplotype), size = 3.5, max.overlaps = 20) +
    scale_color_gradient(low = "darkgreen", high = "red", name = "Stability\nRank") +
    scale_size_continuous(range = c(3, 10), name = "Sample\nSize") +
    labs(title = "Haplotype Stability Analysis",
         subtitle = "Lower CV = More Stable",
         x = "Mean Performance",
         y = "Coefficient of Variation (%)") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(paste0(output_prefix, "_mean_vs_cv.pdf"), p1, width = 12, height = 8)
  ggsave(paste0(output_prefix, "_mean_vs_cv.png"), p1, width = 12, height = 8, dpi = 300)
  
  # 2. Stability Rankings Bar Plot
  top_stable <- stability_result %>%
    arrange(Stability_Rank) %>%
    head(15)
  
  p2 <- ggplot(top_stable, aes(x = reorder(Haplotype, Stability_Rank), y = CV)) +
    geom_bar(stat = "identity", aes(fill = Mean_trait), color = "black") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Mean\nPerformance") +
    geom_text(aes(label = sprintf("CV=%.1f%%", CV)), 
              hjust = -0.1, size = 3.5) +
    coord_flip() +
    labs(title = "Most Stable Haplotypes",
         subtitle = "Ranked by Coefficient of Variation (lower = better)",
         x = "Haplotype",
         y = "Coefficient of Variation (%)") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      panel.grid.major.x = element_line(color = "gray90")
    )
  
  ggsave(paste0(output_prefix, "_rankings.pdf"), p2, width = 10, height = 10)
  ggsave(paste0(output_prefix, "_rankings.png"), p2, width = 10, height = 10, dpi = 300)
  
  cat("  Saved:", paste0(output_prefix, "_mean_vs_cv.*"), "\n")
  cat("  Saved:", paste0(output_prefix, "_rankings.*"), "\n\n")
  
  return(list(p1 = p1, p2 = p2))
}


# ==============================================================================
# FEATURE 4 PLOTS: Haplotype Network
# ==============================================================================

#' Plot haplotype relationship network
#' Based on: Bandelt et al. (1999) median-joining networks
#' @param distance_matrix Distance matrix from calculate_haplotype_distances()
#' @param haplotype_sizes Named vector of haplotype sizes
#' @param output_prefix Prefix for output files
plot_haplotype_network <- function(distance_matrix, haplotype_sizes = NULL,
                                  output_prefix = "network") {
  
  cat("\n========================================\n")
  cat("PLOTTING HAPLOTYPE NETWORK\n")
  cat("========================================\n\n")
  
  # Create network from distance matrix
  # Use only edges with distance < threshold to avoid clutter
  threshold <- quantile(distance_matrix[upper.tri(distance_matrix)], 0.5)
  
  # Create edge list
  edges <- data.frame()
  for(i in 1:(nrow(distance_matrix)-1)) {
    for(j in (i+1):ncol(distance_matrix)) {
      if(distance_matrix[i,j] < threshold) {
        edges <- rbind(edges, data.frame(
          from = rownames(distance_matrix)[i],
          to = colnames(distance_matrix)[j],
          weight = 1 - distance_matrix[i,j]  # Convert distance to similarity
        ))
      }
    }
  }
  
  # Create igraph object
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Set node sizes
  if(!is.null(haplotype_sizes)) {
    V(g)$size <- haplotype_sizes[V(g)$name]
  } else {
    V(g)$size <- 5
  }
  
  # Set edge weights
  E(g)$width <- E(g)$weight * 5
  
  # Layout
  layout <- layout_with_fr(g)
  
  # Plot
  pdf(paste0(output_prefix, "_graph.pdf"), width = 12, height = 12)
  plot(g, 
       vertex.size = sqrt(V(g)$size) * 3,
       vertex.color = "lightblue",
       vertex.frame.color = "darkblue",
       vertex.label.cex = 1.2,
       vertex.label.color = "black",
       vertex.label.font = 2,
       edge.width = E(g)$width,
       edge.color = "gray70",
       layout = layout,
       main = "Haplotype Relationship Network",
       sub = "Node size = number of samples, Edge width = genetic similarity")
  dev.off()
  
  png(paste0(output_prefix, "_graph.png"), width = 12, height = 12, 
      units = "in", res = 300)
  plot(g, 
       vertex.size = sqrt(V(g)$size) * 3,
       vertex.color = "lightblue",
       vertex.frame.color = "darkblue",
       vertex.label.cex = 1.2,
       vertex.label.color = "black",
       vertex.label.font = 2,
       edge.width = E(g)$width,
       edge.color = "gray70",
       layout = layout,
       main = "Haplotype Relationship Network",
       sub = "Node size = number of samples, Edge width = genetic similarity")
  dev.off()
  
  # Also create heatmap
  pdf(paste0(output_prefix, "_heatmap.pdf"), width = 10, height = 10)
  pheatmap(distance_matrix,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Pairwise Genetic Distances Between Haplotypes",
           fontsize = 12,
           fontsize_row = 10,
           fontsize_col = 10)
  dev.off()
  
  png(paste0(output_prefix, "_heatmap.png"), width = 10, height = 10, 
      units = "in", res = 300)
  pheatmap(distance_matrix,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Pairwise Genetic Distances Between Haplotypes",
           fontsize = 12,
           fontsize_row = 10,
           fontsize_col = 10)
  dev.off()
  
  cat("  Saved:", paste0(output_prefix, "_graph.*"), "\n")
  cat("  Saved:", paste0(output_prefix, "_heatmap.*"), "\n\n")
}


# ==============================================================================
# FEATURE 5: Enhanced Haplotype Distribution Plots
# ==============================================================================

#' Create customized haplotype distribution plots
#' @param haplotype_table Data frame with Accession and Haplotype
#' @param colors Named vector of colors for haplotypes
#' @param output_prefix Prefix for output files
plot_haplotype_distribution <- function(haplotype_table, colors = NULL,
                                       output_prefix = "distribution") {
  
  cat("\n========================================\n")
  cat("PLOTTING HAPLOTYPE DISTRIBUTION\n")
  cat("========================================\n\n")
  
  # Calculate frequencies
  freq_data <- haplotype_table %>%
    count(Haplotype) %>%
    mutate(Percentage = n / sum(n) * 100) %>%
    arrange(desc(n))
  
  # Set colors if not provided
  if(is.null(colors)) {
    n_haps <- nrow(freq_data)
    colors <- setNames(
      brewer.pal(max(3, min(n_haps, 12)), "Set3")[1:n_haps],
      freq_data$Haplotype
    )
  }
  
  # 1. Bar plot
  p1 <- ggplot(freq_data, aes(x = reorder(Haplotype, -n), y = n, fill = Haplotype)) +
    geom_bar(stat = "identity", color = "black", size = 0.8) +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = paste0(n, "\n(", sprintf("%.1f", Percentage), "%)")),
              vjust = -0.3, size = 3.5, fontface = "bold") +
    labs(title = "Haplotype Frequency Distribution",
         x = "Haplotype",
         y = "Number of Accessions") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "gray90")
    )
  
  ggsave(paste0(output_prefix, "_barplot.pdf"), p1, width = 10, height = 8)
  ggsave(paste0(output_prefix, "_barplot.png"), p1, width = 10, height = 8, dpi = 300)
  
  # 2. Pie chart
  p2 <- ggplot(freq_data, aes(x = "", y = n, fill = Haplotype)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 2) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = paste0(Haplotype, "\n", n, " (", 
                                 sprintf("%.1f", Percentage), "%)")),
              position = position_stack(vjust = 0.5),
              size = 4, fontface = "bold") +
    labs(title = "Haplotype Composition") +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold")
    )
  
  ggsave(paste0(output_prefix, "_piechart.pdf"), p2, width = 10, height = 8)
  ggsave(paste0(output_prefix, "_piechart.png"), p2, width = 10, height = 8, dpi = 300)
  
  cat("  Saved:", paste0(output_prefix, "_barplot.*"), "\n")
  cat("  Saved:", paste0(output_prefix, "_piechart.*"), "\n\n")
  
  return(list(p1 = p1, p2 = p2))
}

cat("\n")
cat("================================================================================\n")
cat("  SuperHap Visualization Module Loaded\n")
cat("================================================================================\n")
cat("\nAvailable plotting functions:\n")
cat("  - plot_ml_results()              : ML model performance & rankings\n")
cat("  - plot_selection_index()         : Multi-trait selection index\n")
cat("  - plot_stability_analysis()      : Environmental stability\n")
cat("  - plot_haplotype_network()       : Genetic distance network\n")
cat("  - plot_haplotype_distribution()  : Frequency distributions\n")
cat("================================================================================\n\n")
