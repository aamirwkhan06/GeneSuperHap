#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap Pipeline - Step 05: Network Analysis & Final Visualizations
# Author: Aamir Khan
# Description: Genetic distance network and haplotype distribution plots
# ==============================================================================

# Load required functions
source("core_functions.R")
source("unique_features.R")
source("visualization.R")

cat("\n")
cat("================================================================================\n")
cat("  SUPERHAP PIPELINE - STEP 05: NETWORK ANALYSIS & VISUALIZATION\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

HAPLOTYPE_DATA <- "results/haplotypes_data.RData"
OUTPUT_DIR <- "results"

# Custom colors (optional - set to NULL for automatic)
CUSTOM_COLORS <- NULL
# Example custom colors:
# CUSTOM_COLORS <- c("H1" = "#FF0000", "H2" = "#00FF00", "H3" = "#0000FF")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading haplotype data...\n")
load(HAPLOTYPE_DATA)

cat("  Haplotypes:", length(haplotypes$haplotypes), "\n")
cat("  Total accessions:", nrow(hap_table), "\n\n")

# ==============================================================================
# PART 1: GENETIC DISTANCE NETWORK
# ==============================================================================

cat("PART 1: Calculating Genetic Distances\n")
cat("========================================\n")

# Calculate distance matrix
distances <- calculate_haplotype_distances(
  variant_patterns = haplotypes$variant_patterns
)

# Save distance matrix
write.table(
  distances,
  file.path(OUTPUT_DIR, "haplotype_distances.txt"),
  sep = "\t",
  quote = FALSE
)

cat("  Distance matrix saved\n\n")

# Create haplotype sizes for network plot
haplotype_sizes <- sapply(haplotypes$haplotypes, length)

# Generate network plots
cat("Generating network visualizations...\n")
plot_haplotype_network(
  distance_matrix = distances,
  haplotype_sizes = haplotype_sizes,
  output_prefix = file.path(OUTPUT_DIR, "network")
)

cat("  Network plots saved\n\n")

# ==============================================================================
# PART 2: HAPLOTYPE DISTRIBUTION PLOTS
# ==============================================================================

cat("PART 2: Creating Distribution Plots\n")
cat("========================================\n")

# Create custom colors if not provided
if(is.null(CUSTOM_COLORS)) {
  n_haplotypes <- length(unique(hap_table$Haplotype))
  CUSTOM_COLORS <- create_haplotype_colors(
    n_haplotypes = n_haplotypes,
    palette_name = "Set3"
  )
}

# Generate distribution plots
cat("Generating distribution visualizations...\n")
dist_plots <- plot_haplotype_distribution(
  haplotype_table = hap_table,
  colors = CUSTOM_COLORS,
  output_prefix = file.path(OUTPUT_DIR, "distribution")
)

cat("  Distribution plots saved\n\n")

# ==============================================================================
# PART 3: COMPREHENSIVE SUMMARY REPORT
# ==============================================================================

cat("PART 3: Creating Comprehensive Summary Report\n")
cat("========================================\n")

# Create summary report
report_file <- file.path(OUTPUT_DIR, "ANALYSIS_SUMMARY.txt")
report <- file(report_file, "w")

writeLines("================================================================================", report)
writeLines("                    SUPERHAP ANALYSIS SUMMARY REPORT", report)
writeLines("================================================================================", report)
writeLines("", report)
writeLines(paste("Analysis Date:", Sys.Date()), report)
writeLines("", report)

writeLines("DATA SUMMARY:", report)
writeLines("----------------------------------------", report)
writeLines(paste("Total variants analyzed:", nrow(vcf_data$variants)), report)
writeLines(paste("Chromosome(s):", unique(vcf_data$variants$CHROM)), report)
writeLines(paste("Position range:", min(vcf_data$variants$POS), "-", 
                max(vcf_data$variants$POS)), report)
writeLines(paste("Initial samples:", length(vcf_data$samples)), report)
writeLines(paste("Excluded samples:", length(filtered$removed_samples)), report)
writeLines(paste("Retained samples:", length(filtered$kept_samples)), report)
writeLines("", report)

writeLines("HAPLOTYPE SUMMARY:", report)
writeLines("----------------------------------------", report)
writeLines(paste("Number of haplotypes:", length(haplotypes$haplotypes)), report)
writeLines("", report)

hap_freq <- hap_table %>%
  count(Haplotype) %>%
  mutate(Percentage = round(n / sum(n) * 100, 2)) %>%
  arrange(desc(n))

writeLines("Haplotype Frequencies:", report)
for(i in 1:nrow(hap_freq)) {
  writeLines(sprintf("  %s: %d samples (%.2f%%)", 
                    hap_freq$Haplotype[i], 
                    hap_freq$n[i], 
                    hap_freq$Percentage[i]), report)
}
writeLines("", report)

writeLines("GENETIC DIVERSITY METRICS:", report)
writeLines("----------------------------------------", report)
writeLines(sprintf("Average pairwise distance: %.4f", 
                  mean(distances[upper.tri(distances)])), report)
writeLines(sprintf("Minimum distance: %.4f", 
                  min(distances[distances > 0])), report)
writeLines(sprintf("Maximum distance: %.4f", 
                  max(distances)), report)
writeLines("", report)

writeLines("OUTPUT FILES GENERATED:", report)
writeLines("----------------------------------------", report)
writeLines("Haplotype Calling:", report)
writeLines("  - haplotypes_accessions.txt", report)
writeLines("  - haplotypes_summary.txt", report)
writeLines("  - haplotypes_data.RData", report)
writeLines("", report)

writeLines("Statistical Analysis:", report)
writeLines("  - superior_haplotypes_summary.txt", report)
writeLines("  - superior_haplotypes_*_stats.txt (per trait)", report)
writeLines("", report)

if(file.exists(file.path(OUTPUT_DIR, "ml_analysis_rankings.txt"))) {
  writeLines("Machine Learning:", report)
  writeLines("  - ml_analysis_rankings.txt", report)
  writeLines("  - ml_analysis_predictions.*", report)
  writeLines("  - ml_analysis_model.RData", report)
  writeLines("", report)
}

if(file.exists(file.path(OUTPUT_DIR, "selection_index_summary.txt"))) {
  writeLines("Selection Index:", report)
  writeLines("  - selection_index_summary.txt", report)
  writeLines("  - selection_index_rankings.*", report)
  writeLines("", report)
}

if(file.exists(file.path(OUTPUT_DIR, "stability_analysis.txt"))) {
  writeLines("Stability Analysis:", report)
  writeLines("  - stability_analysis.txt", report)
  writeLines("  - stability_mean_vs_cv.*", report)
  writeLines("", report)
}

writeLines("Network & Visualization:", report)
writeLines("  - haplotype_distances.txt", report)
writeLines("  - network_graph.*", report)
writeLines("  - network_heatmap.*", report)
writeLines("  - distribution_barplot.*", report)
writeLines("  - distribution_piechart.*", report)
writeLines("", report)

writeLines("================================================================================", report)
writeLines("                          ANALYSIS COMPLETE", report)
writeLines("================================================================================", report)

close(report)

cat("  Summary report saved:", report_file, "\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  SUPERHAP PIPELINE COMPLETE!\n")
cat("================================================================================\n")
cat("\nAll analyses finished successfully.\n")
cat("\nGenerated files in 'results/' directory:\n")
cat("  - Haplotype calling results\n")
cat("  - Statistical comparisons\n")
cat("  - Machine learning rankings\n")
cat("  - Multi-trait selection index\n")
cat("  - Stability analysis (if multi-environment data)\n")
cat("  - Genetic distance network\n")
cat("  - Distribution visualizations\n")
cat("  - Comprehensive summary report\n")
cat("\nCheck ANALYSIS_SUMMARY.txt for complete details.\n")
cat("================================================================================\n\n")
