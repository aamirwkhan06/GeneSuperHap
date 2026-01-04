#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap Pipeline - Step 04: Multi-Trait Selection & Stability
# Author: Aamir Khan
# Description: Calculate selection indices and stability analysis
# ==============================================================================

# Load required functions
source("core_functions.R")
source("unique_features.R")
source("visualization.R")

cat("\n")
cat("================================================================================\n")
cat("  SUPERHAP PIPELINE - STEP 04: SELECTION INDEX & STABILITY\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

HAPLOTYPE_DATA <- "results/haplotypes_data.RData"
PHENOTYPE_FILE <- "input_data/phenotypes.txt"
OUTPUT_DIR <- "results"

# Selection index parameters
# Customize these for your breeding goals!
TRAIT_WEIGHTS <- c(
  Yield = 2.0,      # Highest weight
  Protein = 1.5,    # Medium-high weight
  Oil = 1.0         # Standard weight
)

TRAIT_DIRECTION <- c(
  Yield = 1,        # Maximize
  Protein = 1,      # Maximize
  Oil = -1          # Minimize (example)
)

# Stability analysis
STABILITY_TRAIT <- "Yield"  # Change to your trait of interest

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")
load(HAPLOTYPE_DATA)
phenotypes <- read.table(PHENOTYPE_FILE, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE)

cat("  Haplotypes:", nrow(hap_table), "accessions\n")
cat("  Phenotypes:", nrow(phenotypes), "accessions\n\n")

# ==============================================================================
# PART 1: MULTI-TRAIT SELECTION INDEX
# ==============================================================================

cat("PART 1: Calculating Multi-Trait Selection Index\n")
cat("========================================\n")

selection_result <- calculate_selection_index(
  haplotype_data = hap_table,
  phenotype_data = phenotypes,
  trait_weights = TRAIT_WEIGHTS,
  trait_direction = TRAIT_DIRECTION
)

# Save results
write.table(
  selection_result$summary,
  file.path(OUTPUT_DIR, "selection_index_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  selection_result$individual_values,
  file.path(OUTPUT_DIR, "selection_index_individual.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Generate plots
cat("\nGenerating selection index plots...\n")
selection_plots <- plot_selection_index(
  selection_result = selection_result,
  output_prefix = file.path(OUTPUT_DIR, "selection_index")
)

cat("  Selection index results saved\n\n")

# ==============================================================================
# PART 2: HAPLOTYPE STABILITY ANALYSIS
# ==============================================================================

cat("PART 2: Haplotype Stability Analysis\n")
cat("========================================\n")

# Check if Environment column exists
if("Environment" %in% colnames(phenotypes)) {
  
  stability_result <- calculate_haplotype_stability(
    haplotype_data = hap_table,
    phenotype_data = phenotypes,
    trait = STABILITY_TRAIT
  )
  
  if(!is.null(stability_result)) {
    # Save results
    write.table(
      stability_result,
      file.path(OUTPUT_DIR, "stability_analysis.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    # Generate plots
    cat("\nGenerating stability plots...\n")
    stability_plots <- plot_stability_analysis(
      stability_result = stability_result,
      output_prefix = file.path(OUTPUT_DIR, "stability")
    )
    
    cat("  Stability results saved\n\n")
  }
  
} else {
  cat("  Note: No 'Environment' column found in phenotype data.\n")
  cat("  Skipping stability analysis.\n")
  cat("  To enable this analysis, add an 'Environment' column to your phenotype file.\n\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STEP 04 COMPLETE\n")
cat("================================================================================\n")
cat("\nTop 5 Haplotypes by Selection Index:\n")
print(head(selection_result$summary, 5))

if(exists("stability_result") && !is.null(stability_result)) {
  cat("\nMost Stable Haplotypes:\n")
  print(head(stability_result[, c("Haplotype", "Mean_trait", "CV", "Stability_Rank")], 5))
}

cat("\nNext step: Run 05_network_analysis.R\n")
cat("================================================================================\n\n")
