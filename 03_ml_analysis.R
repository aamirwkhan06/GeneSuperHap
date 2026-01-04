#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap Pipeline - Step 03: Machine Learning Haplotype Ranking
# Author: Aamir Khan
# Description: Train ML models to predict and rank haplotypes
# ==============================================================================

# Load required functions
source("core_functions.R")
source("unique_features.R")
source("visualization.R")

cat("\n")
cat("================================================================================\n")
cat("  SUPERHAP PIPELINE - STEP 03: MACHINE LEARNING ANALYSIS\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

HAPLOTYPE_DATA <- "results/haplotypes_data.RData"
PHENOTYPE_FILE <- "input_data/phenotypes.txt"
OUTPUT_DIR <- "results"
OUTPUT_PREFIX <- "ml_analysis"

# ML parameters
TARGET_TRAIT <- "Yield"        # Change this to your trait of interest
TEST_FRACTION <- 0.2           # 20% for testing
N_TREES <- 500                 # Number of trees in Random Forest

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
# TRAIN ML MODEL
# ==============================================================================

cat("Training Machine Learning model...\n")
cat("Target trait:", TARGET_TRAIT, "\n")
cat("----------------------------------------\n")

ml_model <- train_haplotype_ranker(
  haplotype_data = hap_table,
  phenotype_data = phenotypes,
  target_trait = TARGET_TRAIT,
  test_fraction = TEST_FRACTION,
  ntree = N_TREES
)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\nSaving ML results...\n")

# Save haplotype rankings
write.table(
  ml_model$haplotype_rankings,
  file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_rankings.txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Save predictions
write.table(
  ml_model$predictions,
  file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_predictions.txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Save model object
save(
  ml_model,
  file = file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_model.RData"))
)

cat("  Results saved\n\n")

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

cat("Generating visualizations...\n")

plots <- plot_ml_results(
  ml_result = ml_model,
  output_prefix = file.path(OUTPUT_DIR, OUTPUT_PREFIX)
)

cat("  Plots saved\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STEP 03 COMPLETE\n")
cat("================================================================================\n")
cat("\nMachine Learning Model Performance:\n")
cat("  R-squared:", round(ml_model$performance$r_squared, 4), "\n")
cat("  RMSE:", round(ml_model$performance$rmse, 4), "\n")
cat("  MAE:", round(ml_model$performance$mae, 4), "\n")
cat("\nTop 5 Haplotypes by ML Prediction:\n")
print(head(ml_model$haplotype_rankings[, c("Haplotype", "N", "Mean_Predicted", "ML_Rank")], 5))
cat("\nNext step: Run 04_selection_index.R\n")
cat("================================================================================\n\n")
