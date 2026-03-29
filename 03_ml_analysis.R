#!/usr/bin/env Rscript

source("core_functions.R")
source("unique_features.R")
source("visualization.R")

cat("\n")
cat("================================================================================\n")
cat("  GeneSuperHap Pipeline - STEP 03: MACHINE LEARNING ANALYSIS\n")
cat("================================================================================\n\n")

HAPLOTYPE_DATA <- "results/haplotypes_data.RData"
PHENOTYPE_FILE <- "input_data/phenotypes.txt"
OUTPUT_DIR <- "results"
OUTPUT_PREFIX <- "ml_analysis"

TARGET_TRAIT <- "Yield"        # Change this to your trait of interest
TEST_FRACTION <- 0.2           # 20% for testing
N_TREES <- 500                 # Number of trees in Random Forest

cat("Loading data...\n")
load(HAPLOTYPE_DATA)
phenotypes <- read.table(PHENOTYPE_FILE, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE)

cat("  Haplotypes:", nrow(hap_table), "accessions\n")
cat("  Phenotypes:", nrow(phenotypes), "accessions\n\n")

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

cat("\nSaving ML results...\n")

write.table(
  ml_model$haplotype_rankings,
  file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_rankings.txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  ml_model$predictions,
  file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_predictions.txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

save(
  ml_model,
  file = file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_model.RData"))
)

cat("  Results saved\n\n")

cat("Generating visualizations...\n")

plots <- plot_ml_results(
  ml_result = ml_model,
  output_prefix = file.path(OUTPUT_DIR, OUTPUT_PREFIX)
)

cat("  Plots saved\n\n")
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
