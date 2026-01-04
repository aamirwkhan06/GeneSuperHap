#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap: Master Pipeline Runner
# Author: Aamir Khan
# Description: Run complete pipeline with validation
# ==============================================================================

cat("
################################################################################
                            SUPERHAP PIPELINE
                     Machine Learning Haplotype Analysis
################################################################################

Author: Aamir Khan
Version: 1.0
Date: ", format(Sys.Date(), "%B %d, %Y"), "

################################################################################
\n")

# ==============================================================================
# SETUP
# ==============================================================================

cat("Checking setup...\n")

# Check if running from correct directory
if(!file.exists("core_functions.R")) {
  stop("ERROR: Please run this script from the superhap directory")
}

# Check for required packages
required_packages <- c("dplyr", "tidyr", "ggplot2", "RColorBrewer", 
                      "randomForest", "caret", "pheatmap", "igraph", 
                      "ggrepel", "gridExtra")

missing_packages <- c()
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if(length(missing_packages) > 0) {
  cat("\nERROR: Missing required packages:\n")
  cat("  ", paste(missing_packages, collapse = ", "), "\n\n")
  cat("Please install them using:\n")
  cat("  install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n\n")
  stop("Missing packages")
}

cat("  ✓ All required packages installed\n")

# Check for input files
if(!file.exists("input_data/variants.vcf") && 
   !file.exists("input_data/variants.vcf.gz")) {
  cat("\nWARNING: No VCF file found in input_data/\n")
  cat("Using test data instead...\n")
  
  # Copy test data to input_data
  if(!dir.exists("input_data")) {
    dir.create("input_data")
  }
  
  file.copy("test_data/variants.vcf", "input_data/variants.vcf", overwrite = TRUE)
  file.copy("test_data/phenotypes.txt", "input_data/phenotypes.txt", overwrite = TRUE)
  
  cat("  ✓ Test data copied to input_data/\n")
}

cat("  ✓ Input files ready\n\n")

# ==============================================================================
# RUN PIPELINE
# ==============================================================================

# Track execution time
pipeline_start <- Sys.time()

scripts <- c(
  "01_haplotype_calling.R",
  "02_superior_haplotypes.R",
  "03_ml_analysis.R",
  "04_selection_index.R",
  "05_network_analysis.R"
)

results <- data.frame(
  Script = scripts,
  Status = character(length(scripts)),
  Runtime = numeric(length(scripts)),
  stringsAsFactors = FALSE
)

for(i in seq_along(scripts)) {
  script <- scripts[i]
  
  cat("\n")
  cat("################################################################################\n")
  cat("STEP", i, "of", length(scripts), ":", script, "\n")
  cat("################################################################################\n\n")
  
  start_time <- Sys.time()
  
  tryCatch({
    source(script)
    results$Status[i] <- "SUCCESS"
  }, error = function(e) {
    cat("\n\nERROR in", script, ":\n")
    cat(as.character(e), "\n\n")
    results$Status[i] <- "FAILED"
  }, warning = function(w) {
    cat("\nWARNING in", script, ":\n")
    cat(as.character(w), "\n\n")
  })
  
  end_time <- Sys.time()
  results$Runtime[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("\nStep", i, "completed in", round(results$Runtime[i], 2), "seconds\n")
  
  # Stop if step failed
  if(results$Status[i] == "FAILED") {
    cat("\nPipeline stopped due to error.\n")
    break
  }
}

pipeline_end <- Sys.time()
total_time <- as.numeric(difftime(pipeline_end, pipeline_start, units = "secs"))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n\n")
cat("################################################################################\n")
cat("                         PIPELINE EXECUTION SUMMARY\n")
cat("################################################################################\n\n")

for(i in 1:nrow(results)) {
  status_symbol <- ifelse(results$Status[i] == "SUCCESS", "✓", "✗")
  cat(sprintf("  %s %-35s [%7s] %6.2fs\n", 
              status_symbol, 
              results$Script[i], 
              results$Status[i],
              results$Runtime[i]))
}

cat("\n")
cat("Total execution time:", round(total_time, 2), "seconds\n")
cat("Successful steps:", sum(results$Status == "SUCCESS"), "/", nrow(results), "\n")

if(any(results$Status == "FAILED")) {
  cat("\nWARNING: Some steps failed. Check error messages above.\n")
} else {
  cat("\n✓ PIPELINE COMPLETED SUCCESSFULLY!\n")
  cat("\nResults are in the 'results/' directory.\n")
  cat("\nKey files to check:\n")
  cat("  - results/ANALYSIS_SUMMARY.txt (complete overview)\n")
  cat("  - results/ml_analysis_rankings.txt (ML predictions)\n")
  cat("  - results/selection_index_summary.txt (breeding recommendations)\n")
  cat("  - results/*.pdf (visualizations)\n")
}

cat("\n################################################################################\n\n")

# Save execution log
log_file <- file.path("results", paste0("pipeline_log_", 
                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
writeLines(capture.output({
  cat("SuperHap Pipeline Execution Log\n")
  cat("================================\n\n")
  cat("Start time:", format(pipeline_start), "\n")
  cat("End time:", format(pipeline_end), "\n")
  cat("Total runtime:", round(total_time, 2), "seconds\n\n")
  cat("Step-by-step results:\n")
  print(results)
}), log_file)

cat("Execution log saved:", log_file, "\n\n")
