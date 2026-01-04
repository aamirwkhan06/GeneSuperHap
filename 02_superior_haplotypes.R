#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap Pipeline - Step 02: Superior Haplotype Identification
# Author: Aamir Khan
# Description: Statistical analysis to identify best haplotypes per trait
# ==============================================================================

# Load required functions
source("core_functions.R")
library(dplyr)

cat("\n")
cat("================================================================================\n")
cat("  SUPERHAP PIPELINE - STEP 02: SUPERIOR HAPLOTYPE IDENTIFICATION\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

HAPLOTYPE_DATA <- "results/haplotypes_data.RData"
PHENOTYPE_FILE <- "input_data/phenotypes.txt"
OUTPUT_DIR <- "results"
OUTPUT_PREFIX <- "superior_haplotypes"
MIN_SAMPLES <- 3  # Minimum samples per haplotype for analysis

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading haplotype data...\n")
load(HAPLOTYPE_DATA)

cat("Loading phenotype data...\n")
phenotypes <- read.table(PHENOTYPE_FILE, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE)

cat("  Phenotype samples:", nrow(phenotypes), "\n")
cat("  Traits:", paste(setdiff(colnames(phenotypes), "Accession"), collapse = ", "), "\n\n")

# ==============================================================================
# MERGE DATA
# ==============================================================================

cat("Merging haplotype and phenotype data...\n")
merged <- inner_join(hap_table, phenotypes, by = "Accession")
cat("  Merged samples:", nrow(merged), "\n\n")

# ==============================================================================
# ANALYZE EACH TRAIT
# ==============================================================================

# Get trait columns
trait_cols <- setdiff(colnames(phenotypes), c("Accession", "Environment"))

cat("Analyzing", length(trait_cols), "traits:\n")
cat("----------------------------------------\n\n")

results_list <- list()

for(trait in trait_cols) {
  
  cat("Trait:", trait, "\n")
  
  # Filter valid data
  trait_data <- merged %>%
    filter(!is.na(.data[[trait]])) %>%
    select(Accession, Haplotype, all_of(trait))
  
  if(nrow(trait_data) == 0) {
    cat("  No valid data. Skipping.\n\n")
    next
  }
  
  # Calculate haplotype statistics
  hap_stats <- trait_data %>%
    group_by(Haplotype) %>%
    summarise(
      N = n(),
      Mean = mean(.data[[trait]], na.rm = TRUE),
      SD = sd(.data[[trait]], na.rm = TRUE),
      SE = SD / sqrt(N),
      Min = min(.data[[trait]], na.rm = TRUE),
      Max = max(.data[[trait]], na.rm = TRUE),
      Median = median(.data[[trait]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(N >= MIN_SAMPLES) %>%
    arrange(desc(Mean))
  
  if(nrow(hap_stats) == 0) {
    cat("  No haplotypes with sufficient samples. Skipping.\n\n")
    next
  }
  
  # Identify superior haplotype
  superior_hap <- hap_stats$Haplotype[1]
  superior_mean <- hap_stats$Mean[1]
  
  cat("  Haplotypes analyzed:", nrow(hap_stats), "\n")
  cat("  Superior haplotype:", superior_hap, "\n")
  cat("  Mean value:", round(superior_mean, 3), "\n")
  
  # Pairwise t-tests
  haplotypes_list <- unique(trait_data$Haplotype)
  pairwise_results <- data.frame()
  
  if(length(haplotypes_list) >= 2) {
    for(i in 1:(length(haplotypes_list)-1)) {
      for(j in (i+1):length(haplotypes_list)) {
        hap1 <- haplotypes_list[i]
        hap2 <- haplotypes_list[j]
        
        vals1 <- trait_data %>% filter(Haplotype == hap1) %>% pull(trait)
        vals2 <- trait_data %>% filter(Haplotype == hap2) %>% pull(trait)
        
        if(length(vals1) >= MIN_SAMPLES && length(vals2) >= MIN_SAMPLES) {
          test_result <- t.test(vals1, vals2)
          
          pairwise_results <- rbind(pairwise_results, data.frame(
            Hap1 = hap1,
            Hap2 = hap2,
            N1 = length(vals1),
            N2 = length(vals2),
            Mean1 = mean(vals1),
            Mean2 = mean(vals2),
            Difference = mean(vals1) - mean(vals2),
            P_value = test_result$p.value,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Apply Benjamini-Hochberg FDR correction for multiple testing
    if(nrow(pairwise_results) > 0) {
      pairwise_results$P_adjusted <- p.adjust(pairwise_results$P_value, method = "BH")
      pairwise_results$Significant <- ifelse(pairwise_results$P_adjusted < 0.05, "***", 
                                            ifelse(pairwise_results$P_adjusted < 0.1, "*", "ns"))
      cat("  Multiple testing correction applied (Benjamini-Hochberg FDR)\n")
    }
  }
  
  # Store results
  results_list[[trait]] <- list(
    statistics = hap_stats,
    superior_haplotype = superior_hap,
    superior_mean = superior_mean,
    pairwise = pairwise_results
  )
  
  # Save individual trait results
  write.table(
    hap_stats,
    file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_", trait, "_stats.txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  if(nrow(pairwise_results) > 0) {
    write.table(
      pairwise_results,
      file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_", trait, "_pairwise.txt")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  cat("  Results saved\n\n")
}

# ==============================================================================
# SUMMARY ACROSS TRAITS
# ==============================================================================

cat("Creating summary across all traits...\n")

summary_table <- data.frame()

for(trait in names(results_list)) {
  result <- results_list[[trait]]
  
  summary_table <- rbind(summary_table, data.frame(
    Trait = trait,
    Superior_Haplotype = result$superior_haplotype,
    Mean_Value = round(result$superior_mean, 3),
    N_Haplotypes = nrow(result$statistics),
    Total_Samples = sum(result$statistics$N),
    stringsAsFactors = FALSE
  ))
}

write.table(
  summary_table,
  file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_summary.txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STEP 02 COMPLETE\n")
cat("================================================================================\n")
cat("\nSummary of Superior Haplotypes:\n")
print(summary_table)
cat("\nNext step: Run 03_ml_analysis.R\n")
cat("================================================================================\n\n")
