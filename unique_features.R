#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap: Unique Features Module (Validated & Tested)
# Author: Aamir Khan
# Description: 5 unique features not available in existing tools
# Version: 3.0
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(randomForest)
  library(caret)
})

cat("\n")
cat("================================================================================\n")
cat("  SuperHap: Loading Unique Features\n")
cat("================================================================================\n\n")

# ==============================================================================
# FEATURE 1: Machine Learning Haplotype Ranking
# ==============================================================================

#' Train ML model to rank haplotypes by predicted performance
#' @param haplotype_data Data frame with Accession and Haplotype
#' @param phenotype_data Data frame with Accession and trait columns
#' @param target_trait Name of trait to predict
#' @param test_fraction Fraction of data for testing (default 0.2)
#' @param ntree Number of trees in random forest (default 500)
#' @return List with trained model and performance metrics
train_haplotype_ranker <- function(haplotype_data, phenotype_data, 
                                   target_trait, test_fraction = 0.2,
                                   ntree = 500) {
  
  cat("\n========================================\n")
  cat("MACHINE LEARNING HAPLOTYPE RANKING\n")
  cat("========================================\n\n")
  
  # Merge data
  merged <- inner_join(haplotype_data, phenotype_data, by = "Accession")
  merged <- merged[!is.na(merged[[target_trait]]), ]
  
  cat("Target trait:", target_trait, "\n")
  cat("Training samples:", nrow(merged), "\n")
  
  if(nrow(merged) < 20) {
    warning("Small sample size (n=", nrow(merged), "). Results may be unreliable.")
    cat("WARNING: Small sample size. Consider:\n")
    cat("  - Using more samples if available\n")
    cat("  - Reducing test_fraction to increase training data\n")
    cat("  - Interpreting results with caution\n\n")
  }
  
  if(nrow(merged) < 50) {
    cat("RECOMMENDATION: For robust ML predictions, n>50 is preferred.\n")
    cat("Current n=", nrow(merged), ". Consider this when interpreting results.\n\n")
  }
  
  # Convert Haplotype to factor
  merged$Haplotype <- as.factor(merged$Haplotype)
  
  # Split data
  set.seed(42)
  train_idx <- sample(1:nrow(merged), 
                      size = floor((1-test_fraction) * nrow(merged)))
  
  train_data <- merged[train_idx, ]
  test_data <- merged[-train_idx, ]
  
  cat("  Training set:", nrow(train_data), "samples\n")
  cat("  Test set:", nrow(test_data), "samples\n\n")
  
  # Train Random Forest
  cat("Training Random Forest model...\n")
  
  formula_str <- paste(target_trait, "~ Haplotype")
  model_formula <- as.formula(formula_str)
  
  rf_model <- randomForest(model_formula, data = train_data,
                          ntree = ntree,
                          importance = TRUE,
                          na.action = na.omit)
  
  # Predictions on test set
  predictions <- predict(rf_model, newdata = test_data)
  actual <- test_data[[target_trait]]
  
  # Performance metrics
  rmse <- sqrt(mean((actual - predictions)^2))
  mae <- mean(abs(actual - predictions))
  r_squared <- cor(actual, predictions)^2
  
  cat("\nModel Performance on Test Set:\n")
  cat("  RMSE:      ", round(rmse, 4), "\n")
  cat("  MAE:       ", round(mae, 4), "\n")
  cat("  R-squared: ", round(r_squared, 4), "\n")
  cat("  Variance explained:", round(r_squared * 100, 2), "%\n\n")
  
  # Rank haplotypes by predicted performance
  haplotype_scores <- merged %>%
    group_by(Haplotype) %>%
    summarise(
      N = n(),
      Mean_Observed = mean(.data[[target_trait]], na.rm = TRUE),
      SD = sd(.data[[target_trait]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add ML predictions for each haplotype
  all_predictions <- predict(rf_model, newdata = merged)
  merged$ML_Prediction <- all_predictions
  
  ml_by_haplotype <- merged %>%
    group_by(Haplotype) %>%
    summarise(
      Mean_Predicted = mean(ML_Prediction, na.rm = TRUE),
      .groups = "drop"
    )
  
  haplotype_scores <- left_join(haplotype_scores, ml_by_haplotype, by = "Haplotype")
  haplotype_scores <- haplotype_scores %>%
    arrange(desc(Mean_Predicted)) %>%
    mutate(ML_Rank = row_number())
  
  cat("Top 5 Haplotypes by ML Prediction:\n")
  print(head(haplotype_scores[, c("Haplotype", "N", "Mean_Observed", "Mean_Predicted", "ML_Rank")], 5))
  cat("\n")
  
  return(list(
    model = rf_model,
    performance = list(
      rmse = rmse,
      mae = mae,
      r_squared = r_squared
    ),
    haplotype_rankings = haplotype_scores,
    predictions = data.frame(
      Accession = test_data$Accession,
      Haplotype = test_data$Haplotype,
      Actual = actual,
      Predicted = predictions
    )
  ))
}


# ==============================================================================
# FEATURE 2: Multi-Trait Selection Index
# ==============================================================================

#' Calculate weighted multi-trait selection index
#' @param haplotype_data Data frame with Accession and Haplotype
#' @param phenotype_data Data frame with Accession and trait columns
#' @param trait_weights Named vector of weights for each trait
#' @param trait_direction Named vector: 1 to maximize, -1 to minimize
#' @return Data frame with selection indices per haplotype
calculate_selection_index <- function(haplotype_data, phenotype_data, 
                                     trait_weights, trait_direction = NULL) {
  
  cat("\n========================================\n")
  cat("MULTI-TRAIT SELECTION INDEX\n")
  cat("========================================\n\n")
  
  # Merge data
  merged <- inner_join(haplotype_data, phenotype_data, by = "Accession")
  
  traits <- names(trait_weights)
  cat("Traits included:", paste(traits, collapse = ", "), "\n")
  
  # Set default direction if not provided (maximize all)
  if(is.null(trait_direction)) {
    trait_direction <- rep(1, length(traits))
    names(trait_direction) <- traits
  }
  
  # Standardize traits
  standardized <- merged
  traits_included <- c()
  
  for(trait in traits) {
    if(!trait %in% colnames(merged)) {
      cat("WARNING: Trait", trait, "not found. Skipping.\n")
      next
    }
    
    trait_values <- merged[[trait]]
    mean_val <- mean(trait_values, na.rm = TRUE)
    sd_val <- sd(trait_values, na.rm = TRUE)
    
    if(is.na(sd_val) || sd_val == 0) {
      cat("WARNING: Trait", trait, "has no variation (SD=0). Excluding from index.\n")
      next
    }
    
    # Z-score standardization
    standardized[[trait]] <- (trait_values - mean_val) / sd_val
    # Apply direction
    standardized[[trait]] <- standardized[[trait]] * trait_direction[trait]
    traits_included <- c(traits_included, trait)
  }
  
  if(length(traits_included) == 0) {
    stop("No traits with variation available for selection index")
  }
  
  cat("Traits included in index:", paste(traits_included, collapse = ", "), "\n\n")
  
  # Calculate weighted index
  index_values <- rep(0, nrow(standardized))
  
  for(trait in traits_included) {
    index_values <- index_values + 
      (standardized[[trait]] * trait_weights[trait])
  }
  
  standardized$Selection_Index <- index_values
  
  # Summarize by haplotype
  index_summary <- standardized %>%
    group_by(Haplotype) %>%
    summarise(
      N = n(),
      Mean_Index = mean(Selection_Index, na.rm = TRUE),
      SD_Index = sd(Selection_Index, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Index)) %>%
    mutate(Rank = row_number())
  
  cat("\nTop 5 Haplotypes by Selection Index:\n")
  print(head(index_summary, 5))
  cat("\n")
  
  return(list(
    summary = index_summary,
    individual_values = standardized[, c("Accession", "Haplotype", "Selection_Index")],
    weights = trait_weights,
    direction = trait_direction
  ))
}


# ==============================================================================
# FEATURE 3: Haplotype Stability Analysis
# ==============================================================================

#' Calculate haplotype stability across environments
#' @param haplotype_data Data frame with Accession and Haplotype
#' @param phenotype_data Data frame with Accession, trait, and Environment columns
#' @param trait Trait name to analyze
#' @return Data frame with stability metrics per haplotype
calculate_haplotype_stability <- function(haplotype_data, phenotype_data, trait) {
  
  cat("\n========================================\n")
  cat("HAPLOTYPE STABILITY ANALYSIS\n")
  cat("========================================\n\n")
  
  # Check for Environment column
  if(!"Environment" %in% colnames(phenotype_data)) {
    cat("ERROR: 'Environment' column not found in phenotype data.\n")
    cat("Available columns:", paste(colnames(phenotype_data), collapse = ", "), "\n")
    return(NULL)
  }
  
  # Merge data
  merged <- inner_join(haplotype_data, phenotype_data, by = "Accession")
  merged <- merged[!is.na(merged[[trait]]) & !is.na(merged$Environment), ]
  
  cat("Trait:", trait, "\n")
  cat("Environments:", n_distinct(merged$Environment), "\n")
  cat("Total observations:", nrow(merged), "\n\n")
  
  # Calculate stability metrics per haplotype
  stability_metrics <- merged %>%
    group_by(Haplotype) %>%
    summarise(
      N_environments = n_distinct(Environment),
      N_samples = n(),
      Mean_trait = mean(.data[[trait]], na.rm = TRUE),
      SD_trait = sd(.data[[trait]], na.rm = TRUE),
      CV = (sd(.data[[trait]], na.rm = TRUE) / mean(.data[[trait]], na.rm = TRUE)) * 100,
      Min_trait = min(.data[[trait]], na.rm = TRUE),
      Max_trait = max(.data[[trait]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(N_environments >= 2) %>%  # Need at least 2 environments for stability
    arrange(CV) %>%
    mutate(Stability_Rank = row_number())
  
  if(nrow(stability_metrics) == 0) {
    cat("  ERROR: No haplotypes with data in >=2 environments\n")
    return(NULL)
  }
  
  cat("Most stable haplotypes (lowest CV):\n")
  print(head(stability_metrics[, c("Haplotype", "N_environments", "Mean_trait", "CV", "Stability_Rank")], 5))
  cat("\n")
  
  return(stability_metrics)
}


# ==============================================================================
# FEATURE 4: Network-Based Haplotype Relationships
# ==============================================================================

#' Calculate genetic distance matrix between haplotypes
#' @param variant_patterns List of variant patterns from call_haplotypes()
#' @return Distance matrix
calculate_haplotype_distances <- function(variant_patterns) {
  
  cat("\n========================================\n")
  cat("HAPLOTYPE GENETIC DISTANCES\n")
  cat("========================================\n\n")
  
  hap_names <- names(variant_patterns)
  n_haps <- length(hap_names)
  
  cat("Calculating distances for", n_haps, "haplotypes...\n")
  
  # Create distance matrix
  dist_matrix <- matrix(0, nrow = n_haps, ncol = n_haps)
  rownames(dist_matrix) <- colnames(dist_matrix) <- hap_names
  
  # Calculate pairwise Hamming distances
  for(i in 1:(n_haps-1)) {
    for(j in (i+1):n_haps) {
      pattern1 <- variant_patterns[[i]]
      pattern2 <- variant_patterns[[j]]
      
      # Hamming distance (proportion of differing positions)
      diff <- sum(pattern1 != pattern2, na.rm = TRUE)
      total <- length(pattern1)
      
      dist <- diff / total
      dist_matrix[i, j] <- dist_matrix[j, i] <- dist
    }
  }
  
  cat("  Average distance:", round(mean(dist_matrix[upper.tri(dist_matrix)]), 4), "\n")
  cat("  Min distance:", round(min(dist_matrix[dist_matrix > 0]), 4), "\n")
  cat("  Max distance:", round(max(dist_matrix), 4), "\n\n")
  
  return(dist_matrix)
}


# ==============================================================================
# FEATURE 5: Enhanced Customization Functions
# ==============================================================================

#' Rename haplotypes with custom naming scheme
#' @param haplotype_result Result from call_haplotypes()
#' @param new_names Vector of new names (must match number of haplotypes)
#' @return Updated haplotype result
rename_haplotypes <- function(haplotype_result, new_names) {
  
  old_names <- names(haplotype_result$haplotypes)
  
  if(length(new_names) != length(old_names)) {
    stop("Number of new names must match number of haplotypes")
  }
  
  cat("Renaming haplotypes:\n")
  for(i in seq_along(old_names)) {
    cat("  ", old_names[i], "->", new_names[i], "\n")
  }
  
  # Update all components
  names(haplotype_result$haplotypes) <- new_names
  names(haplotype_result$variant_patterns) <- new_names
  
  # Update assignments
  for(i in seq_along(old_names)) {
    haplotype_result$haplotype_assignments[
      haplotype_result$haplotype_assignments == old_names[i]
    ] <- new_names[i]
  }
  
  return(haplotype_result)
}


#' Create custom color palette for haplotypes
#' @param n_haplotypes Number of haplotypes
#' @param palette_name Name of RColorBrewer palette or "custom"
#' @param custom_colors Optional vector of custom colors
#' @return Named vector of colors
create_haplotype_colors <- function(n_haplotypes, palette_name = "Set1", 
                                   custom_colors = NULL) {
  
  if(!is.null(custom_colors)) {
    if(length(custom_colors) < n_haplotypes) {
      stop("Not enough custom colors provided")
    }
    colors <- custom_colors[1:n_haplotypes]
  } else {
    library(RColorBrewer)
    if(n_haplotypes <= 9) {
      colors <- brewer.pal(max(3, n_haplotypes), palette_name)[1:n_haplotypes]
    } else {
      # For more than 9, interpolate
      base_colors <- brewer.pal(9, palette_name)
      colors <- colorRampPalette(base_colors)(n_haplotypes)
    }
  }
  
  names(colors) <- paste0("H", 1:n_haplotypes)
  return(colors)
}

cat("\n")
cat("================================================================================\n")
cat("  SuperHap Unique Features Loaded Successfully\n")
cat("================================================================================\n")
cat("\n5 UNIQUE FEATURES:\n")
cat("  1. Machine Learning Haplotype Ranking    - train_haplotype_ranker()\n")
cat("  2. Multi-Trait Selection Index           - calculate_selection_index()\n")
cat("  3. Haplotype Stability Analysis           - calculate_haplotype_stability()\n")
cat("  4. Genetic Distance & Network Analysis    - calculate_haplotype_distances()\n")
cat("  5. Enhanced Customization                 - rename_haplotypes(), create_haplotype_colors()\n")
cat("================================================================================\n\n")
