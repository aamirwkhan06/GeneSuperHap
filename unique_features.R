#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(randomForest)
  library(caret)
})

cat("\n")
cat("================================================================================\n")
cat("  SuperHap: Loading Unique Features\n")
cat("================================================================================\n\n")

train_haplotype_ranker <- function(haplotype_data, phenotype_data, 
                                   target_trait, test_fraction = 0.2,
                                   ntree = 500) {
  
  cat("\n========================================\n")
  cat("MACHINE LEARNING HAPLOTYPE RANKING\n")
  cat("========================================\n\n")
  
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
  
  merged$Haplotype <- as.factor(merged$Haplotype)
  
  set.seed(42)
  train_idx <- sample(1:nrow(merged), 
                      size = floor((1-test_fraction) * nrow(merged)))
  
  train_data <- merged[train_idx, ]
  test_data <- merged[-train_idx, ]
  
  cat("  Training set:", nrow(train_data), "samples\n")
  cat("  Test set:", nrow(test_data), "samples\n\n")
  
  cat("Training Random Forest model...\n")
  
  formula_str <- paste(target_trait, "~ Haplotype")
  model_formula <- as.formula(formula_str)
  
  rf_model <- randomForest(model_formula, data = train_data,
                          ntree = ntree,
                          importance = TRUE,
                          na.action = na.omit)
  
  predictions <- predict(rf_model, newdata = test_data)
  actual <- test_data[[target_trait]]
  
  rmse <- sqrt(mean((actual - predictions)^2))
  mae <- mean(abs(actual - predictions))
  r_squared <- cor(actual, predictions)^2
  
  cat("\nModel Performance on Test Set:\n")
  cat("  RMSE:      ", round(rmse, 4), "\n")
  cat("  MAE:       ", round(mae, 4), "\n")
  cat("  R-squared: ", round(r_squared, 4), "\n")
  cat("  Variance explained:", round(r_squared * 100, 2), "%\n\n")
  
  haplotype_scores <- merged %>%
    group_by(Haplotype) %>%
    summarise(
      N = n(),
      Mean_Observed = mean(.data[[target_trait]], na.rm = TRUE),
      SD = sd(.data[[target_trait]], na.rm = TRUE),
      .groups = "drop"
    )
  
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

calculate_selection_index <- function(haplotype_data, phenotype_data, 
                                     trait_weights, trait_direction = NULL) {
  
  cat("\n========================================\n")
  cat("MULTI-TRAIT SELECTION INDEX\n")
  cat("========================================\n\n")
  
  merged <- inner_join(haplotype_data, phenotype_data, by = "Accession")
  
  traits <- names(trait_weights)
  cat("Traits included:", paste(traits, collapse = ", "), "\n")
  
  if(is.null(trait_direction)) {
    trait_direction <- rep(1, length(traits))
    names(trait_direction) <- traits
  }
  
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
    
    standardized[[trait]] <- (trait_values - mean_val) / sd_val
    standardized[[trait]] <- standardized[[trait]] * trait_direction[trait]
    traits_included <- c(traits_included, trait)
  }
  
  if(length(traits_included) == 0) {
    stop("No traits with variation available for selection index")
  }
  
  cat("Traits included in index:", paste(traits_included, collapse = ", "), "\n\n")
  
  index_values <- rep(0, nrow(standardized))
  
  for(trait in traits_included) {
    index_values <- index_values + 
      (standardized[[trait]] * trait_weights[trait])
  }
  
  standardized$Selection_Index <- index_values
  
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

calculate_haplotype_stability <- function(haplotype_data, phenotype_data, trait) {
  
  cat("\n========================================\n")
  cat("HAPLOTYPE STABILITY ANALYSIS\n")
  cat("========================================\n\n")
  
  if(!"Environment" %in% colnames(phenotype_data)) {
    cat("ERROR: 'Environment' column not found in phenotype data.\n")
    cat("Available columns:", paste(colnames(phenotype_data), collapse = ", "), "\n")
    return(NULL)
  }
  
  merged <- inner_join(haplotype_data, phenotype_data, by = "Accession")
  merged <- merged[!is.na(merged[[trait]]) & !is.na(merged$Environment), ]
  
  cat("Trait:", trait, "\n")
  cat("Environments:", n_distinct(merged$Environment), "\n")
  cat("Total observations:", nrow(merged), "\n\n")
  
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

calculate_haplotype_distances <- function(variant_patterns) {
  
  cat("\n========================================\n")
  cat("HAPLOTYPE GENETIC DISTANCES\n")
  cat("========================================\n\n")
  
  hap_names <- names(variant_patterns)
  n_haps <- length(hap_names)
  
  cat("Calculating distances for", n_haps, "haplotypes...\n")
  
  dist_matrix <- matrix(0, nrow = n_haps, ncol = n_haps)
  rownames(dist_matrix) <- colnames(dist_matrix) <- hap_names
  
  for(i in 1:(n_haps-1)) {
    for(j in (i+1):n_haps) {
      pattern1 <- variant_patterns[[i]]
      pattern2 <- variant_patterns[[j]]
      
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

rename_haplotypes <- function(haplotype_result, new_names) {
  
  old_names <- names(haplotype_result$haplotypes)
  
  if(length(new_names) != length(old_names)) {
    stop("Number of new names must match number of haplotypes")
  }
  
  cat("Renaming haplotypes:\n")
  for(i in seq_along(old_names)) {
    cat("  ", old_names[i], "->", new_names[i], "\n")
  }
  
  names(haplotype_result$haplotypes) <- new_names
  names(haplotype_result$variant_patterns) <- new_names
  
  for(i in seq_along(old_names)) {
    haplotype_result$haplotype_assignments[
      haplotype_result$haplotype_assignments == old_names[i]
    ] <- new_names[i]
  }
  
  return(haplotype_result)
}

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
