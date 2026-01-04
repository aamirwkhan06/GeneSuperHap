#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap: Core Functions (Validated & Tested)
# Author: Aamir Khan
# Description: Independent VCF parsing and haplotype calling
# Version: 3.0 (Realistic features only)
# ==============================================================================

#' Read and parse VCF file
#' @param vcf_file Path to VCF file (can be gzipped)
#' @return List containing genotype matrix and variant information
read_vcf <- function(vcf_file) {
  
  cat("Reading VCF file:", vcf_file, "\n")
  
  # Detect if file is gzipped
  if(grepl("\\.gz$", vcf_file)) {
    con <- gzfile(vcf_file, "r")
  } else {
    con <- file(vcf_file, "r")
  }
  
  # Read header lines
  header_lines <- c()
  sample_line <- NULL
  
  while(TRUE) {
    line <- readLines(con, n = 1)
    if(length(line) == 0) break
    
    if(grepl("^##", line)) {
      header_lines <- c(header_lines, line)
    } else if(grepl("^#CHROM", line)) {
      sample_line <- line
      break
    }
  }
  
  # Parse sample names
  fields <- strsplit(sample_line, "\t")[[1]]
  samples <- fields[10:length(fields)]
  cat("  Samples found:", length(samples), "\n")
  
  # Read variant data
  variant_data <- read.table(con, sep = "\t", stringsAsFactors = FALSE,
                            col.names = c("CHROM", "POS", "ID", "REF", "ALT", 
                                        "QUAL", "FILTER", "INFO", "FORMAT", samples),
                            comment.char = "")
  close(con)
  
  cat("  Variants loaded:", nrow(variant_data), "\n")
  
  # Extract genotype matrix
  gt_cols <- 10:ncol(variant_data)
  gt_matrix <- as.matrix(variant_data[, gt_cols])
  rownames(gt_matrix) <- paste(variant_data$CHROM, variant_data$POS, sep = "_")
  
  # Parse genotypes (extract GT field)
  gt_matrix <- apply(gt_matrix, 2, function(x) {
    sapply(strsplit(as.character(x), ":"), function(y) y[1])
  })
  
  result <- list(
    header = header_lines,
    variants = variant_data[, 1:9],
    genotypes = gt_matrix,
    samples = samples
  )
  
  return(result)
}


#' Convert genotypes to numeric coding
#' @param gt_matrix Genotype matrix from VCF
#' @return Numeric genotype matrix (0=Ref, 1=Het, 2=Alt, NA=Missing)
encode_genotypes <- function(gt_matrix) {
  
  cat("Encoding genotypes...\n")
  
  # Convert to numeric: 0/0 -> 0, 0/1 or 1/0 -> 1, 1/1 -> 2, ./. -> NA
  numeric_gt <- apply(gt_matrix, 2, function(x) {
    sapply(x, function(g) {
      g <- as.character(g)
      if(grepl("\\./\\.", g) || grepl("^\\.$", g) || is.na(g)) {
        return(NA)
      } else if(g == "0/0" || g == "0|0") {
        return(0)
      } else if(g == "0/1" || g == "1/0" || g == "0|1" || g == "1|0") {
        return(1)
      } else if(g == "1/1" || g == "1|1") {
        return(2)
      } else if(g == "2/2" || g == "2|2") {
        # Multiallelic homozygous alt
        return(2)
      } else if(grepl("[0-9]/[0-9]", g) || grepl("[0-9]\\|[0-9]", g)) {
        # Other multiallelic variants - treat as missing with warning
        if(!exists("multiallelic_warned")) {
          warning("Multiallelic variants detected. Treating complex genotypes as missing data.")
          assign("multiallelic_warned", TRUE, envir = .GlobalEnv)
        }
        return(NA)
      } else {
        return(NA)  # Unknown format
      }
    })
  })
  
  # Clean up warning flag
  if(exists("multiallelic_warned")) {
    rm(multiallelic_warned, envir = .GlobalEnv)
  }
  
  rownames(numeric_gt) <- rownames(gt_matrix)
  colnames(numeric_gt) <- colnames(gt_matrix)
  
  cat("  Genotypes encoded\n")
  cat("  Missing data:", sum(is.na(numeric_gt)), "calls\n")
  cat("  Heterozygous:", sum(numeric_gt == 1, na.rm = TRUE), "calls\n")
  
  return(numeric_gt)
}


#' Filter samples based on quality metrics
#' @param gt_matrix Numeric genotype matrix
#' @param max_heterozygosity Maximum proportion of heterozygous calls allowed (default 0)
#' @param max_missing Maximum proportion of missing data allowed (default 0.1)
#' @return List with cleaned matrix and statistics
filter_samples <- function(gt_matrix, max_heterozygosity = 0, max_missing = 0.1) {
  
  cat("\nFiltering samples...\n")
  
  n_variants <- nrow(gt_matrix)
  n_samples <- ncol(gt_matrix)
  
  # Calculate statistics per sample
  sample_stats <- data.frame(
    Sample = colnames(gt_matrix),
    N_heterozygous = colSums(gt_matrix == 1, na.rm = TRUE),
    N_missing = colSums(is.na(gt_matrix)),
    Prop_heterozygous = colSums(gt_matrix == 1, na.rm = TRUE) / n_variants,
    Prop_missing = colSums(is.na(gt_matrix)) / n_variants,
    stringsAsFactors = FALSE
  )
  
  # Identify samples to remove
  to_remove <- sample_stats$Prop_heterozygous > max_heterozygosity | 
               sample_stats$Prop_missing > max_missing
  
  removed_samples <- sample_stats$Sample[to_remove]
  kept_samples <- sample_stats$Sample[!to_remove]
  
  cat("  Initial samples:", n_samples, "\n")
  cat("  Removed:", length(removed_samples), "\n")
  cat("    - Heterozygous:", sum(sample_stats$Prop_heterozygous[to_remove] > max_heterozygosity), "\n")
  cat("    - Missing data:", sum(sample_stats$Prop_missing[to_remove] > max_missing), "\n")
  cat("  Retained:", length(kept_samples), "\n\n")
  
  # Filter matrix
  filtered_gt <- gt_matrix[, !to_remove, drop = FALSE]
  
  return(list(
    genotypes = filtered_gt,
    kept_samples = kept_samples,
    removed_samples = removed_samples,
    sample_stats = sample_stats
  ))
}


#' Call haplotypes based on unique genotype patterns
#' @param gt_matrix Numeric genotype matrix (variants x samples)
#' @param min_samples Minimum number of samples required for a haplotype (default 1)
#' @param hap_prefix Prefix for haplotype names (default "H")
#' @return List of haplotypes with assignments and patterns
call_haplotypes <- function(gt_matrix, min_samples = 1, hap_prefix = "H") {
  
  cat("Calling haplotypes...\n")
  
  # Transpose to samples x variants
  gt_transpose <- t(gt_matrix)
  
  # Convert each row to a string pattern, handling NA
  haplotype_strings <- apply(gt_transpose, 1, function(x) {
    if(any(is.na(x))) {
      return(NA_character_)
    }
    paste(x, collapse = "_")
  })
  
  # Remove samples with missing data
  samples_with_na <- is.na(haplotype_strings)
  if(sum(samples_with_na) > 0) {
    cat("  Excluding", sum(samples_with_na), "samples with missing genotype data\n")
    haplotype_strings <- haplotype_strings[!samples_with_na]
  }
  
  # Find unique patterns
  unique_patterns <- unique(haplotype_strings)
  cat("  Unique patterns found:", length(unique_patterns), "\n")
  
  # Group samples by haplotype
  haplotype_list <- list()
  variant_patterns <- list()
  haplotype_assignments <- character(length(haplotype_strings))
  names(haplotype_assignments) <- names(haplotype_strings)
  
  hap_counter <- 1
  for(i in seq_along(unique_patterns)) {
    pattern <- unique_patterns[i]
    members <- names(haplotype_strings)[haplotype_strings == pattern]
    
    if(length(members) >= min_samples) {
      hap_name <- paste0(hap_prefix, hap_counter)
      haplotype_list[[hap_name]] <- members
      
      # Store variant pattern
      variant_pattern <- as.numeric(strsplit(pattern, "_")[[1]])
      variant_patterns[[hap_name]] <- variant_pattern
      
      # Store assignments
      haplotype_assignments[members] <- hap_name
      
      hap_counter <- hap_counter + 1
    }
  }
  
  cat("  Haplotypes called:", length(haplotype_list), "\n")
  
  # Print haplotype sizes
  for(hap in names(haplotype_list)) {
    cat("   ", hap, ":", length(haplotype_list[[hap]]), "samples\n")
  }
  
  return(list(
    haplotypes = haplotype_list,
    variant_patterns = variant_patterns,
    haplotype_assignments = haplotype_assignments[haplotype_assignments != ""]
  ))
}


#' Create accession-haplotype data frame
#' @param haplotype_result Result from call_haplotypes()
#' @return Data frame with Accession and Haplotype columns
create_haplotype_table <- function(haplotype_result) {
  
  assignments <- haplotype_result$haplotype_assignments
  
  hap_df <- data.frame(
    Accession = names(assignments),
    Haplotype = as.character(assignments),
    stringsAsFactors = FALSE
  )
  
  return(hap_df)
}

cat("\n")
cat("================================================================================\n")
cat("  SuperHap Core Functions Loaded\n")
cat("================================================================================\n")
cat("Available functions:\n")
cat("  - read_vcf()           : Parse VCF files\n")
cat("  - encode_genotypes()   : Convert to numeric\n")
cat("  - filter_samples()     : Quality control\n")
cat("  - call_haplotypes()    : Identify haplotypes\n")
cat("  - create_haplotype_table() : Create accession mapping\n")
cat("================================================================================\n\n")
