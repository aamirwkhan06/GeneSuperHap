#!/usr/bin/env Rscript
# ==============================================================================
# SuperHap Pipeline - Step 01: Haplotype Calling
# Author: Aamir Khan
# Description: Read VCF, filter samples, and call haplotypes
# ==============================================================================

# Load required functions
source("core_functions.R")

cat("\n")
cat("================================================================================\n")
cat("  SUPERHAP PIPELINE - STEP 01: HAPLOTYPE CALLING\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

VCF_FILE <- "input_data/variants.vcf.gz"
OUTPUT_DIR <- "results"
OUTPUT_PREFIX <- "haplotypes"

# Quality control parameters
MAX_HETEROZYGOSITY <- 0      # Remove samples with any heterozygous calls
MAX_MISSING <- 0.1            # Remove samples with >10% missing data
MIN_SAMPLES_PER_HAP <- 3      # Minimum samples to call a haplotype
HAPLOTYPE_PREFIX <- "H"       # Prefix for haplotype names (H1, H2, ...)

# ==============================================================================
# VALIDATE INPUT FILES
# ==============================================================================

if(!file.exists(VCF_FILE)) {
  VCF_FILE_ALT <- gsub("\\.gz$", "", VCF_FILE)
  if(file.exists(VCF_FILE_ALT)) {
    VCF_FILE <- VCF_FILE_ALT
  } else {
    stop("VCF file not found: ", VCF_FILE, 
         "\nPlease ensure file exists or update VCF_FILE path")
  }
}

# ==============================================================================
# CREATE OUTPUT DIRECTORY
# ==============================================================================

if(!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR, "\n\n")
}

# ==============================================================================
# STEP 1: READ VCF FILE
# ==============================================================================

cat("STEP 1: Reading VCF file\n")
cat("----------------------------------------\n")

vcf_data <- read_vcf(VCF_FILE)

cat("\nVCF Summary:\n")
cat("  Chromosome(s):", unique(vcf_data$variants$CHROM), "\n")
cat("  Position range:", min(vcf_data$variants$POS), "-", max(vcf_data$variants$POS), "\n")
cat("  Total variants:", nrow(vcf_data$variants), "\n")
cat("  Total samples:", length(vcf_data$samples), "\n\n")

# ==============================================================================
# STEP 2: ENCODE GENOTYPES
# ==============================================================================

cat("STEP 2: Encoding genotypes to numeric format\n")
cat("----------------------------------------\n")

gt_numeric <- encode_genotypes(vcf_data$genotypes)

# ==============================================================================
# STEP 3: FILTER SAMPLES
# ==============================================================================

cat("STEP 3: Quality control - filtering samples\n")
cat("----------------------------------------\n")

filtered <- filter_samples(
  gt_numeric,
  max_heterozygosity = MAX_HETEROZYGOSITY,
  max_missing = MAX_MISSING
)

# Save exclusion report
if(length(filtered$removed_samples) > 0) {
  exclusion_report <- filtered$sample_stats %>%
    filter(Sample %in% filtered$removed_samples) %>%
    mutate(
      Reason = case_when(
        Prop_heterozygous > MAX_HETEROZYGOSITY ~ "Heterozygous",
        Prop_missing > MAX_MISSING ~ "Missing_data",
        TRUE ~ "Both"
      )
    )
  
  write.table(
    exclusion_report,
    file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_excluded_samples.txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("Exclusion report saved\n\n")
}

# ==============================================================================
# STEP 4: CALL HAPLOTYPES
# ==============================================================================

cat("STEP 4: Calling haplotypes\n")
cat("----------------------------------------\n")

haplotypes <- call_haplotypes(
  filtered$genotypes,
  min_samples = MIN_SAMPLES_PER_HAP,
  hap_prefix = HAPLOTYPE_PREFIX
)

# ==============================================================================
# STEP 5: CREATE OUTPUT FILES
# ==============================================================================

cat("\nSTEP 5: Generating output files\n")
cat("----------------------------------------\n")

# 5.1 Haplotype-Accession mapping table
hap_table <- create_haplotype_table(haplotypes)

write.table(
  hap_table,
  file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_accessions.txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("  Saved:", file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_accessions.txt")), "\n")

# 5.2 Haplotype summary
hap_summary <- hap_table %>%
  group_by(Haplotype) %>%
  summarise(
    N_samples = n(),
    Percentage = round(n() / nrow(hap_table) * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(N_samples))

write.table(
  hap_summary,
  file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_summary.txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("  Saved:", file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_summary.txt")), "\n")

# 5.3 Save R object for later use
save(
  haplotypes,
  hap_table,
  hap_summary,
  filtered,
  vcf_data,
  file = file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_data.RData"))
)

cat("  Saved:", file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_data.RData")), "\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STEP 01 COMPLETE\n")
cat("================================================================================\n")
cat("\nSummary:\n")
cat("  Input samples:", length(vcf_data$samples), "\n")
cat("  Excluded samples:", length(filtered$removed_samples), "\n")
cat("  Retained samples:", length(filtered$kept_samples), "\n")
cat("  Haplotypes identified:", nrow(hap_summary), "\n")
cat("\nTop 5 haplotypes:\n")
print(head(hap_summary, 5))
cat("\nNext step: Run 02_superior_haplotypes.R\n")
cat("================================================================================\n\n")
