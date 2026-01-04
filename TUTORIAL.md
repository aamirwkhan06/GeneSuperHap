# SuperHap: Complete Tutorial

## Step-by-Step Guide to Haplotype Analysis with Machine Learning

---

## Table of Contents
1. [Setup and Installation](#setup)
2. [Data Preparation](#data-prep)
3. [Running the Pipeline](#pipeline)
4. [Understanding the Results](#results)
5. [Customization](#customization)
6. [Troubleshooting](#troubleshooting)

---

## 1. Setup and Installation

### 1.1 Install R Packages

Open R or RStudio and run:

```r
# Core packages
install.packages(c("dplyr", "tidyr", "ggplot2", "RColorBrewer"))

# Machine learning
install.packages(c("randomForest", "caret"))

# Visualization
install.packages(c("pheatmap", "igraph", "ggrepel", "gridExtra"))
```

### 1.2 Download SuperHap

```bash
git clone https://github.com/aamirwkhan06/GeneSuperHap.git
cd GeneSuperHap
```

### 1.3 Verify Installation

```r
# Test if packages load
library(dplyr)
library(randomForest)
library(ggplot2)

# Load SuperHap functions
source("core_functions.R")
source("unique_features.R")
source("visualization.R")

# Should see confirmation messages
```

---

## 2. Data Preparation

### 2.1 Prepare Your VCF File

**Requirements:**
- Standard VCF format (v4.x)
- Can be gzipped (.vcf.gz) or uncompressed (.vcf)
- Sample names in header must match phenotype file

**Example VCF structure:**
```
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
Chr01   1000  .   A    G    100   PASS    .     GT      0/0      1/1
Chr01   2000  .   C    T    100   PASS    .     GT      0/0      0/0
```

**Place your VCF in:**
```
input_data/variants.vcf.gz
```

### 2.2 Prepare Phenotype File

**Requirements:**
- Tab-delimited text file
- Must have "Accession" column matching VCF sample names
- Numeric trait values
- Optional: "Environment" column for stability analysis

**Example phenotypes.txt:**
```
Accession  Yield  Protein  Oil  Environment
Sample1    3500   42.5     18.2  Location1
Sample2    3200   40.1     19.5  Location1
```

**Place your phenotypes in:**
```
input_data/phenotypes.txt
```

### 2.3 Directory Structure

Your project should look like:
```
superhap/
├── core_functions.R
├── unique_features.R
├── visualization.R
├── 01_haplotype_calling.R
├── 02_superior_haplotypes.R
├── 03_ml_analysis.R
├── 04_selection_index.R
├── 05_network_analysis.R
├── input_data/
│   ├── variants.vcf.gz
│   └── phenotypes.txt
└── results/  (will be created automatically)
```

---

## 3. Running the Pipeline

### 3.1 Quick Start (Run Everything)

```r
# Run all steps sequentially
source("01_haplotype_calling.R")
source("02_superior_haplotypes.R")
source("03_ml_analysis.R")
source("04_selection_index.R")
source("05_network_analysis.R")
```

### 3.2 Step-by-Step Execution

#### Step 1: Haplotype Calling

```r
source("01_haplotype_calling.R")
```

**What it does:**
- Reads your VCF file
- Removes samples with heterozygous calls or excessive missing data
- Identifies unique haplotypes
- Creates accession-haplotype mapping

**Output files:**
- `results/haplotypes_accessions.txt` - Main mapping file
- `results/haplotypes_summary.txt` - Frequency table
- `results/haplotypes_data.RData` - R object for next steps

**Customization in script:**
```r
MAX_HETEROZYGOSITY <- 0      # Change if you want to allow some heterozygosity
MAX_MISSING <- 0.1            # Allow up to 10% missing data
MIN_SAMPLES_PER_HAP <- 3      # Minimum samples to call a haplotype
```

---

#### Step 2: Superior Haplotype Identification

```r
source("02_superior_haplotypes.R")
```

**What it does:**
- Statistical comparisons between haplotypes for each trait
- Identifies superior haplotype per trait
- Performs pairwise t-tests

**Output files:**
- `results/superior_haplotypes_summary.txt` - Overall summary
- `results/superior_haplotypes_[TRAIT]_stats.txt` - Per-trait statistics
- `results/superior_haplotypes_[TRAIT]_pairwise.txt` - Pairwise comparisons

**Look for:**
- Which haplotype has highest mean for each trait
- Significant differences (P < 0.05)

---

#### Step 3: Machine Learning Analysis

```r
source("03_ml_analysis.R")
```

**What it does:**
- Trains Random Forest model
- Predicts haplotype performance
- Ranks haplotypes by ML prediction
- Cross-validates results

**Output files:**
- `results/ml_analysis_rankings.txt` - ML-based rankings
- `results/ml_analysis_predictions.pdf` - Observed vs Predicted plot
- `results/ml_analysis_model.RData` - Trained model

**Key metrics to check:**
- **R-squared**: How well model predicts (>0.5 is good)
- **RMSE**: Prediction error (lower is better)

**Customization in script:**
```r
TARGET_TRAIT <- "Yield"       # Change to your trait
TEST_FRACTION <- 0.2          # 20% for testing
N_TREES <- 500                # More trees = better but slower
```

---

#### Step 4: Multi-Trait Selection Index

```r
source("04_selection_index.R")
```

**What it does:**
- Combines multiple traits into single breeding criterion
- Applies custom weights to each trait
- Identifies best multi-trait performers

**Output files:**
- `results/selection_index_summary.txt` - Rankings
- `results/selection_index_rankings.pdf` - Visualization
- `results/stability_analysis.txt` - If Environment data available

**IMPORTANT - Customize for your goals:**
```r
TRAIT_WEIGHTS <- c(
  Yield = 2.0,      # Highest priority
  Protein = 1.5,    
  Oil = 1.0         
)

TRAIT_DIRECTION <- c(
  Yield = 1,        # 1 = maximize
  Protein = 1,      
  Oil = -1          # -1 = minimize
)
```

---

#### Step 5: Network Analysis

```r
source("05_network_analysis.R")
```

**What it does:**
- Calculates genetic distances between haplotypes
- Creates network visualization
- Distribution plots

**Output files:**
- `results/network_graph.pdf` - Network diagram
- `results/network_heatmap.pdf` - Distance heatmap
- `results/distribution_barplot.pdf` - Frequency bar chart
- `results/ANALYSIS_SUMMARY.txt` - Complete summary

---

## 4. Understanding the Results

### 4.1 Which Haplotype is Best?

**Check these files in order:**

1. **ML Rankings** (`ml_analysis_rankings.txt`)
   - Look at `ML_Rank` column
   - Lower rank = better predicted performance

2. **Selection Index** (`selection_index_summary.txt`)
   - Look at `Mean_Index` column
   - Higher value = better overall

3. **Stability** (`stability_analysis.txt`) - if available
   - Look at `CV` column
   - Lower CV = more stable

**Example decision:**
```
Best haplotype for breeding:
- High ML_Rank (top 3)
- High Selection_Index (top 5)
- Low CV (stable across environments)
```

### 4.2 Interpreting Plots

**ML Predictions Plot:**
- Points near diagonal line = good predictions
- R² close to 1 = excellent model
- Groups of points = haplotype clusters

**Selection Index Plot:**
- Green bars = positive index (favorable)
- Red bars = negative index (unfavorable)
- Height = breeding value

**Network Graph:**
- Node size = number of samples
- Edges = genetic similarity
- Clusters = related haplotypes

**Stability Plot:**
- Left side = high performance
- Bottom = low variability (stable)
- Ideal = bottom-left corner

---

## 5. Customization

### 5.1 Custom Haplotype Names

Instead of H1, H2, H3, use meaningful names:

```r
# Load your haplotype data
load("results/haplotypes_data.RData")

# Rename based on your knowledge
haplotypes <- rename_haplotypes(
  haplotypes,
  new_names = c("Elite", "Wild", "Landrace")
)

# Update table
hap_table <- create_haplotype_table(haplotypes)
```

### 5.2 Custom Colors

```r
# Define your colors
my_colors <- c(
  "Elite" = "#FF0000",     # Red
  "Wild" = "#00FF00",      # Green
  "Landrace" = "#0000FF"   # Blue
)

# Use in plots
plot_haplotype_distribution(hap_table, colors = my_colors)
```

### 5.3 Different ML Models

For advanced users:

```r
# Modify unique_features.R
# In train_haplotype_ranker() function, try:

# More trees for better accuracy
rf_model <- randomForest(model_formula, data = train_data,
                        ntree = 1000)  # Instead of 500

# Different number of features per split
rf_model <- randomForest(model_formula, data = train_data,
                        mtry = 3)      # Adjust based on features
```

---

## 6. Troubleshooting

### Problem 1: "File not found"

**Solution:**
```r
# Check current directory
getwd()

# Should be in superhap folder
# If not:
setwd("/path/to/superhap")
```

### Problem 2: "Package not found"

**Solution:**
```r
# Reinstall package
install.packages("package_name")

# If that fails, try:
install.packages("package_name", dependencies = TRUE)
```

### Problem 3: "Accession mismatch"

**Error:** "No samples match between VCF and phenotypes"

**Solution:**
```r
# Check VCF samples
vcf <- read_vcf("input_data/variants.vcf.gz")
head(vcf$samples)

# Check phenotype samples
pheno <- read.table("input_data/phenotypes.txt", header = TRUE)
head(pheno$Accession)

# They must match exactly!
```

### Problem 4: "Insufficient data for ML"

**Error:** "Small sample size. Results may be unreliable"

**Solution:**
- You need at least 20 samples for reliable ML
- If you have fewer, skip ML analysis
- Or reduce test_fraction: `TEST_FRACTION <- 0.1`

### Problem 5: "All samples excluded"

**Solution:**
```r
# Relax filtering criteria in 01_haplotype_calling.R
MAX_HETEROZYGOSITY <- 0.05    # Allow 5% heterozygosity
MAX_MISSING <- 0.2             # Allow 20% missing
```

### Problem 6: ML model R² is very low

**Possible causes:**
1. Haplotype doesn't affect the trait
2. Environment effects dominate
3. Not enough variants in region

**Solution:**
- Check if you're analyzing the right gene/region
- Try different traits
- Add more variants to VCF if possible

---

## 7. Example Workflow

### Complete Example from Start to Finish

```r
# ==============================================================================
# COMPLETE GENESUPERHAP ANALYSIS EXAMPLE
# ==============================================================================

# 1. Setup
setwd("/path/to/superhap")
source("core_functions.R")
source("unique_features.R")
source("visualization.R")

# 2. Run pipeline
cat("Starting SuperHap pipeline...\n")

# Step 1: Call haplotypes
source("01_haplotype_calling.R")
# Output: 5 haplotypes identified

# Step 2: Statistical analysis
source("02_superior_haplotypes.R")
# Output: H3 superior for Yield, H1 for Protein

# Step 3: Machine Learning
source("03_ml_analysis.R")
# Output: R² = 0.72, H3 ranked #1

# Step 4: Selection index
source("04_selection_index.R")
# Output: H3 best overall (combined traits)

# Step 5: Network & visualization
source("05_network_analysis.R")
# Output: Beautiful plots generated

# 3. Review results
cat("\nReview these files:\n")
cat("  results/ml_analysis_rankings.txt\n")
cat("  results/selection_index_summary.txt\n")
cat("  results/ANALYSIS_SUMMARY.txt\n")

# 4. Breeding decision
cat("\nBREEDING RECOMMENDATION:\n")
cat("  Use H3 as parent in crosses\n")
cat("  High yield + stable + good protein\n")
```

---

## 8. Next Steps

### After Running the Pipeline

1. **Validate top haplotypes**
   - Check in field trials
   - Verify with additional markers

2. **Use in breeding**
   - Select accessions with top haplotypes
   - Design crosses
   - Track haplotype segregation

3. **Expand analysis**
   - Add more traits
   - Include more environments
   - Analyze additional genes

### Advanced Analyses

```r
# Cross-validation for ML model
source("unique_features.R")

# Run 5-fold cross-validation
cv_results <- cross_validate_haplotype_model(
  hap_table, phenotypes,
  target_trait = "Yield",
  k_folds = 5
)
```

---

## 9. Citation

If you use GeneSuperHap in your research:

```
GeneSuperHap: A machine learning-enhanced haplotype analysis pipeline for plant breeding. 
GitHub: [https://github.com/aamirwkhan06/GeneSuperHap]
```

---

## 10. Getting Help

- **GitHub Issues:** Report bugs or request features
- **Email:** maky74@missouri.edu
- **Documentation:** Check README.md

---

**You're all set! Happy haplotyping! 🧬🌱**
