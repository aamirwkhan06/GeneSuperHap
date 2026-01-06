# GeneSuperHap: Superior Haplotype Analysis for Plant Breeding

**A breeding-focused haplotype analysis toolkit with machine learning capabilities**

## Overview

GeneSuperHap is an independent R-based pipeline for haplotype analysis specifically designed for plant breeding programs. Unlike existing tools that focus primarily on descriptive analysis, GeneSuperHap integrates **machine learning** and **multi-trait optimization** to provide actionable breeding recommendations.

### Key Differentiator
**GeneSuperHap not only identifies haplotypes - it predicts which ones will perform best and tells you which to use in your breeding program.**

---

## Why GeneSuperHap?

- ✅ GeneSuperHap: Predictive ML-based ranking + breeding decisions
- ✅ GeneSuperHap: Multi-trait selection indices

---

## Features

### 1. **Machine Learning Haplotype Ranking** 🤖

Train Random Forest models to predict haplotype performance based on phenotypic data.

```r
# Train ML model
ml_model <- train_haplotype_ranker(
  haplotype_data = hap_table,
  phenotype_data = phenotypes,
  target_trait = "Yield"
)

# View top-ranked haplotypes
print(ml_model$haplotype_rankings)
```

**What it does:**
- Trains predictive models on your data
- Ranks haplotypes by predicted performance
- Provides cross-validation metrics
- Quantifies prediction uncertainty

**Why it matters:**
- Prioritize germplasm for phenotyping
- Reduce expensive field trials
- Identify hidden patterns

---

### 2. **Multi-Trait Selection Index** 📊

Combine multiple traits into a single breeding criterion with customizable weights.

```r
# Calculate selection index
selection_idx <- calculate_selection_index(
  haplotype_data = hap_table,
  phenotype_data = phenotypes,
  trait_weights = c(Yield = 2, Protein = 1.5, Oil = 1),
  trait_direction = c(Yield = 1, Protein = 1, Oil = -1)  # Minimize oil
)
```

**What it does:**
- Weighted combination of traits
- Standardized z-scores for fair comparison
- Maximize or minimize each trait
- Breeding-ready rankings

**Why it matters:**
- Real breeding involves multiple objectives
- Single-trait selection often fails
- Accounts for trait importance

---

### 3. **Haplotype Stability Analysis** 🌍

Evaluate performance across environments to identify robust haplotypes.

```r
# Analyze stability
stability <- calculate_haplotype_stability(
  haplotype_data = hap_table,
  phenotype_data = multi_env_phenotypes,  # Must have 'Environment' column
  trait = "Yield"
)
```

**What it does:**
- Coefficient of variation across environments
- Identifies stable performers
- GxE interaction assessment

**Why it matters:**
- Environmental variability is reality
- Wide adaptation is valuable
- Reduces deployment risk

---

### 4. **Network-Based Haplotype Relationships** 🕸️

Visualize genetic distances and relationships between haplotypes.

```r
# Calculate distances
distances <- calculate_haplotype_distances(
  variant_patterns = haplotypes$variant_patterns
)

# Use for network visualization (see visualization module)
```

**What it does:**
- Hamming distance matrix
- Quantifies genetic diversity
- Input for network plots

**Why it matters:**
- Understand haplotype relationships
- Guide diversity management
- Aesthetic visualizations

---

### 5. **Enhanced Customization** 🎨

Full control over haplotype names, colors, and plot aesthetics.

```r
# Rename haplotypes
haplotypes <- rename_haplotypes(
  haplotypes,
  new_names = c("Elite1", "Elite2", "Wild", "Landrace")
)

# Custom color scheme
colors <- create_haplotype_colors(
  n_haplotypes = 4,
  custom_colors = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00")
)
```

**What it does:**
- Custom haplotype naming
- Flexible color palettes
- Publication-quality plots

**Why it matters:**
- Professional presentations
- Journal requirements
- Branding consistency

---

## Installation

### Prerequisites

```r
# Required packages
install.packages(c("dplyr", "tidyr", "ggplot2", "RColorBrewer"))

# Machine learning
install.packages(c("randomForest", "caret"))

# Visualization (optional)
install.packages(c("igraph", "pheatmap", "ggrepel"))
```

### Download

```bash
git clone https://github.com/aamirwkhan06/GeneSuperHap.git
cd GeneSuperHap
```

---

## Quick Start

```r
# Load GeneSuperHap
source("core_functions.R")
source("unique_features.R")

# 1. Read VCF and call haplotypes
vcf_data <- read_vcf("variants.vcf.gz")
gt_numeric <- encode_genotypes(vcf_data$genotypes)
filtered <- filter_samples(gt_numeric, max_heterozygosity = 0, max_missing = 0.1)
haplotypes <- call_haplotypes(filtered$genotypes, min_samples = 3)

# 2. Create haplotype table
hap_table <- create_haplotype_table(haplotypes)

# 3. Load phenotypes
phenotypes <- read.table("phenotypes.txt", header = TRUE, sep = "\t")

# 4. Train ML model
ml_model <- train_haplotype_ranker(hap_table, phenotypes, "Yield")

# 5. Calculate selection index
index <- calculate_selection_index(
  hap_table, phenotypes,
  trait_weights = c(Yield = 2, Protein = 1),
  trait_direction = c(Yield = 1, Protein = 1)
)

# 6. Check stability (if multi-environment data)
stability <- calculate_haplotype_stability(hap_table, phenotypes, "Yield")
```

---

## Input Files

### 1. VCF File
Standard VCF format (v4.x), can be gzipped

```
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
Chr01   1000  .   A    G    100   PASS    .     GT      0/0      1/1
```

### 2. Phenotype Data
Tab-delimited with Accession column

```
Accession  Yield  Protein  Oil
Sample1    3500   42.5     18.2
Sample2    3200   40.1     19.5
```

### 3. Multi-Environment Data (for stability)
Include 'Environment' column

```
Accession  Yield  Environment
Sample1    3500   Location1
Sample1    3300   Location2
Sample2    3200   Location1
```

---

## Complete Workflow Example

```r
# ==============================================================================
# Complete GeneSuperHap Analysis
# ==============================================================================

# Load functions
source("core_functions.R")
source("unique_features.R")
source("visualization.R")  # For plots

# Read and process VCF
vcf <- read_vcf("soybean_protein_gene.vcf.gz")
gt_numeric <- encode_genotypes(vcf$genotypes)
filtered <- filter_samples(gt_numeric)
haplotypes <- call_haplotypes(filtered$genotypes)
hap_table <- create_haplotype_table(haplotypes)

# Load phenotypes
pheno <- read.table("protein_phenotypes.txt", header = TRUE, sep = "\t")

# Feature 1: Machine Learning Ranking
ml_results <- train_haplotype_ranker(hap_table, pheno, "Protein_Content")

# Feature 2: Multi-Trait Selection
selection <- calculate_selection_index(
  hap_table, pheno,
  trait_weights = c(Protein_Content = 2, Oil_Content = 1),
  trait_direction = c(Protein_Content = 1, Oil_Content = 1)
)

# Feature 3: Stability Analysis (if multi-env data)
if("Environment" %in% colnames(pheno)) {
  stability <- calculate_haplotype_stability(hap_table, pheno, "Protein_Content")
}

# Feature 4: Genetic Distances
distances <- calculate_haplotype_distances(haplotypes$variant_patterns)

# Feature 5: Custom naming and colors
haplotypes_renamed <- rename_haplotypes(haplotypes, c("High", "Medium", "Low"))
colors <- create_haplotype_colors(3, custom_colors = c("#FF0000", "#00FF00", "#0000FF"))

# Save results
write.table(ml_results$haplotype_rankings, "ml_rankings.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(selection$summary, "selection_index.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

---

## Applications

### 1. Breeding Programs
- Rank haplotypes for crossing
- Multi-trait selection decisions
- Parent selection

### 2. Germplasm Management
- Prioritize accessions for phenotyping
- Assess diversity
- Guide collections

### 3. QTL Follow-up
- Fine-map causal haplotypes
- Test across environments
- Deployment strategies

---

## Citation

If you use GeneSuperHap, please cite:

```
GitHub: https://github.com/aamirwkhan06/GeneSuperHap
```

---

## Contact

**Aamir W. Khan, PhD**  
Research Scientist, Bioinformatics  

---

## License

MIT License - Free for research and commercial use

---

## Future Development

Planned features:
- Additional ML algorithms (XGBoost, neural networks)
- Web interface
- More visualization types
- Integration with genomic selection

---

**GeneSuperHap: From Haplotypes to Breeding Decisions**
