# Compartment Covariate Analysis

## Purpose
Explore the effect of adding compartment covariates (tumor/stroma) to SHADE models to control for spatial confounding. The simulation study showed SHADE can incorrectly attribute compartment effects to cell-cell interactions when both cell types are abundant. This analysis investigates whether adding compartment covariates mitigates this issue.

## Biological Motivation
Tissue compartmentalization (tumor islands vs stroma) is ubiquitous in solid tumors. Cell types may appear to interact simply because they co-localize within the same compartment. By adding compartment indicators as covariates, we can ask: **"Within the same tissue compartment, do cells truly interact?"**

## Workflow

### 01_identify_compartments_binary.R
- **Input**: CRC dataset images
- **Output**: Binary compartment fields (tumor-enriched vs stromal) saved to `data/binary_compartment_fields.rds`
- **Methods**:
  - Define compartments based on local tumor cell density (>75th percentile = tumor-enriched)
  - Create spatial fields for compartment lookup at any location

### 02_fit_models_with_compartments.R
- **Input**: CRC data + compartment fields
- **Output**: SHADE model fits with compartment covariates for CTLs, Memory CD4+ T, and Granulocytes
- **Methods**:
  - Add binary compartment covariate to SHADE logistic regression
  - Fit models for each target cell type with truncated RBFs (5 basis functions, sigma=12μm, min=30μm)
  - Save fits to `data/fit_{target}_with_compartment.rds` and metadata to `data/metadata_{target}_with_compartment.rds`

### 03_analyze_compartment_models.R (Exploratory)
- **Input**: Compartment-adjusted model fits
- **Output**: Exploratory visualizations in `compartment_analysis_figures/` (not used in manuscript)
- **Methods**:
  - Extract and visualize SICs from compartment-adjusted models
  - Generate exploratory plots for CLR vs DII comparisons

### 03_compare_with_without.R
- **Input**: Original model fits (from crc_analysis/) + compartment-adjusted fits
- **Output**: Comparison data saved to `data/comparison_sics.rds`
- **Methods**:
  - Extract SICs from both original and compartment-adjusted models
  - Compute differences in estimated interactions
  - Identify which interactions change substantially when controlling for compartments
  - Used to generate Supplement Figure (compartment comparison)

### 04_rbf_sensitivity_analysis.R (Exploratory)
- **Input**: Compartment-adjusted model fits
- **Output**: Sensitivity analysis results in `rbf_sensitivity/` (not used in manuscript)
- **Methods**:
  - Test robustness of results to RBF hyperparameter choices
  - Exploratory analysis for methodological validation

### 05_generate_paper_figures.R
- **Input**: Compartment-adjusted model fits and comparison data
- **Output**: Manuscript-ready figures and tables in `manuscript_figures/`
- **Figures**:
  - `compartment_effects_table.csv` - Compartment covariate effects (for main text table)
  - `compartment_comparison_supplement.pdf` - Original vs adjusted SIC comparison (for supplement)
  - `sic_change_summary.csv` and `clr_dii_concordance.csv` - Summary statistics

## Expected Insights

1. **Which interactions are confounded by compartmentalization?**
   - Large changes when adjusting = compartment-driven co-localization
   - Small changes = true cell-cell interaction

2. **Do biological conclusions change?**
   - Are CLR vs DII differences explained by compartment structure rather than cell interactions?

3. **Guidance for future analyses**
   - Recommend always checking for compartment effects
   - Provide diagnostics for when compartment adjustment is necessary

## Usage

```r
# Run analysis pipeline for manuscript
Sys.setenv(SYSTEM_ENV="laptop")
source("compartment_analysis/01_identify_compartments_binary.R")
source("compartment_analysis/02_fit_models_with_compartments.R")
source("compartment_analysis/03_compare_with_without.R")
source("compartment_analysis/05_generate_paper_figures.R")

# Optional exploratory analyses (not used in manuscript)
source("compartment_analysis/03_analyze_compartment_models.R")
source("compartment_analysis/04_rbf_sensitivity_analysis.R")
```

## Notes

- Script 01 defines compartments based on local tumor cell density (>75th percentile threshold)
- Scripts 03_analyze and 04_rbf are exploratory and not used in the final manuscript
- The main pipeline (01 → 02 → 03_compare → 05) generates all figures/tables for the paper
- Compartment adjustment uses the same truncated RBF basis as the main analysis (5 basis functions, sigma=12μm, min=30μm)
