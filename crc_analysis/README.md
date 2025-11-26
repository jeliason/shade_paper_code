# Colorectal Cancer (CRC) Analysis

Analysis of colorectal cancer tissue data from Schürch et al. (2020) using the SHADE model.

## Purpose

Uncover how immune and tumor cells are spatially organized in the tumor microenvironment:

- Model directional spatial interactions between key cell types
- Identify hierarchical spatial patterns across images, patients, and cohorts (CLR vs DII)
- Characterize spatial interaction curves (SICs) at multiple scales
- Compare with classical spatial statistics methods (G-cross, mFPCA)

## Workflow

### Main Analysis Pipeline

Run these scripts in order:

```
00_preprocess_data.R    → Download and preprocess single-cell CRC data
01_generate_data.R      → Create point patterns and quadrature approximations
02_fit_models.R         → Fit SHADE models (saves to data/)
03_analyze_results.R    → Extract SICs, compute uncertainty, generate main figures
04_g_cross_comparison.R → Compare SHADE with G-cross (generates gx_comparison.pdf)
05_mad_comparison.R     → Analyze spatial heterogeneity across patients/images
```

### mFPCA Comparison (Supplement)

For the SHADE vs mFPCA comparison in the supplement, run:

```r
source("crc_analysis/04b_method_comparison_new.R")
```

This orchestrator script runs modular analyses:
- `04b_00_setup.R` - Load data and prepare spatial patterns
- `04b_01_shade_analysis.R` - SHADE group-level detections
- `04b_02_mfpca.R` - mFPCA for G-cross and L-cross by group
- `04b_03_comparison.R` - Combine results into summary tables
- `04b_04_visualizations.R` - Create SHADE SIC figures

Output: `gcross_mfpca_grid_supplement.pdf`, `lcross_mfpca_grid_supplement.pdf`

## Cell Types

**Targets** (outcome cells):
- CTLs (CD8+ T cells)
- Memory CD4+ T cells
- Granulocytes

**Sources** (predictor cells):
- TAMs (tumor-associated macrophages)
- CAFs (cancer-associated fibroblasts)
- Vasculature
- Hybrid E/M cells
- Tumor cells

## Usage

```r
# Run main analysis
source("crc_analysis/00_preprocess_data.R")
source("crc_analysis/01_generate_data.R")
source("crc_analysis/02_fit_models.R")
source("crc_analysis/03_analyze_results.R")
source("crc_analysis/04_g_cross_comparison.R")
source("crc_analysis/05_mad_comparison.R")

# Run mFPCA comparison for supplement
source("crc_analysis/04b_method_comparison_new.R")
```

## Output

- `data/` - Processed point patterns, model fits, analysis summaries
- `CRC_analysis_paper/` - Generated figures:
  - `sic_CRC_all_ci.pdf` - All SICs with credible bands
  - `sic_CRC_diff.pdf` - CLR vs DII difference curves
  - `gx_comparison.pdf` - G-cross comparison heatmap
  - `sic_mad.pdf` - Heterogeneity visualization
  - `*_mfpca_grid_supplement.pdf` - mFPCA supplement figures
