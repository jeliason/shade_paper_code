# Colorectal Cancer (CRC) Analysis

This directory contains code for analyzing real CRC tissue data using the SHADE model.

## Purpose

The goal is to uncover how immune and tumor cells are spatially organized in the tumor microenvironment and what this reveals about cancer biology and patient outcomes. Specifically, the analysis:

-   Models directional spatial interactions between key cell types\
-   Identifies hierarchical spatial patterns across images, patients, and cohorts\
-   Characterizes spatial interaction curves (SICs) at multiple scales\
-   Links spatial patterns to biological mechanisms in CRC

## Workflow

The analysis is structured into sequential scripts:

0. **00_preprocess_data.R**\
    Downloads and preprocesses single-cell CRC data as well as clinical patient annotations.

1.  **01_generate_data.R**\
    Prepares point pattern data, creates quadrature approximations, organizes it hierarchically, and performs quality checks.

2.  **02_fit_models.R**\
    Fits SHADE models using MCMC or variational inference, accounting for patient- and image-level variation. Outputs model fits as `.rds` files.

3.  **03_analyze_results.R**\
    Extracts and interprets spatial interaction patterns, quantifies uncertainty, and generates summaries and visualizations.

4.  **04_gcross_comparison.R**\
    Compares results from SHADE with those from a G-cross point pattern analysis.

5.  **05_mad_comparison.R**\
    Analyzes heterogeneity in spatial patterns across patients and images.

6.  **06_compartment_sensitivity.R**\
    Tests robustness to spatial compartments by adding density-based compartments as covariates and refitting models.

7.  **07_compartment_comparison.R**\
    Compares original vs compartment-adjusted SIC estimates to assess sensitivity to unmeasured spatial structure.

## Usage

To run the main analysis locally:

``` r
source("crc_analysis/00_preprocess_data.R")
source("crc_analysis/01_generate_data.R")
source("crc_analysis/02_fit_models.R")
source("crc_analysis/03_analyze_results.R")
source("crc_analysis/04_gcross_comparison.R")
source("crc_analysis/05_mad_comparison.R")
```

To run the compartment sensitivity analysis (for revision):

``` r
# First ensure main analysis is complete, then:
source("crc_analysis/06_compartment_sensitivity.R")  # Refit with compartments
source("crc_analysis/07_compartment_comparison.R")   # Compare results
```

This will:
1. Add density-based compartments to each image (median split of local tumor density)
2. Refit all models including compartment as a covariate
3. Compare original vs compartment-adjusted SIC estimates
4. Generate comparison figures and tables for the supplement

## Output

-   **./data/**: Contains processed CRC point pattern data
-   **./CRC_analysis_paper/**: Contains visualizations of spatial interaction patterns in CRC tissue
