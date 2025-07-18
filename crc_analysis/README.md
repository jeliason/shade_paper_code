# Colorectal Cancer (CRC) Analysis

This directory contains code for analyzing real CRC tissue data using the SHADE model.

## Purpose

The goal is to uncover how immune and tumor cells are spatially organized in the tumor microenvironment and what this reveals about cancer biology and patient outcomes. Specifically, the analysis:

-   Models directional spatial interactions between key cell types\
-   Identifies hierarchical spatial patterns across images, patients, and cohorts\
-   Characterizes spatial interaction curves (SICs) at multiple scales\
-   Links spatial patterns to biological mechanisms in CRC

## Workflow

The analysis is structured into three scripts:

1.  **01_generate_data.R**\
    Prepares point pattern data, creates quadrature approximations, organizes it hierarchically, and performs quality checks.

2.  **02_fit_models.R**\
    Fits SHADE models using MCMC or variational inference, accounting for patient- and image-level variation. Outputs model fits as `.rds` files.

3.  **03_analyze_results.R**\
    Extracts and interprets spatial interaction patterns, quantifies uncertainty, and generates summaries and visualizations.
    
3.  **04_gcross_comparison.R**\
    Compares results from SHADE with those from a G-cross point pattern analysis.

## Usage

To run locally:

``` r
source("01_generate_data.R")
source("02_fit_models.R")
source("03_analyze_results.R")
source("04_gcross_comparison.R")
```

## Output

-   **./data/**: Contains processed CRC point pattern data
-   **./CRC_analysis_paper/**: Contains visualizations of spatial interaction patterns in CRC tissue
