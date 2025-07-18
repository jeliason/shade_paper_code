# SHADE Paper Code

This repository contains code for reproducing the analyses and figures from the [SHADE (Spatial Hierarchical Asymmetry via Directional Estimation) manuscript](https://doi.org/10.1101/2025.06.24.661393). The codebase includes simulation studies and real data analyses that demonstrate the capabilities of the SHADE method.

## Overview

SHADE is a statistical method for modeling asymmetric spatial associations between cell types in tissue images using a hierarchical Bayesian framework. It captures directional spatial interactions using nonparametric spatial interaction curves (SICs) and models variation across biological levels (images, patients, cohorts).

-   **crc_analysis/**: Analysis of colorectal cancer tissue data
-   **sim_dummy_points/**: Parameter recovery simulations using synthetic point patterns with dummy points
-   **sim_sample_size/**: Simulations evaluating the effect of sample size on model performance
-   **sim_flat_model/**: Comparisons between hierarchical and flat (non-hierarchical) models
-   **sim_shade_gcross/**: Comparisons of SHADE model with existing Gcross methods to evaluate recovery of true interactions

See the README.md in each subdirectory for specific details about each simulation or analysis.

SLURM scripts (`*.slurm`) are provided for running computationally intensive parts on a high-performance computing cluster.

## Utility Files

-   **utils.R**: Shared utility functions used across analyses (density smoothing)
-   **summary_figures/summary_figure_subplots.R**: Code for generating Figure 1
-   **summary_figures/other_summary_figures.R**: Script to reproduce other summary paper figures - Figure 2 and Supplementary Figures S1 and S2

## Requirements

-   R (\>= 4.1.0)
-   SHADE package (available at <https://github.com/jeliason/SHADE>)
-   Additional R packages:
    -   **tidyverse**: Data manipulation and visualization
    -   **spatstat**: Spatial point pattern analysis
    -   **Matrix**: Sparse matrix operations
    -   **posterior**: Working with posterior samples
    -   **cmdstanr**: Interface to Stan
    -   **ggdist**: Visualizing distributions

## License

This code is released under the MIT license.
