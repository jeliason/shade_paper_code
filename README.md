# SHADE Paper Code

This repository contains code for reproducing the analyses and figures from the SHADE (Spatial Hierarchical Asymmetry via Directional Estimation) paper. The codebase includes simulation studies and real data analyses that demonstrate the capabilities of the SHADE method.

## Overview

SHADE is a statistical method for modeling asymmetric spatial associations between cell types in tissue images using a hierarchical Bayesian framework. It captures directional spatial interactions using nonparametric spatial interaction curves (SICs) and models variation across biological levels (images, patients, cohorts).

Key features of SHADE:
- Hierarchical modeling of spatial interactions at multiple biological scales
- Nonparametric spatial interaction curves using radial basis functions
- Support for asymmetric/directional spatial relationships
- Bayesian inference via Stan (MCMC and variational methods)
- Robust parameter estimation with uncertainty quantification

## Repository Structure

- **sim_dummy_points/**: Parameter recovery simulations using synthetic point patterns with dummy points
- **sim_sample_size/**: Simulations evaluating the effect of sample size on model performance
- **sim_flat_model/**: Comparisons between hierarchical and flat (non-hierarchical) models
- **crc_analysis/**: Analysis of colorectal cancer tissue data
- **figures/**: Scripts for generating summary figures across simulations
- **utils.R**: Shared utility functions for spatial analysis

## Simulation Workflow

Each simulation directory follows a consistent workflow:

1. **01_generate_data.R**: Generate synthetic data with known parameters
2. **02_fit_models.R**: Fit SHADE models to the generated data
3. **03_analyze_results.R**: Analyze model performance and parameter recovery
4. **04_create_figures.R**: Create publication-quality figures

SLURM scripts (`*.slurm`) are provided for running computationally intensive parts on a high-performance computing cluster.

## Utility Files

- **utils.R**: Shared utility functions used across analyses (density smoothing, radial basis functions)
- **summary_figure_subplots.R**: Code for generating summary figures
- **demo_paper_figures.R**: Script to reproduce key paper figures

## Requirements

- R (>= 4.1.0)
- SHADE package (available at https://github.com/jeliason/SHADE)
- Additional R packages:
  - **tidyverse**: Data manipulation and visualization
  - **spatstat**: Spatial point pattern analysis
  - **Matrix**: Sparse matrix operations
  - **posterior**: Working with posterior samples
  - **cmdstanr**: Interface to Stan
  - **ggdist**: Visualizing distributions

## Usage

For local execution:
```r
# Set environment to local
Sys.setenv(SYSTEM_ENV="laptop")

# Run scripts in sequence
source("01_generate_data.R")
source("02_fit_models.R")
source("03_analyze_results.R")
source("04_create_figures.R")
```

For HPC execution, submit the corresponding SLURM scripts.

See the README.md in each subdirectory for specific details about each simulation or analysis.

## License

This code is released under the MIT license.