# SHADE Paper Code

This repository contains code for reproducing the analyses and figures from the [SHADE (Spatial Hierarchical Asymmetry via Directional Estimation) paper](https://doi.org/10.1371/journal.pcbi.1013930), published in *PLOS Computational Biology*. The codebase includes simulation studies and real data analyses that demonstrate the capabilities of the SHADE method.

## Overview

SHADE is a statistical method for modeling asymmetric spatial associations between cell types in tissue images using a hierarchical Bayesian framework. It captures directional spatial interactions using nonparametric spatial interaction curves (SICs) and models variation across biological levels (images, patients, cohorts).

## Directory Structure

### Simulation Studies

Each simulation follows a standardized 4-step workflow (`01_generate_data.R` → `02_fit_models.R` → `03_analyze_results.R` → `04_create_figures.R`):

-   **sim_shade_comparison/**: Comparisons of SHADE detection power against G-cross and K-cross envelope tests
-   **sim_flat_model/**: Comparisons between hierarchical and flat (non-hierarchical) models
-   **sim_timing/**: Computational scaling experiments across cell counts
-   **sim_sample_size/**: Simulations evaluating sample size effects on model performance
-   **sim_dummy_points/**: Parameter recovery simulations with varying quadrature densities

### Real Data Analysis

-   **crc_analysis/**: Analysis of colorectal cancer tissue data from Schürch et al. (2020)
-   **compartment_analysis/**: Sensitivity analysis with compartment covariates (tumor vs stroma)

### Supporting Files

-   **scripts/**: HPC workflow automation tools (see `scripts/README.md` for detailed usage)
-   **summary_figures/**: Scripts for generating overview figures (Figures 1, 2, S1, S2)
-   **utils.R**: Shared utility functions (density smoothing, path handling)

See the `README.md` in each subdirectory for specific details.

## Running the Code

### On HPC (Recommended for Simulations)

The `scripts/` folder contains Python tools for automating HPC workflows. See `scripts/README.md` for setup and usage:

```bash
# Setup (one-time)
./scripts/setup_env.sh
# Edit .env with your HPC details

# Typical workflow
source .venv/bin/activate
python scripts/hpc.py install              # Install R packages on HPC
python scripts/hpc.py submit sim_timing 01_generate_data
python scripts/hpc.py status
python scripts/hpc.py fetch sim_timing
```

### Running Locally

To run analyses on your local machine, set the environment variable before executing scripts:

```r
Sys.setenv(SYSTEM_ENV = "laptop")
```

This adjusts data paths and may reduce the scope of simulations for tractability.

## Requirements

-   R (>= 4.1.0)
-   [SHADE package](https://github.com/jeliason/SHADE)
-   Additional R packages:
    -   **tidyverse**: Data manipulation and visualization
    -   **spatstat**: Spatial point pattern analysis
    -   **Matrix**: Sparse matrix operations
    -   **posterior**: Working with posterior samples
    -   **cmdstanr**: Interface to Stan
    -   **ggdist**: Visualizing distributions
    -   **patchwork**: Combining plots

For HPC workflow tools:
-   Python 3.8+
-   See `scripts/requirements.txt`

## License

This code is released under the MIT license.
