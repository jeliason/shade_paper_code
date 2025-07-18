# Flat Model Comparison Simulation

This directory contains code for simulating and analyzing the performance of flat (non-hierarchical) models compared to hierarchical models using the SHADE framework.

## Purpose

Hierarchical models account for multiple levels of variation (e.g., within and between images/patients/cohorts), while flat models treat all observations as independent. This simulation:

1.  Compares parameter recovery between hierarchical and flat SHADE models
2.  Evaluates how well each model type captures the true spatial interaction patterns

## Workflow

1.  **01_generate_data.R**\
    Simulates hierarchical point patterns with known group- and individual-level parameters. Saves data for both hierarchical and flat model comparisons.

2.  **02_fit_models.R**\
    Fits SHADE models with hierarchical and flat structures using MCMC and variational inference. Saves results as `.rds` files.

3.  **03_analyze_results.R**\
    Compares recovery accuracy, uncertainty, and computational efficiency between model types. Summarizes performance metrics.

4.  **04_create_figures.R**\
    Produces figures highlighting differences in recovery and uncertainty between hierarchical and flat models.

## Usage

First, compile the flat model:

``` r
cmdstanr::cmdstan_model("sim_flat_model/SHADE_flat.stan",cpp_options = list(stan_threads = TRUE))
```

For local execution:

``` r
# Set environment to local
Sys.setenv(SYSTEM_ENV="laptop")

# Run scripts in sequence
source("01_generate_data.R")
source("02_fit_models.R")
source("03_analyze_results.R")
source("04_create_figures.R")
```

For HPC execution:

``` bash
SYSTEM_ENV="HPC"
# Submit jobs to SLURM
sbatch 01_generate_data.slurm
sbatch 02_fit_models.slurm
sbatch 03_analyze_results.slurm
```

## Output

-   **./data/**: Contains generated hierarchical point patterns and true parameters
-   **./sim_flat_model_figures/**: Contains comparative visualizations between hierarchical and flat models
