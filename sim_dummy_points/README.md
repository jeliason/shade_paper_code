# Parameter Recovery Simulation with Dummy Points

This directory contains code for simulating and analyzing parameter recovery using synthetic point patterns with dummy points. The simulation assesses the ability of the SHADE model to recover known parameters from synthetic data.

## Purpose

Dummy points in this context are used to approximate quadrature for estimating continuous spatial interaction functions. This simulation suite investigates:

1. The effect of the ratio of dummy points to real points (quadrature density)
2. How varying the number of points per type impacts model performance
3. Parameter recovery accuracy under different quadrature approximation settings

## Workflow

1. **01_generate_data.R**: 
   - Generates synthetic point patterns with specific spatial interaction parameters
   - Creates dummy points with varying densities relative to real points
   - Saves the data with ground truth parameters for later comparison

2. **02_fit_models.R**: 
   - Fits SHADE models to the generated data using radial basis functions
   - Uses `run_SHADE_model()` with MCMC sampling for full posterior inference
   - Evaluates models with different dummy point ratios
   - Saves model fits as RDS files

3. **03_analyze_results.R**: 
   - Loads fitted models and ground truth
   - Computes parameter recovery metrics (RMSE, coverage)
   - Compares performance across different dummy point densities
   - Creates summary statistics for figures

4. **04_create_figures.R**: 
   - Generates publication-quality figures showing parameter recovery
   - Creates visualizations of spatial interaction curves
   - Compares true vs. estimated parameters across different quadrature densities

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

For HPC execution:
```bash
# Submit jobs to SLURM
sbatch 01_generate_data.slurm
sbatch 02_fit_models.slurm
sbatch 03_analyze_results.slurm
```

## Output

- **./data/**: Contains generated point patterns, dummy point ratios, and true parameters
- **./figures/**: Contains visualizations of results showing parameter recovery across different quadrature approximation settings