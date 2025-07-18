# Parameter Recovery Simulation with Dummy Points

This directory contains code for simulating and analyzing parameter recovery using synthetic point patterns with dummy points. The simulation assesses the ability of the SHADE model to recover known parameters from synthetic data.

## Purpose

Dummy points in this context are used to approximate quadrature for estimating continuous spatial interaction functions. This simulation suite investigates:

1. The effect of the ratio of dummy points to real points (quadrature density)
2. How varying the number of points per type impacts model performance
3. Parameter recovery accuracy under different quadrature approximation settings

## Workflow

1. **01_generate_data.R**  
   Simulates spatial point patterns with known parameters, adds dummy points with varying densities, and saves ground truth for comparison.

2. **02_fit_models.R**  
   Fits SHADE models using MCMC or variational inference, tests different dummy point ratios, and saves results as `.rds` files.

3. **03_analyze_results.R**  
   Loads model fits and ground truth, computes recovery metrics (e.g. RMSE, coverage), and summarizes performance across conditions.

4. **04_create_figures.R**  
   Generates figures of parameter recovery and interaction curves, comparing estimated vs. true parameters.

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
SYSTEM_ENV="HPC"

sbatch 01_generate_data.slurm
sbatch 02_fit_models.slurm
sbatch 03_analyze_results.slurm
```

## Output

- **./data/**: Contains generated point patterns, dummy point ratios, and true parameters
- **./sim_dummy_points_figures/**: Contains visualizations of results showing parameter recovery across different quadrature approximation settings