# Sample Size Effect Simulation

This directory contains code for simulating and analyzing the effect of sample size (number of point patterns) on SHADE model performance. The simulation explores how model accuracy and uncertainty change with varying sample sizes.

## Purpose

The number of point patterns (e.g., tissue images) available for analysis can significantly impact model performance and inference. This simulation:

1. Evaluates how parameter recovery changes with increasing sample sizes
2. Quantifies uncertainty reduction as more data becomes available

## Workflow

1. **01_generate_data.R**  
   Simulates point patterns with fixed interaction parameters across varying sample sizes. Saves datasets with ground truth for comparison.

2. **02_fit_models.R**  
   Fits SHADE models using MCMC and variational inference across sample sizes, with consistent settings. Saves model outputs as `.rds` files.

3. **03_analyze_results.R**  
   Assesses parameter recovery, uncertainty, and computation time across sample sizes. Summarizes results for plotting.

4. **04_create_figures.R**  
   Generates figures showing effects of sample size on accuracy, uncertainty, and convergence.

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
SYSTEM_ENV="HPC"
# Submit jobs to SLURM
sbatch 01_generate_data.slurm
sbatch 02_fit_models.slurm
sbatch 03_analyze_results.slurm
```

## Output

- **./data/**: Contains generated point patterns of varying sample sizes
- **./sim_sample_size_figures/**: Contains visualizations showing the effect of sample size on model performance