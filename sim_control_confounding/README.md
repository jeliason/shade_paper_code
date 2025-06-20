# Sample Size Effect Simulation

This directory contains code for simulating and analyzing the effect of sample size (number of point patterns) on SHADE model performance. The simulation explores how model accuracy and uncertainty change with varying sample sizes.

## Purpose

The number of point patterns (e.g., tissue images) available for analysis can significantly impact model performance and inference. This simulation:

1. Evaluates how parameter recovery changes with increasing sample sizes
2. Quantifies uncertainty reduction as more data becomes available
3. Determines minimum sample sizes needed for reliable inference
4. Assesses computational scaling with increased sample sizes

## Workflow

1. **01_generate_data.R**: 
   - Generates synthetic point patterns with known spatial interaction parameters
   - Creates datasets with varying numbers of point patterns (small to large sample sizes)
   - Uses consistent generation parameters to isolate sample size effects
   - Saves the data with ground truth parameters for later comparison

2. **02_fit_models.R**: 
   - Fits SHADE models to datasets of different sample sizes
   - Uses `run_SHADE_model()` with consistent settings across sample sizes
   - Implements both MCMC sampling and variational inference methods
   - Applies appropriate adaptation for different sample sizes
   - Saves model fits as RDS files

3. **03_analyze_results.R**: 
   - Quantifies parameter recovery accuracy across sample sizes
   - Measures posterior uncertainty (credible interval width) as a function of sample size
   - Evaluates computational requirements for different sample sizes
   - Creates summary statistics for figures

4. **04_create_figures.R**: 
   - Creates figures showing parameter recovery vs. sample size
   - Generates uncertainty quantification visualizations
   - Plots convergence metrics across sample sizes
   - Creates publication-quality figures demonstrating sample size effects

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

- **./data/**: Contains generated point patterns of varying sample sizes
- **./figures/**: Contains visualizations showing the effect of sample size on model performance