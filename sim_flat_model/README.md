# Flat Model Comparison Simulation

This directory contains code for simulating and analyzing the performance of flat (non-hierarchical) models compared to hierarchical models using the SHADE framework. The simulation evaluates the tradeoffs between model complexity and performance.

## Purpose

Hierarchical models account for multiple levels of variation (e.g., within and between images/patients/cohorts), while flat models treat all observations as independent. This simulation:

1. Compares parameter recovery between hierarchical and flat SHADE models
2. Evaluates how well each model type captures the true spatial interaction patterns
3. Assesses computational efficiency and statistical performance tradeoffs
4. Demonstrates when hierarchical modeling provides meaningful advantages

## Workflow

1. **01_generate_data.R**: 
   - Generates hierarchical synthetic point patterns with group-level and individual-level parameters
   - Creates data with known hierarchical structure to test model performance
   - Saves the data with ground truth parameters for both hierarchical and flat interpretations

2. **02_fit_models.R**: 
   - Fits both hierarchical and flat SHADE models to the same datasets
   - Uses `run_SHADE_model()` with appropriate hierarchical/flat specifications
   - Implements both MCMC sampling and variational inference methods
   - Saves model fits as RDS files

3. **03_analyze_results.R**: 
   - Compares parameter recovery between hierarchical and flat models
   - Evaluates computational efficiency (runtime, memory usage)
   - Assesses statistical performance (RMSE, coverage, posterior uncertainty)
   - Creates summary statistics for figures

4. **04_create_figures.R**: 
   - Generates comparison figures between hierarchical and flat models
   - Visualizes differences in parameter recovery and uncertainty
   - Creates publication-quality figures showing when hierarchical modeling matters

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

- **./data/**: Contains generated hierarchical point patterns and true parameters
- **./figures/**: Contains comparative visualizations between hierarchical and flat models