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
   Assesses parameter recovery, uncertainty, and computation time across sample sizes. Saves standardized `analysis_summary.rds` containing all processed data needed for figures.

4. **04_create_figures.R**
   Loads `analysis_summary.rds` and generates figures showing effects of sample size on accuracy, uncertainty, and convergence.

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
# Submit jobs
python scripts/hpc.py submit sim_sample_size 01_generate_data
python scripts/hpc.py submit sim_sample_size 02_fit_models
python scripts/hpc.py submit sim_sample_size 03_analyze_results

# Fetch results (only downloads analysis_summary.rds, not raw data)
python scripts/hpc.py fetch sim_sample_size

# Generate figures locally
Rscript sim_sample_size/04_create_figures.R
```

## Output

- **./data/** (HPC only): Raw simulation data (point patterns, model fits)
- **./data/analysis_summary.rds**: Standardized analysis output containing:
  - `rmse_tb`: RMSE metrics by coefficient and scale
  - `sic_tb`: Spatial interaction curve estimates for plotting
  - `metadata`: Analysis metadata (date, simulation count, etc.)
- **manuscript/images/sim_sample_size_figures/**: Visualizations showing the effect of sample size on model performance