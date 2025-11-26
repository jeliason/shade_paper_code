# Computational Scaling Study

This directory contains code for measuring SHADE's computational runtime as a function of cell count. The simulation evaluates both feature construction time and full model fitting time across datasets ranging from 1,000 to 100,000 cells.

## Purpose

As multiplexed imaging technologies scale to hundreds of thousands of cells per sample, it is important to understand how SHADE's computational requirements grow with dataset size. This simulation:

1. Measures feature construction time (distance matrix computation and feature engineering) across cell counts
2. Measures full model fitting time (end-to-end) using variational inference
3. Characterizes empirical scaling behavior and compares to theoretical expectations

## Workflow

1. **01_generate_data.R**
   Simulates spatial point patterns with 3 cell types across varying cell counts (1K, 2.5K, 5K, 10K, 25K, 50K, 100K). Generates 20 replicates per condition with fixed window size so density increases with cell count. Saves datasets for model fitting.

2. **02_fit_models.R**
   Fits SHADE models using variational inference (1000 draws) and records timing for both feature construction and full fitting process. Saves model outputs and timing data as `.rds` files.

3. **03_analyze_results.R**
   Aggregates timing results across replicates, computes mean ± SD for each cell count, and calculates empirical scaling exponents via log-log regression. Saves standardized `analysis_summary.rds` containing all processed data needed for figures.

4. **04_create_figures.R**
   Loads `analysis_summary.rds` and generates scaling plots showing feature construction time and total fitting time vs cell count. Saves figures to `manuscript/images/sim_timing/`.

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
python scripts/hpc.py submit sim_timing 01_generate_data
python scripts/hpc.py submit sim_timing 02_fit_models
python scripts/hpc.py submit sim_timing 03_analyze_results

# Fetch results (only downloads analysis_summary.rds, not raw data)
python scripts/hpc.py fetch sim_timing

# Generate figures locally
Rscript sim_timing/04_create_figures.R
```

## Simulation Design

- **Cell counts**: 1000, 2500, 5000, 10000, 25000, 50000, 100000 (total across all images)
- **Hierarchical structure**: 10 patients, 4 images per patient (40 total images)
- **Cell types**: 3 (2 source types, 1 target type)
- **Basis functions**: 3 radial basis functions
- **Window size**: 1500 × 1500 (fixed)
- **Replicates**: 20 per cell count condition
- **Inference method**: Variational inference with 1000 draws

## Output

- **./data/** (HPC only): Raw simulation data (point patterns, Stan data, timing results)
- **./data/analysis_summary.rds**: Standardized analysis output containing:
  - `timing_data`: Individual timing measurements for each replicate
  - `timing_summary`: Aggregated statistics by cell count
  - `scaling_exponents`: Empirical scaling exponents from log-log regression
  - `fit_feature`, `fit_total`: Linear model objects
  - `metadata`: Analysis metadata (date, simulation count, etc.)
- **manuscript/images/sim_timing/**: Scaling visualizations for manuscript supplement
