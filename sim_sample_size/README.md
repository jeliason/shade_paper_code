# Sample Size Effect Simulation

This directory contains code for simulating and analyzing the effect of sample size (number of patterns) on model performance. The simulation explores how model accuracy changes with varying sample sizes.

## Workflow

1. `01_generate_data.R`: Generate synthetic data with different sample sizes
2. `02_fit_models.R`: Fit SHADE models to the generated data
3. `03_analyze_results.R`: Analyze model performance across sample sizes
4. `04_create_figures.R`: Create figures for publication

SLURM scripts are provided for running on a high-performance computing cluster.