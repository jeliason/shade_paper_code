# Flat Model Comparison Simulation

This directory contains code for simulating and analyzing the performance of flat (non-hierarchical) models compared to hierarchical models. The simulation evaluates the tradeoffs between model complexity and performance.

## Workflow

1. `01_generate_data.R`: Generate synthetic data for model comparison
2. `02_fit_models.R`: Fit both hierarchical and flat models to the data
3. `03_analyze_results.R`: Compare model performance
4. `04_create_figures.R`: Create figures for publication

SLURM scripts are provided for running on a high-performance computing cluster.