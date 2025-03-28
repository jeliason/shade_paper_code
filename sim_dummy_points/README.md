# Parameter Recovery Simulation with Dummy Points

This directory contains code for simulating and analyzing parameter recovery using dummy point patterns. The simulation assesses the ability of the SHADE model to recover known parameters from synthetic data.

## Workflow

1. `01_generate_data.R`: Generate synthetic point patterns with known parameters
2. `02_fit_models.R`: Fit SHADE models to the generated data
3. `03_analyze_results.R`: Analyze parameter recovery performance
4. `04_create_figures.R`: Create figures for publication

SLURM scripts are provided for running on a high-performance computing cluster.