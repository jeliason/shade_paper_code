# SHADE vs G-cross Comparison Simulation

This directory contains code for comparing the performance of SHADE models against traditional G-cross function analysis for detecting spatial interactions between cell types. The simulation evaluates detection power across different density conditions and sample sizes.

## Purpose

The G-cross function is a classical spatial statistics method for analyzing cross-type spatial interactions, while SHADE provides a hierarchical modeling approach. This simulation:

1.  Compares detection power between SHADE, flat SHADE, and G-cross methods
2.  Evaluates performance across varying T cell and tumor cell density conditions
3.  Assesses the effect of sample size (number of images per patient) on detection accuracy

## Workflow

1.  **simulation_runner.R**\
    Main simulation script that generates hierarchical spatial patterns with known group-level differences, fits both SHADE and flat models, runs G-cross analysis, and calculates detection power metrics.

2.  **gx_analysis.R**\
    Analyzes simulation results across all conditions, creates comparative visualizations showing proportion of correct detections for each method.

3.  **slurm.sh**\
    SLURM job array script for HPC execution across 360 simulation conditions.

## Usage

For local execution:

``` r
# Set environment to local
Sys.setenv(SYSTEM_ENV="laptop")

# Run single simulation
source("simulation_runner.R")

# Analyze results
source("gx_analysis.R")
```

For HPC execution:

``` bash
# Submit job array to SLURM
sbatch slurm.sh
```

## Simulation Design

The simulation varies: - **T cell density**: high (150 cells) vs low (15 cells) - **Tumor cell density**: high (150 cells) vs low (15 cells)\
- **Number of images per patient**: 1, 2, or 3 - **Replications**: 30 per condition

Each condition generates hierarchical spatial patterns where: - Responders show clustering (negative spatial interaction coefficients) - Non-responders show repulsion (positive spatial interaction coefficients)

## Output

-   **./data/**: Contains simulation results with detection power metrics for each method
-   **./gx_figures/**: Contains comparative visualizations including example patterns, G-cross functions, and spatial interaction curves
