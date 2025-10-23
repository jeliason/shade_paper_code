# SHADE vs Classical Methods Comparison Simulation

This directory contains code for comparing the performance of SHADE models against classical spatial statistics methods (G-cross and Ripley's K-function) for detecting spatial interactions between cell types. The simulation evaluates detection power across different density conditions and sample sizes.

## Purpose

The G-cross function and Ripley's K-function are classical spatial statistics methods for analyzing cross-type spatial interactions, while SHADE provides a hierarchical modeling approach. This simulation:

1.  Compares detection power between SHADE, flat SHADE, G-cross, and K-cross methods
2.  Evaluates performance across varying T cell and tumor cell density conditions
3.  Assesses the effect of sample size (number of images per patient) on detection accuracy

## Workflow

1.  **comparison_functions.R**\
    Shared functions used by both power and null calibration simulations. Contains data generation, SHADE fitting, G-cross, and K-cross analysis functions. Source this file to reuse functions.

2.  **00_setup.R**\
    Pre-compiles the `simple_shade.stan` model using the HPC toolchain. Must be run before simulation_runner.R on HPC.

3.  **simulation_runner.R**\
    Main power simulation script that generates hierarchical spatial patterns with known group-level differences, fits both SHADE and flat models, runs G-cross and K-cross analyses, and calculates detection power metrics.

4.  **null_calibration.R**\
    Null scenario simulation to test type I error calibration. Generates data under the null hypothesis (no spatial association) and measures false positive rates for each method. Expected rate: ~5% at α=0.05.

5.  **comparison_analysis.R**\
    Analyzes simulation results across all conditions, creates comparative visualizations showing proportion of correct detections for each method.

## Usage

For local execution:

``` r
# Set environment to local
Sys.setenv(SYSTEM_ENV="laptop")

# Run power simulation
source("sim_shade_comparison/simulation_runner.R")

# Run null calibration check
source("sim_shade_comparison/null_calibration.R")

# Analyze results
source("sim_shade_comparison/comparison_analysis.R")
```

For HPC execution:

``` bash
# Step 1: Compile Stan model (run once)
python scripts/hpc.py submit sim_shade_comparison 00_setup

# Step 2: Wait for setup to complete, then submit simulations
python scripts/hpc.py submit sim_shade_comparison simulation_runner

# Step 3: Check status
python scripts/hpc.py status sim_shade_comparison

# Step 4: Fetch results
python scripts/hpc.py fetch sim_shade_comparison

# Step 5: Run analysis
python scripts/hpc.py submit sim_shade_comparison comparison_analysis
```

## Simulation Design

The simulation varies: - **T cell density**: high (150 cells) vs low (15 cells) - **Tumor cell density**: high (150 cells) vs low (15 cells)\
- **Number of images per patient**: 1, 2, or 3 - **Replications**: 30 per condition

Each condition generates hierarchical spatial patterns where: - Responders show clustering (negative spatial interaction coefficients) - Non-responders show repulsion (positive spatial interaction coefficients)

## Output

-   **./data/**: Contains simulation results with detection power metrics for each method
-   **./comparison_figures/**: Contains comparative visualizations including example patterns, G-cross/K-cross functions, and spatial interaction curves
