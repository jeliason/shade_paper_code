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
    Pre-compiles the `simple_shade.stan` model using the HPC toolchain. Must be run before 01_simulation_runner.R on HPC.

3.  **01_simulation_runner.R**\
    Main simulation script that, for each grid point (density × sample size), generates two sets of patterns: (1) power patterns with true spatial clustering (cohort_mean = [1.5, 1.0, 0.5]), and (2) null patterns with no spatial association (cohort_mean = [0, 0, 0]). Runs all 4 methods (SHADE hierarchical, SHADE flat, G-cross, K-cross) on both sets. Calculates power, coverage, and Type I error for each condition. Saves combined results to `sim_{idx}.rds`.

4.  **02_analyze_results.R**\
    Loads simulation results from all grid points, computes power and coverage summaries, computes Type I error rates by condition, and saves `analysis_summary.rds`.

5.  **03_create_figures.R**\
    Loads `analysis_summary.rds` and creates comparative visualizations showing detection power, coverage, and Type I error rates by data regime for each method.

6.  **02_robustness_hardcore_power.R**\
    Hard-core robustness check: generates patterns using SHADE model with true spatial clustering, applies hard-core thinning (r_hc ≈ 10 μm), then tests whether all methods recover spatial interactions at distances > r_hc despite Poisson assumption violations.

7.  **03_robustness_hardcore_null.R**\
    Hard-core null calibration: generates null patterns (no spatial association) with hard-core thinning, then verifies all methods maintain proper Type I error control (~5%) and do not mistake geometric crowding for spatial interaction.

## Usage

For local execution:

``` r
# Set environment to local
Sys.setenv(SYSTEM_ENV="laptop")

# Run main simulation (power + null calibration)
source("sim_shade_comparison/01_simulation_runner.R")

# Analyze results and create summary
source("sim_shade_comparison/02_analyze_results.R")

# Create figures
source("sim_shade_comparison/03_create_figures.R")
```

For HPC execution:

``` bash
# Step 1: Compile Stan model (run once)
python scripts/hpc.py submit sim_shade_comparison 00_setup

# Step 2: Wait for setup to complete, then submit simulations
python scripts/hpc.py submit sim_shade_comparison 01_simulation_runner

# Step 3: Check status
python scripts/hpc.py status sim_shade_comparison

# Step 4: Fetch results
python scripts/hpc.py fetch sim_shade_comparison

# Step 5: Run analysis locally (creates analysis_summary.rds and figures)
source("sim_shade_comparison/02_analyze_results.R")
source("sim_shade_comparison/03_create_figures.R")
```

## Simulation Design

### Original Grid (Full Power/Coverage Simulation)

The simulation varies:
- **T cell density**: high (150 cells) vs low (15 cells)
- **Tumor cell density**: high (150 cells) vs low (15 cells)
- **Number of images per patient**: 1, 2, or 3
- **Replications**: 30 per condition
- **Total jobs**: 2 × 2 × 3 × 30 = 360

All patients have clustering pattern (single group). Uses centered parameterization with variational inference.

### Experimental Grid (Model Configuration Testing)

Tests different SHADE model fitting approaches:
- **Method**: variational (ADVI) vs sampling (MCMC)
- **Parameterization**: centered vs non-centered
- **Prior scale**: σ_scale = 1 vs 5
- **T cell density**: high vs low (to test data quantity effect)
- **Tumor cell density**: high only
- **Number of images per patient**: 3 only (best case)
- **Replications**: 3 per condition
- **Total jobs**: 2 × 2 × 2 × 2 × 1 × 1 × 3 = 48

Purpose: Understand how inference method, parameterization, and prior scales affect:
1. Global coefficient estimates (bias/shrinkage)
2. Coverage rates (should be ~95%)
3. Power (detection rates)
4. Computational time

Note: To switch between grids, comment/uncomment the appropriate `grid <- expand.grid(...)` section in `01_simulation_runner.R`.

## Output

-   **./data/**: Contains simulation results with detection power metrics for each method
-   **./comparison_figures/**: Contains comparative visualizations including example patterns, G-cross/K-cross functions, and spatial interaction curves
