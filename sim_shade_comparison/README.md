# SHADE vs Classical Methods Comparison Simulation

Compares SHADE detection power against classical spatial statistics methods (G-cross and Ripley's K-function) for detecting spatial interactions between cell types.

## Purpose

1. Compare detection power between SHADE, flat SHADE, G-cross, and K-cross methods
2. Evaluate performance across varying T cell and tumor cell density conditions
3. Assess the effect of sample size (number of images per patient) on detection accuracy
4. Test robustness to unmeasured spatial confounders (compartment structure)

## Workflow

### Main Comparison

```
00_setup.R               → Pre-compile Stan model (run once on HPC)
01_simulation_runner.R   → Main simulation: power patterns + null patterns
02_analyze_results.R     → Aggregate results, compute power/Type I error
03_create_figures.R      → Create comparison visualizations
04_create_example_figure.R → Create example 2x2 figure for manuscript
```

### Compartment Robustness Simulation

Tests SHADE's behavior when unmeasured spatial structure (compartments) affects cell densities:

```
01_simulation_runner_compartments.R → Simulate with compartment confounders
02_analyze_compartments.R           → Analyze compartment simulation results
03_create_compartment_figures.R     → Create compartment robustness figures
```

### Shared Files

- `comparison_functions.R` - Data generation, SHADE fitting, G-cross, K-cross functions
- `plot_confounder_diagnostics.R` - Diagnostic plotting utilities for compartment simulations
- `super_simple_shade.stan` - Simplified Stan model for simulations
- `slurm_config.yaml` - HPC job configuration

## Usage

For local execution:

```r
Sys.setenv(SYSTEM_ENV="laptop")

# Main comparison
source("sim_shade_comparison/01_simulation_runner.R")
source("sim_shade_comparison/02_analyze_results.R")
source("sim_shade_comparison/03_create_figures.R")
source("sim_shade_comparison/04_create_example_figure.R")

# Compartment robustness (optional)
source("sim_shade_comparison/01_simulation_runner_compartments.R")
source("sim_shade_comparison/02_analyze_compartments.R")
source("sim_shade_comparison/03_create_compartment_figures.R")
```

For HPC execution:

```bash
# Compile Stan model (once)
python scripts/hpc.py submit sim_shade_comparison 00_setup

# Run main simulation
python scripts/hpc.py submit sim_shade_comparison 01_simulation_runner
python scripts/hpc.py status sim_shade_comparison
python scripts/hpc.py fetch sim_shade_comparison

# Analyze locally
Rscript sim_shade_comparison/02_analyze_results.R
Rscript sim_shade_comparison/03_create_figures.R
```

## Simulation Design

### Main Comparison Grid

- **T cell density**: high (150 cells) vs low (15 cells)
- **Tumor cell density**: high (150 cells) vs low (15 cells)
- **Images per patient**: 1, 2, or 3
- **Replications**: 30 per condition
- **Total jobs**: 2 × 2 × 3 × 30 = 360

Each job generates:
1. Power patterns with true clustering (cohort_mean = [1.5, 1.0, 0.5])
2. Null patterns with no association (cohort_mean = [0, 0, 0])

### Compartment Robustness Grid

- **Compartment effect strength**: weak, moderate, strong
- **T cell density**: low, medium, high
- **Tumor cell density**: low, medium, high
- Tests whether SHADE incorrectly attributes compartment effects to cell-cell interactions

## Output

- `data/` - Simulation results (`.rds` files)
- `manuscript/images/comparison_figures/` - Generated figures:
  - `example_comparison_2x2.pdf` - Example pattern with all methods
  - `power_by_regime.pdf` - Detection power comparison
  - `type_i_error_by_regime.pdf` - Type I error rates
  - `coverage_by_regime.pdf` - Coverage comparison
  - `compartment_main_figure.pdf` - Compartment robustness results
