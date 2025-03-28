# Colorectal Cancer (CRC) Analysis

This directory contains code for analyzing real colorectal cancer tissue data using the SHADE model. The analysis examines spatial interactions between different cell types in tumor microenvironments to understand their ecological relationships.

## Purpose

Spatial organization of immune and tumor cells within the tumor microenvironment can provide insights into cancer biology and patient outcomes. This analysis:

1. Models directional/asymmetric spatial interactions between key cell types in CRC tissue
2. Identifies hierarchical patterns across patients and cohorts
3. Characterizes spatial interaction curves (SICs) at different biological scales
4. Relates spatial patterns to biological mechanisms in colorectal cancer

## Workflow

1. **01_generate_data.R**: 
   - Prepares CRC point pattern data for analysis
   - Creates appropriate quadrature approximation using dummy points
   - Organizes data into hierarchical structures (images nested within patients)
   - Performs data quality checks and preprocessing

2. **02_fit_models.R**: 
   - Fits SHADE models to the CRC data using radial basis functions
   - Implements hierarchical models to capture patient-level and image-level variation
   - Uses `run_SHADE_model()` with MCMC sampling for full posterior inference
   - Handles computational demands using high-performance computing
   - Saves model fits as RDS files

3. **03_analyze_results.R**: 
   - Extracts posterior distributions from fitted models
   - Characterizes cell-type interaction patterns and their biological implications
   - Quantifies uncertainty in spatial relationships
   - Creates biologically interpretable summaries of the results
   - Prepares data for visualization and further analysis

## Usage

For local execution (with subsampled data):
```r
# Set environment to local
Sys.setenv(SYSTEM_ENV="laptop")

# Run scripts in sequence
source("01_generate_data.R")
source("02_fit_models.R")
source("03_analyze_results.R")
```

For HPC execution:
```bash
# Submit jobs to SLURM
sbatch 01_generate_data.slurm
sbatch 02_fit_models.slurm
sbatch 03_analyze_results.slurm
```

To retrieve results from compute cluster:
```bash
./get_fit_back.sh
```

## Output

- **./data/**: Contains processed CRC point pattern data
- **./figures/**: Contains visualizations of spatial interaction patterns in CRC tissue

## References

- Original CRC dataset references
- SHADE methodology paper
- Related literature on spatial analysis in tumor microenvironments