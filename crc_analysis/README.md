# Colorectal Cancer (CRC) Analysis

This directory contains code for analyzing real colorectal cancer tissue data using the SHADE model. The analysis examines spatial interactions between different cell types in tumor microenvironments.

## Workflow

1. `01_generate_data.R`: Prepare CRC data for analysis
2. `02_fit_models.R`: Fit SHADE models to the CRC data
3. `03_analyze_results.R`: Analyze model results and biological implications

SLURM scripts are provided for running on a high-performance computing cluster.

The `get_fit_back.sh` script is used to retrieve model fit results from the compute cluster.