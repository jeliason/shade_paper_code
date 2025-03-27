#!/bin/bash
#SBATCH --job-name=fit_simulated_dummy         # Job name
#SBATCH --output=logs/fit_dummy_task_%A_%a.out   # Standard output log
#SBATCH --error=logs/fit_dummy_task_%A_%a.err    # Standard error log
#SBATCH --array=1-125                   # Array range (adjust as needed)
#SBATCH --time=5:00                # Time limit hh:mm:ss
#SBATCH --mem=4G                       # Memory per task
#SBATCH --cpus-per-task=2              # CPUs per task
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joelne@umich.edu

# Load the required R module (if needed)
module load Rgeospatial

# Run the R script with the SLURM_ARRAY_TASK_ID as an argument
Rscript scripts/simulations_dummy_hier_asymm/fit_simulated_dummy.R $SLURM_ARRAY_TASK_ID
