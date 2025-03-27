#!/bin/bash
#SBATCH --job-name=fit_simulated_no_hier        # Job name
#SBATCH --output=logs/fit_no_hier_task_%A_%a.out   # Standard output log
#SBATCH --error=logs/fit_no_hier_task_%A_%a.err    # Standard error log
#SBATCH --array=1-100                   # Array range (adjust as needed)
#SBATCH --time=5:00                # Time limit hh:mm:ss
#SBATCH --mem=4G                       # Memory per task
#SBATCH --cpus-per-task=2              # CPUs per task
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joelne@umich.edu

# Load the required R module (if needed)
module load Rgeospatial

# Run the R script with the SLURM_ARRAY_TASK_ID as an argument
Rscript scripts/simulations_no_hier/fit_simulated_no_hier.R $SLURM_ARRAY_TASK_ID
