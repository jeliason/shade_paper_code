#!/bin/bash
#SBATCH --job-name=gcross_comparison         # Job name
#SBATCH --account=ukarvind99
#SBATCH --output=logs/gx_%A_%a.out   # Standard output log
#SBATCH --error=logs/gx_%A_%a.err    # Standard error log
#SBATCH --array=1-200                   # Array range (adjust as needed)
#SBATCH --time=15:00                # Time limit hh:mm:ss
#SBATCH --mem=4G                       # Memory per task
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joelne@umich.edu

# Load the required R module (if needed)
module load Rgeospatial

# Run the R script with the SLURM_ARRAY_TASK_ID as an argument
Rscript sim_shade_gcross/simulation_runner.R $SLURM_ARRAY_TASK_ID