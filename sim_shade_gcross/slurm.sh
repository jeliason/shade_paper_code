#!/bin/bash
#SBATCH --job-name=gcross_comparison         # Job name
#SBATCH --output=logs/gx_%A_%a.out   # Standard output log
#SBATCH --error=logs/gx_%A_%a.err    # Standard error log
#SBATCH --array=1-360                   # Array range (adjust as needed)
#SBATCH --time=45:00                # Time limit hh:mm:ss
#SBATCH --mem=4G                       # Memory per task
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

# Load the required R module (if needed)
module load R

# Run the R script with the SLURM_ARRAY_TASK_ID as an argument
Rscript sim_shade_gcross/simulation_runner.R $SLURM_ARRAY_TASK_ID