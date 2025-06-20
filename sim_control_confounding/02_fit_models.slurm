#!/bin/bash
#SBATCH --job-name=fit_simulated_dsize         # Job name
#SBATCH --output=logs/fit_dsize_task_%A_%a.out   # Standard output log
#SBATCH --error=logs/fit_dsize_task_%A_%a.err    # Standard error log
#SBATCH --array=1-135                   # Array range (adjust as needed)
#SBATCH --time=5:00                # Time limit hh:mm:ss
#SBATCH --mem=4G                       # Memory per task
#SBATCH --cpus-per-task=2              # CPUs per task
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

# Load the required R module (if needed)
module load R

# Run the R script with the SLURM_ARRAY_TASK_ID as an argument
Rscript sim_sample_size/02_fit_models.R $SLURM_ARRAY_TASK_ID
