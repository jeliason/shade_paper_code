#!/bin/bash
#SBATCH --job-name=make_dsize_data         # Job name
#SBATCH --output=logs/dsize_task_%A_%a.out   # Standard output log
#SBATCH --error=logs/dsize_task_%A_%a.err    # Standard error log
#SBATCH --array=1-135                   # Array range (adjust as needed)
#SBATCH --time=45:00                # Time limit hh:mm:ss
#SBATCH --mem=16G                       # Memory per task
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

# Load the required R module (if needed)
module load R

# Run the R script with the SLURM_ARRAY_TASK_ID as an argument
Rscript sim_sample_size/01_generate_data.R $SLURM_ARRAY_TASK_ID
