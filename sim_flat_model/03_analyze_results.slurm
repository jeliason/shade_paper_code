#!/bin/bash
#SBATCH --job-name=analysis_simulated_no_hier        # Job name
#SBATCH --output=logs/analysis_no_hier_task_%A_%a.out   # Standard output log
#SBATCH --error=logs/analysis_no_hier_task_%A_%a.err    # Standard error log
#SBATCH --time=30:00                # Time limit hh:mm:ss
#SBATCH --mem=32G                       # Memory per task
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

# Load the required R module (if needed)
module load R

# Run the R script with the SLURM_ARRAY_TASK_ID as an argument
Rscript sim_flat_model/03_analyze_results.R