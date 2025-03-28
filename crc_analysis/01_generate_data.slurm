#!/bin/bash
#SBATCH --job-name=make_large_data
#SBATCH --output=logs/make_large_data_%A_%a.out
#SBATCH --error=logs/make_large_data_%A_%a.err
#SBATCH --array=1-8
#SBATCH --time=10:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

module load R

Rscript crc_analysis/01_generate_data.R $SLURM_ARRAY_TASK_ID
