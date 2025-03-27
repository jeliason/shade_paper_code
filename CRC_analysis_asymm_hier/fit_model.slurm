#!/bin/bash
#SBATCH --job-name=fit_CRC
#SBATCH --output=logs/fit_CRC_%A_%a.out
#SBATCH --error=logs/fit_CRC_%A_%a.err
#SBATCH --array=1-8
#SBATCH --time=20:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joelne@umich.edu

module load Rgeospatial

Rscript scripts/CRC_analysis_asymm_hier/fit_model_CRC.R $SLURM_ARRAY_TASK_ID
