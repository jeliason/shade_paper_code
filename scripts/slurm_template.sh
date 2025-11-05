#!/bin/bash
#SBATCH --job-name={{job_name}}
#SBATCH --output=logs/{{log_prefix}}_%A{% if array_range %}_{% endif %}%a.out
#SBATCH --error=logs/{{log_prefix}}_%A{% if array_range %}_{% endif %}%a.err
{% if array_range %}#SBATCH --array={{array_range}}
{% endif %}#SBATCH --time={{time}}
#SBATCH --mem={{mem}}
#SBATCH --cpus-per-task={{cpus}}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={{email}}

# Set environment variables for R scripts
export SYSTEM_ENV="HPC"
export HPC_DATA_PATH="{{data_path}}"
{% if force_rerun %}export FORCE_RERUN="1"
{% endif %}

# Load the required R module (if needed)
module load {{compiler_toolchain}}
module load {{r_module}}

# Run the R script
Rscript {{script_path}}{% if array_range %} $SLURM_ARRAY_TASK_ID{% endif %}
