#!/bin/bash
#SBATCH --array=1-50
#SBATCH -t 4:00:00
#SBATCH --account=def-liang-ab
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=niwanthimuthukuda@gmail.com
#SBATCH --mail-type=ALL

module load julia

julia --threads=1 CAM_50.jl $SLURM_ARRAY_TASK_ID