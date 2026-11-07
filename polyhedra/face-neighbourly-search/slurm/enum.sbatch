#!/bin/bash
#SBATCH --job-name=fn-enum
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --array=0-255
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%A_%a.out
# complete enumeration of the combinatorial candidates (nine faces, genus 3, three extra edge copies), one cube per task
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
export OMP_NUM_THREADS=1
python enumerate_run.py --cubes "$SLURM_ARRAY_TASK_ID" --out "$FN_OUT/enum_g3"
