#!/bin/bash
#SBATCH --job-name=fn-search
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%A_%a.out
# one worker per array task; submitted by slurm/submit_A.sh and slurm/submit_B.sh with FAMILY, BUDGET (s) and RUN
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1      # one core per task, no BLAS oversubscription
python search.py --family "$FAMILY" --seed "$SLURM_ARRAY_TASK_ID" --seconds "$BUDGET" --out "$FN_OUT/$RUN/$FAMILY"
