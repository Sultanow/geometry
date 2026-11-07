#!/bin/bash
#SBATCH --job-name=fn-climb
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%A_%a.out
# one guided-search worker per array task; FAMILY, BUDGET (s), RUN, optional SCORER, EVAL_TIMEOUT, CLIMB_ARGS
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python climb.py --family "$FAMILY" --seed "$SLURM_ARRAY_TASK_ID" --seconds "$BUDGET" --out "$FN_OUT/$RUN/$FAMILY" \
       --scorer "${SCORER:-anytime}" --eval-timeout "${EVAL_TIMEOUT:-30}" ${CLIMB_ARGS:-}
