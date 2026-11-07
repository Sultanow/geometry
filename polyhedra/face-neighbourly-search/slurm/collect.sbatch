#!/bin/bash
#SBATCH --job-name=fn-collect
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%j.out
# final report of a run: sbatch --export=ALL,RUN=A slurm/collect.sbatch   ->  $FN_OUT/report_A/summary.md
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python collect.py --root "$FN_OUT/$RUN" --out "$FN_OUT/report_$RUN"
