#!/bin/bash
#SBATCH --job-name=fn-analyze
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%j.out
# analysis of the best configurations of a guided run: sbatch --export=ALL,RUN=F,MIN=34 slurm/analyze.sbatch
# for run G (three-stage score): sbatch --export=ALL,RUN=G,MIN=34,SCORER=faces,TIMEOUT=180 slurm/analyze.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python analyze_best.py --root "$FN_OUT/$RUN" --out "$FN_OUT/analysis_$RUN" --min-score "${MIN:-34}" --workers 16 --scorer "${SCORER:-genus}" --timeout "${TIMEOUT:-120}"
