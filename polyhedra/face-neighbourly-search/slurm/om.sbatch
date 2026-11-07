#!/bin/bash
#SBATCH --job-name=fn-om
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=21:00:00
#SBATCH --array=0-202
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%A_%a.out
# oriented-matroid SAT test of the 203 candidates (one per task), 20 h limit per candidate
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
python om_run.py --candidates "$FN_OUT/enum_g3/candidates.json" --ids "$SLURM_ARRAY_TASK_ID" --out "$FN_OUT/om_g3" --limit 72000
