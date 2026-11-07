#!/bin/bash
#SBATCH --job-name=fn-enumcollect
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%j.out
# merge the cubes and cross-check against the main case: sbatch slurm/enum_collect.sbatch  ->  $FN_OUT/enum_g3/summary.md
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
python enumerate_collect.py --root "$FN_OUT/enum_g3" --crosscheck
