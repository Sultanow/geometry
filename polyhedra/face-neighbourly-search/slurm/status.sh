#!/bin/bash
# quick status of a run (A or B): queue, and the summary table without the exact step
cd "$(dirname "$0")/.."; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
squeue -u "$USER" -o "%.12i %.14j %.8T %.10M %.6D %R" | head -n 20
python collect.py --root "$FN_OUT/${1:-A}" --out "$FN_OUT/status_${1:-A}" --no-exact
