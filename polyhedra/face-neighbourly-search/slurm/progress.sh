#!/bin/bash
# progress of the guided run D: best score per family (out of all pairs), final scores of the climbs, full hits
cd "$(dirname "$0")/.."; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
squeue -u "$USER" -h -o "%j %T" | sort | uniq -c
python progress.py --root "$FN_OUT/${1:-D}"
