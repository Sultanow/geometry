#!/bin/bash
# one-time setup on the head node: conda env with numpy, shapely, networkx, python-sat; then the controls must pass
set -euo pipefail
cd "$(dirname "$0")/.."; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh
[ -d "$FN_ENV" ] || conda create -y -p "$FN_ENV" python=3.12 numpy shapely networkx pip
conda activate "$FN_ENV"
pip install -q python-sat
mkdir -p "$FN_OUT" "$FN_LOGS"
python tests/controls.py
