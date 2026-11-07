#!/bin/bash
#SBATCH --job-name=fn-omcert
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --array=0-202
#SBATCH --output=/home/esul01/data/outputs/logs/%x_%A_%a.out
# certified decision of the 203 candidates: v2 encoding with DRAT proof checked by drat-trim, alternative encoding as cross-check
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; source slurm/config.sh
source /home/jupyter/anaconda3/etc/profile.d/conda.sh; conda activate "$FN_ENV"
python om_cert.py --candidates "$FN_OUT/enum_g3/candidates.json" --ids "$SLURM_ARRAY_TASK_ID" --out "$FN_OUT/om_cert_g3" \
       --drat-trim "$SLURM_SUBMIT_DIR/tools/drat-trim/drat-trim" --tmp "${TMPDIR:-/tmp}"
