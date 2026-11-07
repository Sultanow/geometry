#!/bin/bash
# Run A: eight faces in four families, plus the Szilassi control with seven faces. 384 single-core tasks x 4 h on std.
cd "$(dirname "$0")/.."; source slurm/config.sh
sub() { sbatch --parsable --job-name=fnA-$1 --array=0-$(($2 - 1)) --time=04:10:00 --export=ALL,FAMILY=$1,BUDGET=14400,RUN=A slurm/search.sbatch; }
for spec in "C2_7 8" "S4_8 32" "D2_8 96" "C2_8 128" "ASYM_8 120"; do set -- $spec; echo "$1: job $(sub $1 $2) ($2 tasks)"; done
