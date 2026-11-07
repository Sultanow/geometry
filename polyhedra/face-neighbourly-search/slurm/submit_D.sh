#!/bin/bash
# Run D: guided search (MaxSAT score + hill climbing). Calibration with eight faces, then nine, ten and eleven faces.
# 256 single-core tasks x 12 h on std (two nodes).
cd "$(dirname "$0")/.."; source slurm/config.sh
sub() { sbatch --parsable --job-name=fnD-$1 --array=0-$(($2 - 1)) --time=12:10:00 --export=ALL,FAMILY=$1,BUDGET=43200,RUN=D slurm/climb.sbatch; }
for spec in "S4_8 8" "D2_8 8" "S4_9 24" "D2_9 24" "C2_9 24" "S4_10 24" "D2_10 24" "C2_10 24" "S4_11 32" "D2_11 32" "C2_11 32"; do
  set -- $spec; echo "$1: job $(sub $1 $2) ($2 tasks)"; done
