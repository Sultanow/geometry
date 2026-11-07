#!/bin/bash
# Run E: guided search with the corrected cell construction and the anytime score, hops near the best point.
# Focus on 9 and 10 faces with little or no symmetry; the perturbed families redone. 128 single-core tasks x 12 h (one node).
cd "$(dirname "$0")/.."; source slurm/config.sh
sub() { sbatch --parsable --job-name=fnE-$1 --array=0-$(($2 - 1)) --time=12:10:00 --export=ALL,FAMILY=$1,BUDGET=43200,RUN=E slurm/climb.sbatch; }
for spec in "C2_9 24" "ASYM_9 24" "C2_10 24" "ASYM_10 16" "S4_10 8" "D2_10 8" "S4_11 8" "D2_11 8" "C2_11 8"; do
  set -- $spec; echo "$1: job $(sub $1 $2) ($2 tasks)"; done
