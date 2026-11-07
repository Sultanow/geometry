#!/bin/bash
# Run B: eleven faces in four families. 384 single-core tasks x 12 h on std.
cd "$(dirname "$0")/.."; source slurm/config.sh
sub() { sbatch --parsable --job-name=fnB-$1 --array=0-$(($2 - 1)) --time=12:10:00 --export=ALL,FAMILY=$1,BUDGET=43200,RUN=B slurm/search.sbatch; }
for spec in "S4_11 96" "D2_11 96" "C2_11 96" "ASYM_11 96"; do set -- $spec; echo "$1: job $(sub $1 $2) ($2 tasks)"; done
