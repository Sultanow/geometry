#!/bin/bash
# Run C: the ladder between eight and eleven faces (9 and 10 faces), 6 families x 32 single-core tasks x 4 h on std.
cd "$(dirname "$0")/.."; source slurm/config.sh
sub() { sbatch --parsable --job-name=fnC-$1 --array=0-$(($2 - 1)) --time=04:10:00 --export=ALL,FAMILY=$1,BUDGET=14400,RUN=C slurm/search.sbatch; }
for spec in "S4_9 32" "D2_9 32" "C2_9 32" "S4_10 32" "D2_10 32" "C2_10 32"; do set -- $spec; echo "$1: job $(sub $1 $2) ($2 tasks)"; done
