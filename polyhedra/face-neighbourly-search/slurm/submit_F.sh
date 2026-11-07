#!/bin/bash
# Run F: growth. Start from Mizhaev's and Roest-Vigh's polyhedra with one corner cut off by a ninth plane (a genus-3 surface
# with nine faces), then climb in the genus-aware score: only genus >= 3 surfaces count, maximise the pairs of faces that
# share an edge. 36/36 would be the first face-neighbourly polyhedron with nine faces. 128 tasks x 12 h (one node).
cd "$(dirname "$0")/.."; source slurm/config.sh
sub() { sbatch --parsable --job-name=fnF-$1 --array=0-$(($2 - 1)) --time=12:10:00 \
        --export=ALL,FAMILY=$1,BUDGET=43200,RUN=F,SCORER=genus,EVAL_TIMEOUT=60,CLIMB_ARGS="--step 0.05 --patience 8" slurm/climb.sbatch; }
for spec in "GROW_M9 64" "GROW_R9 64"; do set -- $spec; echo "$1: job $(sub $1 $2) ($2 tasks)"; done
