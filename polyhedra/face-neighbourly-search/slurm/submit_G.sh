#!/bin/bash
# Run G: second generation. Start from the certified 34/36 near-misses of run F (analysis_F/analysis.jsonl) and climb in the
# three-stage score: genus 3, then every face >= 8 corners (the hexagon bottleneck), then all 36 pairs. 128 tasks x 12 h.
cd "$(dirname "$0")/.."; source slurm/config.sh
SEEDS="$FN_OUT/analysis_F/analysis.jsonl"; [ -s "$SEEDS" ] || { echo "missing $SEEDS - run the analysis of run F first"; exit 1; }
echo "seeds: $(grep -c '"pairs": 34' "$SEEDS") configurations with 34 pairs"
sbatch --parsable --job-name=fnG-SEED_9 --array=0-127 --time=12:10:00 \
  --export=ALL,FAMILY=SEED_9,BUDGET=43200,RUN=G,SCORER=faces,EVAL_TIMEOUT=90,FN_SEEDS="$SEEDS",FN_SEED_MIN=34,CLIMB_ARGS="--step 0.05 --patience 8 --hop 0.7" \
  slurm/climb.sbatch
