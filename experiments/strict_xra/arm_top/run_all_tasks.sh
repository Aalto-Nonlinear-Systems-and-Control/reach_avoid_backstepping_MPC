#!/bin/bash
# Main experiments of the other kitting tasks (after run_prep_tasks.sh), as for red_fwd_sq but without the horizon
# calibration and the period sweep: comparison on 150 states with Theorem 2 checks, Lemma 2, timing, summary, figures.
# Usage: ./run_all_tasks.sh blue_fwd green_fwd ...  (environment rab_mpc; outputs in data/<task>_sq/ and figs/<task>_sq/)
cd "$(dirname "$0")"; PY=${PY:-python}; F="Warning|Ipopt|EPL|coin-or|\*\*\*\*|^$|For more"
for T in "$@"; do
  task=${T}_sq; [ -f data/lam_$task.txt ] || { echo "no prep for $task"; continue; }; export TASK=$task LAM=$(cat data/lam_$task.txt); mkdir -p data/$task
  echo "=== $task (LAM $LAM) $(date +%H:%M:%S)"
  $PY run_main_top.py 2>&1 | grep -v -E "$F" > data/$task/run_main.out; tail -2 data/$task/run_main.out | cut -c1-300
  $PY run_lemma2_top.py 2>&1 | grep -v -E "$F" > data/$task/run_lemma2.out; tail -1 data/$task/run_lemma2.out
  $PY run_timing_top.py 2>&1 | grep -v -E "$F" > data/$task/run_timing.out
  $PY summarize_top.py > data/$task/summary.txt 2>&1
  $PY make_figs_top.py 2>&1 | grep -v -E "$F" > data/$task/make_figs.out; grep -c saved data/$task/make_figs.out
done
echo "=== done $(date +%H:%M:%S)"
