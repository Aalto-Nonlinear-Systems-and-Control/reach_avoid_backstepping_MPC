#!/bin/bash
# RA-MPC with the LP designs of study14_lp.py: task:s[:obj] ...
cd "$(dirname "$0")"; PY=/Users/jianqiang/miniforge3/envs/rab_mpc/bin/python; F="Warning|Ipopt|EPL|coin-or|\*\*\*\*|^$|For more"
for spec in "$@"; do
  IFS=: read task s o <<< "$spec"; o=${o:-coef}; echo "=== $task s=$s $o $(date +%H:%M:%S)"
  S14_DIR=study14lp S14_TAG=${task}_lp_s${s}_D8_$o TASK=$task LAM=$(cat data/lam_$task.txt) $PY eval14_mpc.py $s 0 2>&1 | grep -v -E "$F" | tail -1
done
echo "=== done $(date +%H:%M:%S)"
