#!/bin/bash
# Speed sweep of study C (study14_lp.py) with the objective in OBJ (default coef): task:s ...
cd "$(dirname "$0")"; PY=/Users/jianqiang/miniforge3/envs/rab_mpc/bin/python
for spec in "$@"; do IFS=: read task s <<< "$spec"; OBJ=${OBJ:-coef} TASK=$task LAM=$(cat data/lam_$task.txt) $PY study14_lp.py $s 2>&1 | grep -v -i warn | grep -v " round "; done
echo "=== done $(date +%H:%M:%S)"
