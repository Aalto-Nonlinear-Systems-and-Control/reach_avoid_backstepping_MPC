#!/bin/bash
# Speed sweep of the study of program (14) solved directly (study14.py), kfac = 1 unless given as task:s:kfac.
cd "$(dirname "$0")"; PY=/Users/jianqiang/miniforge3/envs/rab_mpc/bin/python
for spec in "$@"; do
  IFS=: read task s k <<< "$spec"; k=${k:-1}
  TASK=$task LAM=$(cat data/lam_$task.txt) $PY study14.py $s $k 2>&1 | grep -v -i warn
done
echo "=== done $(date +%H:%M:%S)"
