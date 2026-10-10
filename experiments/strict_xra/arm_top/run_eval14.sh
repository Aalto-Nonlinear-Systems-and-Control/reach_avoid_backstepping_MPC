#!/bin/bash
# RA-MPC on the common initial sets with the terminal designs of study14.py: task:s[:kfac] ...
cd "$(dirname "$0")"; PY=/Users/jianqiang/miniforge3/envs/rab_mpc/bin/python; F="Warning|Ipopt|EPL|coin-or|\*\*\*\*|^$|For more"
for spec in "$@"; do
  IFS=: read task s k <<< "$spec"; k=${k:-1}; echo "=== $task s=$s k=$k $(date +%H:%M:%S)"
  TASK=$task LAM=$(cat data/lam_$task.txt) $PY eval14_mpc.py $s $k 2>&1 | grep -v -E "$F" | tail -1
done
echo "=== done $(date +%H:%M:%S)"
