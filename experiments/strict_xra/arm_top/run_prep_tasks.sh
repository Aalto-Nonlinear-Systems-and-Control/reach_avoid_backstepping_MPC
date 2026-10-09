#!/bin/bash
# Preparation of the other kitting tasks (square layout): task spec, k1 by SOS (MATLAB), mu0 study, scale with the common
# torque bounds of the red task. Usage: ./run_prep_tasks.sh blue_fwd green_fwd ...  (environment rab_mpc; logs in data/)
# lambda: 0.1 as for red_fwd_sq, else the largest design lambda that succeeded (written to data/lam_<task>.txt).
cd "$(dirname "$0")"; PY=/Users/jianqiang/miniforge3/envs/rab_mpc/bin/python; U_FIX=${U_FIX:-28.5,8.5}
for T in "$@"; do
  name=${T%_*}; dir=${T#*_}; task=${T}_sq; echo "=== $task $(date +%H:%M:%S)"
  LAYOUT=square BRIDGE=0.8 $PY make_task.py $name $dir > data/make_task_$task.out 2>&1 || { echo "make_task failed"; cat data/make_task_$task.out; continue; }
  [ -f k1/${task}_d8_l0.1.txt ] || matlab -batch "TASK='$task'; DEG=8; DESIGN_LAMS=[0.05 0.0667 0.08 0.1]; run('run_design_top.m')" > data/design_$task.out 2>&1 < /dev/null
  grep -E "__DESIGN__|__ERR__" data/design_$task.out
  LAM=""; for l in 0.1 0.08 0.0667 0.05; do [ -f k1/${task}_d8_l$l.txt ] && { LAM=$l; break; }; done
  [ -z "$LAM" ] && { echo "no k1 for $task"; continue; }; echo $LAM > data/lam_$task.txt
  TASK=$task LAM=$LAM $PY scale_study_top.py > data/scale_study_$task.out 2>&1; grep -E "^k1|mu0 0.1 " data/scale_study_$task.out | cut -c1-220
  TASK=$task LAM=$LAM U_FIX=$U_FIX $PY solve_scale_top.py > data/solve_scale_$task.out 2>&1; grep -v -i warn data/solve_scale_$task.out | cut -c1-400
done
echo "=== done $(date +%H:%M:%S)"
