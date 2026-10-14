# Reach-Avoid Backstepping MPC

<img src="reach_avoid_mpc.png" width="900" alt="Reach-avoid MPC: method (left) and closed loop (right)"/>

Code and data behind the journal version of *Reach-Avoid Model Predictive
Control with Guaranteed Recursive Feasibility via Input Constrained
Backstepping* (J. Ding, N. J. Bhave, S. A. Deka). This branch,
`tac-revision`, contains only what the paper uses. Earlier code is kept on
other branches, see [Branches](#branches).

## Method in brief

For a control-affine system with output `y = h(x)`, a safe set
`{psi(y) > 0}`, a target set `{phi(y) < 0}` and box input bounds:

1. **Offline.** An SOS program gives the virtual controller `k1(y)`, with the
   rate `lambda` fixed and `|k1_i| <= K` for a numerically robust solution.
   Backstepping turns it into the reach-avoid certificate `Psi` and the
   continuous controller `u^c`. The design is then scaled in time so that
   `u^c` respects the input bounds and `Psi` increases over one held sampling
   period at sampled states of the terminal set `X_RA^delta = {Psi >= delta}`.
2. **Online.** A sampled-data MPC with terminal constraint `X_RA^delta` and
   terminal cost `omega (M_Psi - Psi)` steers the system. The shifted
   solution completed by `u^c` is always feasible, so the MPC stays
   recursively feasible and reaches the target. The horizon shrinks by one
   step when that candidate already enters the target.

## Repository layout

```
.
├── README.md, reach_avoid_mpc.png
├── docs/readme_figure/         LaTeX source of reach_avoid_mpc.png (block diagram and Fig. 1)
├── matlab/                     SOS helpers: poly2sym.m, sym2pvar.m
├── python/ra_mpc_dubins.py     Dubins model, vanilla MPC and shared MPC pieces
├── controllers/                L-CSS export imported by ra_mpc_dubins.py (see Notes)
└── experiments/strict_xra/
    ├── solve_vanilla_k1_controller.m, solve_vanilla_k1_core.m
    │                           SOS program for k1 (lambda fixed, |k1_i| <= K, scaled outputs)
    ├── design/                 Dubins certificate: run_design.m, cert.py (Psi, u^c),
    │                           ra_mpc_cert.py (Algorithm 1), k1/W3_l0.06.txt
    ├── final/                  Example 1, Dubins car: config.py, solve_scale.py,
    │                           run_*.py, make_figs.py, summarize.py, data/, figs/
    ├── arm/                    workcell layout of Example 2 (screen_racks.py and helpers)
    ├── arm_top/                Example 2, SCARA arm, all eight kitting tasks:
    │                           tasks/, k1/, data/<task>/, figs/<task>/
    ├── scara/                  MuJoCo scene, stills and videos of the SCARA
    ├── franka/                 Panda-like seven-DoF arm in MuJoCo (supplementary example)
    └── PAPER_DUBINS.md, PAPER_ARM_TOP.md, PAPER_SCARA.md, PAPER_FRANKA.md
                                notes with every number, table and check behind the paper;
                                REPORT.md and STUDY_14_DIRECT.md record the design decisions
```

## Paper items and where they come from

| Paper item | Script | Output |
| --- | --- | --- |
| Example 1, Fig. 4 | `final/make_figs.py` | `final/figs/fig_dubins_{vanilla,uncon,rampc}.png` |
| Example 1, Fig. 5 | `final/make_figs.py` | `final/figs/fig_dubins_compare.png` |
| Example 2, Fig. 6 | `scara/paper_still.py` | `scara/figs/fig_scara_scene.png` |
| Example 2, Fig. 7 | `arm_top/make_figs_top.py` | `arm_top/figs/red_fwd_sq/fig_manip_{vanilla,uncon,rampc}.png` |
| Example 2, Fig. 8 | `arm_top/make_figs_top.py` | `arm_top/figs/red_fwd_sq/fig_manip_compare.png` |
| Table I | `final/run_timing.py`, `arm_top/run_timing_top.py`, the two `summarize*.py` | `final/data/summary.txt`, `arm_top/data/red_fwd_sq/summary.txt` |
| Time scaling, both examples | `final/solve_scale.py`, `arm_top/solve_scale_top.py` | `final/data/scale.json`, `arm_top/data/scale_red_fwd_sq.json` |
| Panda-like arm | `franka/*.py` | `franka/data/red_fwd_sq/`, `franka/figs/red_fwd_sq/` |

All data behind these items are committed, so every figure and summary can be
redrawn without rerunning the experiments: run `make_figs*.py` and
`summarize*.py` alone.

## Requirements

| Part | Software (versions used) |
| --- | --- |
| Certificate synthesis | MATLAB R2025a with the Symbolic Math Toolbox, [SOSTOOLS](https://github.com/oxfordcontrol/SOSTOOLS) 4.03 and MOSEK 11.0 on the MATLAB path |
| Experiments | Python 3.11 with numpy 2.4, scipy 1.17, sympy 1.14, casadi 3.7.2 (IPOPT), joblib 1.5 and matplotlib 3.10; conda environment `rab_mpc` |
| MuJoCo scenes and videos | the packages above plus mujoco 3.15 and imageio; conda environment `rab_franka` |
| Panda meshes | [MuJoCo Menagerie](https://github.com/google-deepmind/mujoco_menagerie) cloned into `experiments/strict_xra/franka/third_party/mujoco_menagerie` |

```bash
conda create -n rab_mpc python=3.11 numpy scipy sympy matplotlib joblib -c conda-forge
conda activate rab_mpc && pip install casadi
```

## Reproduction

Run times are for an Apple M5 with 10 cores; the experiments use all cores
through joblib, except the timing runs, which run alone.

### Example 1, Dubins car (`experiments/strict_xra/`)

```bash
cd design
matlab -batch "DESIGN_TAGS={'W3'}; DESIGN_LAMS=0.06; run('run_design.m')"   # k1/W3_l0.06.txt
cd ../final
python solve_scale.py          # time scale s from the input bounds and T -> data/scale.json
python verify_inputs_grid.py   # input bounds of u^c on a grid of X_RA^delta
python run_main.py             # three methods on 300 common initial states (3 min)
python run_lemma2.py           # sampled-data checks of the terminal controller (5 min)
python run_sweep.py            # success rate over T = 0.025 ... 0.625 s (15 min)
python run_sweep_ext.py        # adds T = 0.417 s and 1.25 s to the sweep
python run_timing.py           # computation times, run alone
python make_figs.py            # all figures -> figs/
python summarize.py            # all numbers, as in PAPER_DUBINS.md
```

### Example 2, SCARA arm (`experiments/strict_xra/`)

```bash
cd arm && LAYOUT=square python screen_racks.py 0.8 0.0     # workcell layout
cd ../arm_top
./run_prep_tasks.sh red_fwd    # task, k1 (MATLAB), mu study and time scale of red_fwd_sq
export TASK=red_fwd_sq LAM=0.1
python calib_horizon.py 30,45,60 30
python run_main_top.py && python run_lemma2_top.py
PERIODS_MS=5,7.5,10,12.5,15,20,25,30,40,50,60,75,100,150,200,300,600 python run_sweep_top.py 150
python run_timing_top.py
python summarize_top.py > data/red_fwd_sq/summary.txt && python make_figs_top.py
# the other seven kitting tasks, under the torque bounds of red_fwd_sq
./run_prep_tasks.sh red_rev blue_fwd blue_rev green_fwd green_rev orange_fwd orange_rev
./run_all_tasks.sh red_rev blue_fwd blue_rev green_fwd green_rev orange_fwd orange_rev
python summarize_all_tasks.py  # data/all_tasks_summary.{json,md}
```

The shell scripts use `$PY` as the interpreter (default `python`).

### MuJoCo scenes and videos (environment `rab_franka`)

```bash
cd experiments/strict_xra/scara
python build_scene_scara.py && python scara_demo_traj.py && python render_scara.py --stills
python paper_still.py          # figs/fig_scara_scene.png (Fig. 6)
python cycle_demo_traj.py && python render_cycle.py      # full kitting cycle video (16 min)
./view_scara.sh 0.5            # interactive replay window (macOS)
```

The replay scripts take the environment from `$MUJOCO_ENV`, else the active
conda environment. The Panda-like example follows the same pattern; its
steps are listed in `PAPER_FRANKA.md`.

## Branches

| Branch | Content |
| --- | --- |
| `tac-revision` | this branch, the code and data of the journal version |
| `tac-revision-archive` | snapshot of `tac-revision` before this cleanup: the studies of solving the input-constrained program directly (SOS with sampled bounds, sampled LP, regularisation), the alternative time scales of the Dubins example, the first SCARA layout, the early strict reach-avoid-set studies, the original plan `TAC_EXPERIMENT_UPDATE.md`, and the run logs |
| `lcss-revision`, `main` | the letter version: MATLAB synthesis of bounded controllers, exported controllers and notebooks |
| `mujoco-franka-reachavoid` | first MuJoCo experiments with the Franka arm |

## Notes

- `controllers/sop_bounded_control_dubins_car_result.py` is the controller
  exported for the letter version. `python/ra_mpc_dubins.py` imports it when
  loaded, but the journal pipeline does not use it. Do not regenerate it; the
  MATLAB scripts that would overwrite it live on the letter branches.
- `*.png` is ignored by `.gitignore`; the figures listed above are added with
  `git add -f`.
- `config.py` of each example holds its setting. Environment variables such
  as `TASK`, `LAM`, `SCALE` and `TAG` select a task or run a variant into its
  own `data<TAG>/` and `figs<TAG>/` folders.
