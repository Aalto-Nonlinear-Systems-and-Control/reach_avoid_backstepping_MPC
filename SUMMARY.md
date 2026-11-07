# Summary

Code and data of the article *Reach-Avoid Model Predictive Control with
Guaranteed Recursive Feasibility via Input Constrained Backstepping* by
Jianqiang Ding, Nishant Jayesh Bhave and Shankar A. Deka (submitted to the
IEEE Open Journal of Control Systems; preprint arXiv:2604.03407).

## Intent

The article proposes a sampled-data model predictive controller (RA-MPC) for
control-affine nonlinear systems that must reach a target set while staying
inside a safe set under input constraints. Its terminal set is a subset of an
input-constrained reach-avoid set, built by backstepping from an exponential
control guidance-barrier function (ECGBF), and its terminal cost is built from
the same function. The code reproduces every certificate, simulation, table,
figure and video of the article, and the committed data allow every figure and
number to be regenerated without rerunning the experiments.

## Components

| Component | Location | Content |
| --- | --- | --- |
| Certificate synthesis (MATLAB) | `experiments/strict_xra/solve_vanilla_k1_*.m`, `experiments/strict_xra/design/`, `matlab/` | sum-of-squares programs (12) and (14) of the article for the virtual controller `k1` and the reach-avoid set, solved with SOSTOOLS and MOSEK; helpers between symbolic and polynomial objects |
| Controller construction (Python) | `experiments/strict_xra/design/cert.py`, `ra_mpc_cert.py`, `python/ra_mpc_dubins.py`, `controllers/` | ECGBF `Psi`, backstepping controller `u^c`, time scaling, Algorithm 1, vanilla MPC and the shared CasADi/IPOPT MPC pieces |
| Example 1, Dubins car | `experiments/strict_xra/final/` | three controllers on 300 common initial states, sampled-data checks of Lemma 2, success rate over the sampling period, timings, figures, `summary.txt`; notes in `PAPER_DUBINS.md` |
| Example 2, SCARA manipulator | `experiments/strict_xra/arm/`, `arm_top/` | workcell layout, eight kitting tasks under common torque bounds, the same experiments as Example 1 for each task; notes in `PAPER_ARM_TOP.md` |
| MuJoCo scene and videos | `experiments/strict_xra/scara/`, `docs/media/` | SCARA scene, Fig. 6 of the article, kitting-cycle and red-task videos; notes in `PAPER_SCARA.md` |
| Documentation | `README.md`, `experiments/strict_xra/README.md` | requirements, reproduction steps and the folder index |

Material that is not needed to reproduce the article, such as a Panda-like
arm in MuJoCo and the design notes behind the experiments, is on the branch
`journal-full`; the studies of the earlier letter version are on
`tac-revision-archive`.

## Reproduction

See `README.md` for the software versions, the conda environments
(`environment.yml`, `environment-mujoco.yml`) and the commands that reproduce
each example in order, with run times on an Apple M5.

## Licence

The code and data are released under the MIT License (`LICENSE`). The
dependencies keep their own licences, listed in `README.md`.
