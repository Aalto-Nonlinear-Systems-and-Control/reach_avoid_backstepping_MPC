# Summary

Code and data of the article *Reach-Avoid Model Predictive Control with
Guaranteed Recursive Feasibility via Input Constrained Backstepping* by
Jianqiang Ding, Nishant Jayesh Bhave and Shankar A. Deka.

## Intent

The code implements the reach-avoid model predictive controller (RA-MPC) of
the article, whose terminal set is an input-constrained reach-avoid set built
by backstepping from an exponential control guidance-barrier function, and
reproduces every certificate, simulation, table and figure of the article.
The committed data allow every figure and number to be regenerated without
rerunning the experiments.

## Components

| Component | Location | Content |
| --- | --- | --- |
| Certificate synthesis (MATLAB) | `experiments/sos/`, `experiments/*/run_design.m` | sum-of-squares program for the virtual controller `k1`, solved with SOSTOOLS and MOSEK |
| Example 1, Dubins car | `experiments/dubins/` | certificate `Psi` and backstepping controller `u^c`, Algorithm 1 and vanilla MPC, the experiments behind Figs. 4 and 5 and Table I, their data and figures |
| Example 2, SCARA arm | `experiments/scara/` | the same for Figs. 7 and 8 and Table I, with the task in `task.json` |
| MuJoCo scene | `experiments/scara/scene/` | SCARA scene and the transport by RA-MPC behind Fig. 6 |
| Documentation | `README.md`, `docs/` | requirements, reproduction steps, repository layout, source of the README figure |

## Reproduction

See `README.md` for the software versions, the conda environments
(`environment.yml`, `environment-mujoco.yml`) and the commands that reproduce
each example in order, with run times on an Apple M5.

## Licence

The code and data are released under the MIT License (`LICENSE`). The
dependencies keep their own licences, listed in `README.md`.
