# Panda-like arm in MuJoCo, kitting task (red, forward): first results

Setting of 2026-10-09, replacing the two-link arm by a Panda-like arm simulated in MuJoCo (scripts in `franka/`,
environment `rab_franka`, a clone of `rab_mpc` with MuJoCo 3.15; the earlier results are untouched).

## 1. Model and task

| item | value |
| --- | --- |
| robot | Franka Panda kinematics and meshes from MuJoCo Menagerie; joints 1, 2, 4, 6 actuated, joints 3, 5 fixed at 0, joint 7 at -pi/4, fingers fixed |
| dynamics (simplified) | lumped masses 6.8 kg at the elbow, 2.9 kg at the wrist, 1.5 kg at the tool centre point, yaw inertia 0.05 kg m^2 of link 1, armature 0.1 per joint; no gravity (compensated, as in Franka's torque interface), no damping |
| simulator | MuJoCo with the same parameters (`franka/model/panda_simple.xml`), RK4 at 0.5 ms, torque held over T; closed-form CasADi model agrees with MuJoCo to 1e-8 |
| outputs | tool centre point (y1, y2, y3) in metres and hand pitch y4 = -q2 + q4 + q6 (A is 4 x 4, theory unchanged) |
| workcell | the square kitting layout scaled by 0.10 m per unit, transport height 0.30 m above the base; all corridor points of the four tasks reachable with the hand pointing down, at least 0.19 rad from every joint limit |
| safe set | psi = psi_xy(y1/0.1, y2/0.1) - (y3 - 0.30)^2/0.03^2 - y4^2/0.1^2 (3 cm height and 0.1 rad pitch tolerance), psi_xy and k1 of the red task reused (exact spatial scaling), k1 = (0.1 k1_xy, -(y3 - 0.30), -y4) |
| target | disc of radius 3.5 cm around the slot in the plane |

## 2. Certificate, scale and bounds

| item | value |
| --- | --- |
| first layer | margin on the planar grid +1.4e-3; height and pitch terms add (2 + lambda) c (.)^2 >= 0 |
| mu0, lambda | 0.001 (counterpart of the arm's 0.1), 0.1 |
| period threshold | T_eps(1) = 41.3 ms from 40 000 terminal-set samples (4000 samples gave 76.9 ms and missed a rare state: more samples are needed in four outputs) |
| scale | s = 0.8 T_eps(1)/T = 3.31 at T = 10 ms; T_eps(s) = 12.5 ms = T_eps(1)/s (exact time scaling) |
| task torque bounds | 2.6, 4.6, 2.2, 0.4 N m (1.25 s^2 max|tau_1| on the terminal set; 2.5 to 5 % of the hardware limits 87/87/87/12) |
| URA | same k1 at s = 6.61 |
| MPC | N = 60 (0.6 s), 2 RK4 sub-steps, state constraint at sub-step nodes, compiled NLP (conda clang, CasADi JIT) with the active horizon as a parameter |
| initial set | 150 states, corridor positions at the transport height with the hand down, planar velocities around the URA field (r0 = 0.3), Dubins convention |

## 3. Results (150 states)

| | reach the slot | remark |
| --- | --- | --- |
| URA (s = 6.61) | 149 | exceeds the bounds in 144 successes, above the bounds 74 % of the run (median), peaks (11.4, 20.1, 6.6, 1.6) N m |
| vanilla MPC | 111 | all 111 feasible starts succeed |
| RA-MPC | 123 | 27 infeasible at k = 0, none afterwards |
| terminal controller alone | 96 continuous, 94 sampled | |

Pairwise: RA-MPC and not vanilla 15, vanilla and not RA-MPC 3 (feasible for vanilla at k = 0 but not for RA-MPC: the terminal set also asks for the velocity of the field), RA-MPC and not URA 0, URA and not RA-MPC 26.
Common successes (123): time to the slot median 0.44 s (max 1.08) for RA-MPC against 0.40 s (max 1.53) for the URA; cost 5.47 against 6.25.

Theorem 2 / Corollary 1: 5530 solves, no infeasible solve after k = 0, decrease >= stage cost in every step (smallest ratio 0.99993), horizon nonincreasing, smallest psi along the closed loops 0.0049, entry-time bound holds 123/123.
Lemma 2 (8000 independent terminal-set samples): one-step increase at every sample up to 15 ms (first failure at 20 ms); from the terminal set all 150 sampled-data closed loops reach the slot with Psi increasing up to 20 ms; second-order bound 4.4 ms.
Solve time with the compiled NLP: median about 0.95 s per step under parallel load (offline; the simulation waits for the solver).

## 4. Figures (`franka/figs/red_fwd_sq/`)

`fig_franka_vanilla`, `fig_franka_uncon`, `fig_franka_rampc` (+ `_far`): top views in the Dubins style, the arm drawn as its top-view footprint (with q3 = q5 = 0 all links project on the ray from the base to the gripper);
`fig_franka_compare` (+ `_far`): torques of URA and RA-MPC against the bounds; `verify_outcomes`.

## 5. Open points

- The gap between RA-MPC and vanilla (123 against 111) is smaller than for the two-link arm (137 against 84): the 0.6 s horizon covers most of the motion. A shorter horizon (N = 40) should widen it; to be calibrated if wanted.
- MuJoCo video of the forward and reverse cycle, and the other seven certificates.
- Reproduction: `python build_model_simple.py`, `python franka_model.py` (model check), `python solve_scale_franka.py`, `python run_main_franka.py`, `python run_lemma2_franka.py`, `python summarize_franka.py`, `python make_figs_franka.py` (environment `rab_franka`).
