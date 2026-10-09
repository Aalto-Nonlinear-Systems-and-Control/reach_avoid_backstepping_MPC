# Top-view arm, kitting workcell: first task (red, forward)

Setting of 2026-10-09. The arm example moves to the horizontal plane ($g = 0$) and
to a kitting workcell: two vertical racks, four coloured objects on the left, four
slots on the right, a crossing assignment, a forward and a reverse phase. This
document covers the task with the largest span, the red object from the lowest
left compartment to the third right compartment, as agreed before running the
other seven. The workcell is the **compact (square) layout**, chosen so that the
paper panels are square like those of the Dubins car; the first, wide layout is
kept in Section 10. Scripts in `arm_top/` (Section 9, task `red_fwd_sq`); the
option-A results of `PAPER_ARM.md` are untouched. The paper and the notebooks are
not yet changed.

Contents: 1 setting, 2 three-method comparison, 3 Lemma 2, 4 Theorem 2 and
Corollary 1, 5 period sweep and timing, 6 figures, 7 observations, 8 places to
update and open points, 9 reproduction, 10 the wide layout.

## 1. Setting

### 1.1 Workcell and task

| item | value |
| --- | --- |
| arm | $m_i = 1$ kg, $l_i = 4$ m, $l_{c,i} = 2$ m, $I_i = 0.02$ kg m$^2$, horizontal plane, $g = 0$; output $y = h(q)$, branch $q_2 > 0$ |
| bench and racks | bench edge $y_2 = 0.8$, base at the origin; racks at $y_1 = \mp 2.4$, compartment heights $2.3, 3.3, 4.3, 5.3$ (pitch 1.0), 1.2 deep; the scene fits $[-3.7, 3.7] \times [-0.7, 6.7]$ |
| assignment | object $k$ to slot $\sigma(k)$, $\sigma = (3, 1, 4, 2)$: red $2.3 \to 4.3$, blue $3.3 \to 2.3$, green $4.3 \to 5.3$, orange $5.3 \to 3.3$; three crossings |
| fixture | box $\lvert y_1\rvert \le 0.5$, $1.2 \le y_2 \le 2.56$ between the racks |
| this task | red, forward: source $(-2.4, 2.3)$, slot $(2.4, 4.3)$, path length 5.5, rise 2.0 (largest span, tied with orange) |

All eight certificates of the layout (four tasks, two directions) pass the checks of Section 1.2 (`arm/data/screen_racks_square.out`, figure `arm/figs/arm_racks_square.png`).

### 1.2 Safe set and target

Reference path $y_2 = c(y_1)$ with $x = y_1/2.7$,
$c = m + d\,(3x - x^3)/2 + b\,(1 - x^2)^2$, bridge $b = 0.8$, $(m, d)$ such that
$c(\mp 2.4)$ equals the two compartment heights; the shape is flattest 0.3
inside the compartment, so the corridor stays level there.

$$\psi(y) = \frac{w(y_1)\,[h^2 - (y_2 - c(y_1))^2] - \varepsilon\,(y_1/2.4)^{10}}{\psi_{\max}},\qquad w = 1 + 0.0356\,(y_1 + 3.375)^2,$$

$h = 0.38$, $\varepsilon = \kappa (X_0 + a) h^2 X_0 / 5$ with $X_0 = 2.4$ (so that $\partial\psi/\partial y_1 = 0$ at the slot),
$\psi_{\max} = 0.3014$, degree 10. Target $\mathcal{X}_T = \{\lvert y - (2.4, 4.3)\rvert < 0.35\}$.
The weight and the closing term are those of the wide layout scaled with $X_0$ ($a \propto X_0$, $\kappa \propto X_0^{-2}$), so the shape of $\psi$ along the path is the same.

Properties: $\partial\psi/\partial y_2 = 0$ only on the path, the remaining condition depends on $y_1$ alone,
so $\psi$ has exactly one critical point, its maximum, exactly at the slot. $\psi(\text{source}) = 0.45\,\psi_{\max}$.
$\lvert y\rvert \in [2.90, 5.59]$ on $\mathcal{X}_S$ (no singularity), corridor inside its compartments (smallest wall margin 0.02),
0.62 clear of the fixture, 0.57 clear of the other objects and slots.

### 1.3 Certificate, scale and bounds

| item | value |
| --- | --- |
| $k_1$ | (12) at the normalised scale, outputs $y/4$, degree 8, $\lvert k_{1,i}\rvert \le 1$ on $\{\psi \ge 0\}$, $\lambda$ fixed; $\delta_s/\lambda = 0$ for every $\lambda \in \{0.05, 0.0667, 0.08, 0.1\}$, $\lambda = 0.1$ used |
| base design check | ECGBF margin on the grid $\ge +2.0\times 10^{-3}$, $\max\lvert k_{1,i}\rvert = 0.997$ |
| $\mu_0$ | 0.1 (Section 7.3) |
| exact time scaling | with $g = 0$: $\tau_s(D_s x) = s^2\tau_1(x)$ and $T_\varepsilon(s) = T_\varepsilon(1)/s$; checked: $T_\varepsilon(1) = 29.8$ ms, $T_\varepsilon(s) = 12.50$ ms $= T_\varepsilon(1)/s$ |
| scale | $s = \min(s_U, 0.8\,s_\varepsilon) = 0.8\,T_\varepsilon(1)/T = 2.39$ at $T = 10$ ms; $\lambda_s = 0.239$, $\mu_s = 0.569$ |
| torque bounds | $\lvert\tau_1\rvert \le 28.5$, $\lvert\tau_2\rvert \le 8.5$ N m ($1.25\,s^2\max\lvert\tau_1\rvert$ on $\mathcal{X}_{RA}^\delta$, rounded); $s_U = 2.67 > s$ |
| terminal controller on $\mathcal{X}_{RA}^\delta$ | $\max\lvert\tau\rvert = (22.8, 6.7)$ N m, inside the bounds |
| second-order bound | $\eta = 0.12$, $M_2 = 49$, $T^\star = 2\eta/M_2 = 4.9$ ms (reported, not used, as for the Dubins car) |
| MPC | $T = 10$ ms, $N = 60$ (0.6 s), 5 RK4 sub-steps, state constraint at the sub-step nodes with back-off 0.005, $\delta = 0.5\,M_\Psi$, $\omega = 1.2\sup\ell/\varepsilon$, $Q_y = 5$, $R_u = 10^{-3}$ |
| URA | the same $k_1$ at $s = 4.77$ (twice the terminal scale, four times the torque) |
| initial set | 150 states, positions uniform in $\mathcal{X}_S \setminus \mathcal{X}_T$, output velocities uniform in a disc of radius $0.3 \cdot 4.77$ around the URA's field $4.77\,k_1(y)$ (Dubins convention) |

Horizon calibration on 30 states (`data/calib_red_fwd_sq.out`): $N = 30$ vanilla 10 / RA-MPC 16 (1 state vanilla only),
$N = 45$ 11 / 22, $N = 60$ 17 / 26 (none vanilla only). $N = 60$ as in the wide layout.

## 2. Three-method comparison (150 states)

| | reach $\mathcal{X}_T$ | remark |
| --- | --- | --- |
| URA ($s = 4.77$) | 149 | 1 unsafe; exceeds the bounds in 144 of the successes, peaks $(103, 28)$ N m, above the bounds for a median 51 % of the run |
| vanilla MPC | 84 | feasible at $k = 0$ from 127; **43 runs lose feasibility later** (median after 19 steps) |
| RA-MPC | 137 | 13 infeasible at $k = 0$, **none afterwards** |
| terminal controller alone, $s = 2.39$ | 96 continuous, 96 sampled | |

Pairwise: RA-MPC and not vanilla 53, **vanilla and not RA-MPC 0**, RA-MPC and not URA 0, URA and not RA-MPC 12, all three 84, none 1.
Five states are feasible for vanilla at $k = 0$ but not for RA-MPC (the terminal set also asks for the velocity of the field); vanilla loses feasibility later from all five.
Seven of the thirteen RA-MPC failures start within 1.2 of the slot with the URA's speed and cannot stop inside the corridor under the bounds; six start more than 3 away.

| | states | time to $\mathcal{X}_T$, median (max) | cost, median |
| --- | --- | --- | --- |
| URA vs RA-MPC | 137 | 0.63 (2.11) vs 0.65 (1.30) s | 10.8 vs 9.8 |
| all three | 84 | URA 0.39, RA-MPC 0.43, vanilla 0.43 s | 3.7, 3.7, 3.5 |

RA-MPC uses the bounds actively: $\tau_1$ within 2 % of its bound 41 % of the time, $\tau_2$ 84 %.
Smallest $\psi$ along the continuous closed loops before entering $\mathcal{X}_T$: RA-MPC 0.0054, vanilla 0.0050.

## 3. Verification of Lemma 2 (held torque, terminal set)

| $T$ [ms] | arcs | smallest one-step increase of $\Psi$ | increasing | end in terminal set | second-order bound |
| --- | --- | --- | --- | --- | --- |
| 1 | 6532 | $+1.2\times 10^{-4}$ | 100 % | 100 % | $+9.7\times 10^{-5}$ |
| 5 | 6506 | $+5.8\times 10^{-4}$ | 100 % | 100 % | $-1.6\times 10^{-6}$ |
| **10** | 6484 | $+6.2\times 10^{-4}$ | 100 % | 100 % | $-1.3\times 10^{-3}$ |
| 12.5 | 6475 | $+1.9\times 10^{-4}$ | 100 % | 100 % | $-2.4\times 10^{-3}$ |
| 15 | 6461 | $-6.9\times 10^{-4}$ | 99.85 % | 99.85 % | |
| 20 | 6439 | $-3.4\times 10^{-3}$ | 99.27 % | 99.41 % | |
| 50 | 6281 | $-4.8\times 10^{-2}$ | 96.86 % | 97.72 % | |

Sampled-data closed loop (400 states each): from $\mathcal{X}_{RA}^\delta$ all reach $\mathcal{X}_T$ safely at every period up to 30 ms, with $\Psi$ increasing at every sample up to 12.5 ms (398/400 at 15 ms); from $\mathcal{X}_{RA}$ all reach $\mathcal{X}_T$ at every period, median 0.67 s.
At $T = 10$ ms the lemma rests on the sampled check (threshold 12.5 ms); the second-order bound (4.9 ms) is conservative by a factor 2.5.

## 4. Verification of Theorem 2 and Corollary 1

137 runs feasible at $k = 0$, all reach $\mathcal{X}_T$; 8363 solves, no infeasible solve after $k = 0$; the shifted candidate was kept in 4896 steps;
$u^c$ at the predicted terminal states within the bounds in every step.
Cost decrease over 8226 steps: smallest 0.71, $J^*_k - J^*_{k+1} \ge \ell(x_k, u_k)$ in every step (smallest ratio 1.000000), never below $\ell_{\min} = 0.6125$.
Horizon nonincreasing in all runs. Corollary 1 bound holds in 137/137 runs (entry median 0.65 s, bound median 384 s: valid, loose because $\omega$ is large).

## 5. Period sweep and timing

40 states, $NT = 0.6$ s, terminal design ($\delta$, $\omega$) fixed:

| $T$ [ms] | $N$ | vanilla feasible (reached) | $u^c$ held | RA-MPC | later infeasible | decrease violations | min $\Delta\Psi$ on terminal set |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 5 | 120 | 37 (25) | 25 | 36 | 0 | 0/4118 | $+5.8\times 10^{-4}$ |
| 7.5 | 80 | 37 (25) | 25 | 36 | 0 | 0/2738 | $+6.9\times 10^{-4}$ |
| 10 | 60 | 37 (25) | 25 | 36 | 0 | 0/2050 | $+5.2\times 10^{-4}$ |
| 12.5 | 48 | 37 (25) | 24 | 36 | 0 | 0/1636 | $-1.3\times 10^{-5}$ |
| 15 | 40 | 37 (25) | 24 | 36 | 0 | 0/1361 | $-8.5\times 10^{-4}$ |
| 20 | 30 | 37 (25) | 24 | 36 | 0 | 0/1018 | $-3.4\times 10^{-3}$ |

The success rate is flat; from 12.5 ms the one-step increase fails on part of the terminal set (12.5 ms is the threshold itself), so the guarantee lapses although the closed loop still succeeds.

Timing (one run at a time, Apple M5): RA-MPC 156 ms per solve (first 425 ms), vanilla 123 ms, $u^c$ 0.24 ms. The MPC is about 16 times slower than the 10 ms period here (Section 8).

## 6. Figures (`arm_top/figs/red_fwd_sq/`)

| file | content |
| --- | --- |
| `fig_manip_vanilla_far`, `fig_manip_uncon_far`, `fig_manip_rampc_far` | square paper panels in the workcell: corridor, slot, successes and failures over the initial positions, trajectory from the farthest common start (URA 1.04 s, peaks (66, 24) N m; RA-MPC 1.27 s) |
| `fig_manip_compare_far` | $\rho$ against $T$, torques of URA and RA-MPC against the bounds |
| `fig_manip_cost_far` | $J^*_k$, decrease against $\ell$, horizon |
| same names without `_far` | shared start where all three succeed (state 32) |
| `verify_outcomes`, `verify_lemma2`, `verify_theorem2`, `verify_sweep` | verification panels of Sections 2 to 5 |
| `../../arm/figs/arm_racks_square.png` | the workcell with all four tasks in both phases |

## 7. Observations

1. The comparison has the intended shape: vanilla-reached is a subset of RA-MPC-reached, and the URA succeeds only by exceeding the bounds. In this layout RA-MPC has no success the URA lacks (the URA has 12 more); the order of the two is not fixed, as intended.
2. Vanilla MPC loses feasibility in 43 of its 127 initially feasible runs (no terminal ingredient); RA-MPC never does. A direct illustration of recursive feasibility, worth a sentence in the paper.
3. $\mu_0 = 0.1$ makes rest states mostly uncertified (3 % in $\mathcal{X}_{RA}$); $\mu_0 = 1$ certifies 55 % but $T_\varepsilon(1)$ drops to 8.6 ms and the runs become several times longer. With $\mu_0 = 0.1$ the initial set follows the Dubins convention. Rest starts remain possible for the demonstration cycle.

## 8. Places to update in the paper, and open points

- Example 2 statement: horizontal plane ($g = 0$), workcell and task description, $\psi$ of Section 1.2, target disc, bounds 28.5 / 8.5 N m, $T = 10$ ms, $N = 60$; the scale rule is the one of the Dubins car, $s = \min(s_U, 0.8 s_\varepsilon)$.
- Table I: re-measure on the paper's machine; here RA-MPC 156 ms per step against a 10 ms period, so the example is not real time on this hardware.
- Done (2026-10-09): the other seven certificates and their full comparison data (Section 11), and the full forward and reverse cycle on the SCARA (`PAPER_SCARA.md` Section 4); transit moves between tasks are outside the planar model and not certified.
- Open, carried over: second-order Lemma 2 in the paper (D7), notebooks not updated, `controllers/*.py` untouched.

## 9. Reproduction (`experiments/strict_xra/arm_top/`, environment `rab_mpc`)

```
(cd ../arm && LAYOUT=square python screen_racks.py 0.8 0.0)   # layout, eight certificates, racks_best_square.pkl
LAYOUT=square BRIDGE=0.8 python make_task.py red fwd            # tasks/red_fwd_sq.json
matlab -batch "TASK='red_fwd_sq'; DEG=8; DESIGN_LAMS=[0.05 0.0667 0.08 0.1]; run('run_design_top.m')"
export TASK=red_fwd_sq LAM=0.1
python scale_study_top.py && python solve_scale_top.py
python calib_horizon.py 30,45,60 30
python run_main_top.py && python run_lemma2_top.py && python run_sweep_top.py 40 && python run_timing_top.py
python summarize_top.py > data/red_fwd_sq/summary.txt && python make_figs_top.py
```

## 10. The wide layout (first run, kept for reference)

Racks at $y_1 = \mp 3.2$, heights $1.8$ to $4.8$, bridge 1.0, $\lambda = 0.05$, $s = 2.68$, bounds $(35.5, 8.5)$ N m, URA at $s = 5.36$, $N = 60$ (task `red_fwd`, data in `data/red_fwd/`, figures in `figs/red_fwd/`).
Counts of 150: URA 146 (142 exceed the bounds), vanilla 66 (95 feasible at $k = 0$, 29 lose feasibility), RA-MPC 139 (none after $k = 0$); vanilla-reached subset of RA-MPC-reached; RA-MPC and not URA 2, URA and not RA-MPC 9.
All Theorem 2 checks pass (10 405 solves); Lemma 2 sampled threshold 12.5 ms, second-order bound 3.4 ms; sweep flat at 38/40 from 5 to 20 ms; RA-MPC 185 ms per solve.
The panels are about 10.6 by 6.5 units wide, hence the square layout.

## 11. All eight tasks under common torque bounds (2026-10-09)

The red forward task above is the paper's example; the other seven tasks of the cycle (blue, green, orange forward and all four reverse tasks) were run with the same pipeline, so that the data cover the whole workcell and the SCARA cycle video (`PAPER_SCARA.md` Section 4) is backed by full results.
Settings as for red: square layout, k1 from the SOS (12) at lambda 0.1 with |k1_i| <= 1 and degree 8 (all seven designs succeeded at 0.1), mu0 = 0.1, delta = 0.5 M_Psi, T = 10 ms, N = 60, URA at twice the terminal scale, 150 initial states per task drawn as for red (positions uniform in the task's corridor, velocities around the URA field).
Not run for the seven new tasks: the horizon calibration and the period sweep.

**Common torque bounds.** One robot has one set of torque limits, so all tasks use red's bounds (28.5, 8.5) N m instead of bounds derived from each task's own terminal set (`U_FIX=28.5,8.5` in `solve_scale_top.py`; without it nothing changes, and red's numbers are untouched). Each task then takes s = min(s_eps, s_bnd), s_eps = 0.8 T_eps(1) / T from the period and s_bnd = min_i sqrt(U_i / (1.25 max|tau_1,i|)) from the bounds with the same 1.25 margin as before. For the three other forward tasks the bounds bind (their own terminal sets would ask for up to (49, 15) N m), so their terminal controllers are slower, s = 1.62 to 1.75 instead of 1.90 to 2.30; this only enlarges the period margin (T_eps(s) = 14 to 16 ms). For the reverse tasks the period binds, as for red.

### Design and scale

| task | lambda | s_eps | s_bnd | s | binding | T_eps(1) [ms] | T_eps(s) [ms] | second-order T* [ms] | own bounds would be | URA scale |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| red_fwd_sq | 0.1 | 2.386 | 2.387 | 2.386 | period | 29.8 | 12.5 | 5.07 | (28.5, 8.5) | 4.772 |
| blue_fwd_sq | 0.1 | 2.299 | 1.747 | 1.747 | bounds | 28.7 | 16.4 | 5.22 | (49.0, 15.0) | 3.495 |
| green_fwd_sq | 0.1 | 2.074 | 1.623 | 1.623 | bounds | 25.9 | 16.0 | 4.27 | (47.0, 7.0) | 3.246 |
| orange_fwd_sq | 0.1 | 1.901 | 1.690 | 1.690 | bounds | 23.8 | 14.1 | 4.67 | (36.5, 10.0) | 3.380 |
| red_rev_sq | 0.1 | 1.175 | 1.894 | 1.175 | period | 14.7 | 12.5 | 6.84 | (9.5, 3.5) | 2.350 |
| blue_rev_sq | 0.1 | 1.767 | 2.234 | 1.767 | period | 22.1 | 12.5 | 6.13 | (18.0, 5.0) | 3.535 |
| green_rev_sq | 0.1 | 1.499 | 1.970 | 1.499 | period | 18.7 | 12.5 | 7.09 | (17.0, 4.5) | 2.998 |
| orange_rev_sq | 0.1 | 1.653 | 1.829 | 1.653 | period | 20.7 | 12.5 | 5.63 | (23.5, 4.5) | 3.306 |

### Comparison

| task | URA (exceed bounds) | vanilla feasible k=0 / reached | RA-MPC (infeasible k=0, later) | RA not van | van not RA | RA not URA | URA not RA | median time URA / RA [s] | median cost URA / RA |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| red_fwd_sq | 149 (144) | 127 / 84 | 137 (13, 0) | 53 | 0 | 0 | 12 | 0.63 / 0.65 | 10.8 / 9.8 |
| blue_fwd_sq | 150 (139) | 61 / 61 | 137 (13, 0) | 77 | 1 | 0 | 13 | 0.94 / 0.94 | 19.6 / 17.8 |
| green_fwd_sq | 150 (126) | 94 / 83 | 144 (6, 0) | 61 | 0 | 0 | 6 | 0.99 / 0.77 | 17.4 / 12.2 |
| orange_fwd_sq | 150 (143) | 64 / 61 | 138 (12, 0) | 77 | 0 | 0 | 12 | 0.89 / 0.92 | 21.6 / 19.8 |
| red_rev_sq | 150 (20) | 71 / 71 | 150 (0, 0) | 79 | 0 | 0 | 0 | 1.18 / 0.83 | 28.0 / 16.6 |
| blue_rev_sq | 150 (128) | 92 / 92 | 142 (8, 0) | 50 | 0 | 0 | 8 | 0.66 / 0.68 | 7.5 / 7.1 |
| green_rev_sq | 150 (48) | 80 / 80 | 149 (1, 0) | 69 | 0 | 0 | 1 | 0.82 / 0.72 | 13.6 / 10.9 |
| orange_rev_sq | 150 (127) | 75 / 75 | 148 (2, 0) | 73 | 0 | 0 | 2 | 0.75 / 0.80 | 10.0 / 9.2 |

### Theorem 2 / Corollary 1 and Lemma 2

| task | solves | decrease violations | below ell_min | horizon nonincreasing | u^c in bounds | min psi RA | entry bound holds | Lemma 2: all increasing up to [ms] | T* [ms] | ZOH closed loop at 10 ms (terminal set, X_RA) | solve ms RA / vanilla |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| red_fwd_sq | 8363 | 0 | 0 | 137/137 | True | 0.0054 | 137/137 | 12.5 | 4.91 | 400/400, 400/400 | 173 / 123 |
| blue_fwd_sq | 11909 | 0 | 0 | 137/137 | True | 0.0045 | 137/137 | 15 | 5.09 | 400/400, 400/400 | 198 / 137 |
| green_fwd_sq | 10239 | 0 | 0 | 144/144 | True | 0.0055 | 144/144 | 15 | 4.54 | 400/400, 400/400 | 215 / 154 |
| orange_fwd_sq | 11846 | 0 | 0 | 138/138 | True | 0.0050 | 138/138 | 12.5 | 4.48 | 400/400, 400/400 | 213 / 105 |
| red_rev_sq | 11259 | 0 | 0 | 150/150 | True | 0.0050 | 150/150 | 12.5 | 6.93 | 400/400, 400/400 | 200 / 86 |
| blue_rev_sq | 9677 | 0 | 0 | 142/142 | True | 0.0090 | 142/142 | 12.5 | 6.23 | 400/400, 400/400 | 217 / 88 |
| green_rev_sq | 10050 | 0 | 0 | 149/149 | True | 0.0047 | 149/149 | 10 | 7.26 | 400/400, 400/400 | 217 / 97 |
| orange_rev_sq | 10663 | 0 | 0 | 148/148 | True | 0.0050 | 148/148 | 10 | 5.59 | 400/400, 400/400 | 200 / 78 |

### Observations

- Theorem 2 and Corollary 1 hold in every task: 84 006 solves, no infeasible solve after k = 0, cost decrease at least the stage cost at every step, horizon nonincreasing, u^c within the bounds at every candidate, smallest psi along the closed loops 0.0045 (back-off 0.005 at the sub-step nodes), entry-time bound satisfied in every run.
- Lemma 2 at T = 10 ms: Psi increases over one held period at every sampled terminal state in all tasks, and all 400 sampled-data closed loops from the terminal set and from X_RA reach the target safely. For green_rev and orange_rev 10 ms is the largest sampled period with no violation (the 0.8 factor of s_eps is the only margin there); the second-order bound gives 4.5 to 7.3 ms, below T, as for red.
- RA-MPC reaches 137 to 150 of 150 states in every task, vanilla MPC 61 to 92. The vanilla successes are a subset of the RA-MPC successes in seven tasks; in blue_fwd one state (121) is reached by vanilla in 0.26 s but is infeasible for RA-MPC at k = 0: it starts fast and close to the slot, and the terminal set also asks for the velocity of the field (the same effect as in the Franka example). RA-MPC succeeds wherever the URA does not need more torque than allowed: RA-MPC and not URA is 0 everywhere, URA and not RA-MPC 0 to 13.
- The URA exceeds the bounds in 126 to 144 of 150 runs in the forward tasks, but only in 20 (red_rev) and 48 (green_rev) runs where the common bounds leave slack (period-limited tasks whose own bounds would be (9.5, 3.5) and (17, 4.5)). The bound-violation argument against the URA is therefore strongest in the forward tasks, the paper's red task among them.
- On the common successes RA-MPC has the lower median cost in all eight tasks and a similar median time (0.65 to 0.94 s against 0.63 to 1.18 s for the URA).

### Data for the paper (all committed)

| what | where |
| --- | --- |
| task specs (psi, path, source, target, box) | `arm_top/tasks/<task>.json` |
| k1 designs, SOS logs | `arm_top/k1/<task>_d8_l<lambda>.txt`, `arm_top/data/design_<task>.out` |
| mu0 study, scale and bounds | `arm_top/data/scale_study_<task>.{json,out}`, `arm_top/data/scale_<task>.json`, `arm_top/data/solve_scale_<task>.out`, `arm_top/data/lam_<task>.txt` |
| comparison runs (initial set, all closed loops, URA, vanilla, RA-MPC logs) | `arm_top/data/<task>/main.pkl`, `run_main.out` |
| Lemma 2, timing, per-task summary | `arm_top/data/<task>/{lemma2,timing}.pkl`, `run_lemma2.out`, `run_timing.out`, `summary.txt` |
| red only: horizon calibration, period sweep | `arm_top/data/red_fwd_sq/{calib_horizon,sweep}.pkl`, `run_sweep.out` |
| all-task tables | `arm_top/data/all_tasks_summary.{json,md}` |
| figures per task (paper panels pdf/png, verification png) | `arm_top/figs/<task>/` |
| batch logs | `arm_top/data/prep_tasks.log`, `arm_top/data/all_tasks.log` |

`<task>` is one of `red_fwd_sq, blue_fwd_sq, green_fwd_sq, orange_fwd_sq, red_rev_sq, blue_rev_sq, green_rev_sq, orange_rev_sq`.

### Reproduction (environment `rab_mpc`, about 1.5 h on 10 cores)

```
./run_prep_tasks.sh blue_fwd green_fwd orange_fwd red_rev blue_rev green_rev orange_rev     # task specs, SOS designs (MATLAB), mu0 study, scale with U_FIX
./run_all_tasks.sh blue_fwd green_fwd orange_fwd red_rev blue_rev green_rev orange_rev      # comparison, Lemma 2, timing, summary, figures
python summarize_all_tasks.py                                                               # data/all_tasks_summary.{json,md}
```
