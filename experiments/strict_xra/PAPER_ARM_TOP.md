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
- Open: the other seven certificates (same pipeline, `LAYOUT=square`, `TASK=<colour>_<fwd|rev>_sq`), and one RA-MPC run of the full forward and reverse cycle for the illustration and the later simulation; transit moves between tasks are outside the planar model.
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
