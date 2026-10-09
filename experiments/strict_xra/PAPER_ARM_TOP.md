# Top-view arm, kitting workcell: first task (red, forward)

Setting of 2026-10-09. The arm example moves to the horizontal plane ($g = 0$) and
to a kitting workcell: two vertical racks, four coloured objects on the left, four
slots on the right, a crossing assignment, a forward and a reverse phase. This
document covers the task with the largest span, the red object from the lowest
left compartment to the third right compartment, as agreed before running the
other seven. Scripts in `arm_top/` (Section 9); the option-A results of
`PAPER_ARM.md` are untouched. The paper and the notebooks are not yet changed.

Contents: 1 setting, 2 three-method comparison, 3 Lemma 2, 4 Theorem 2 and
Corollary 1, 5 period sweep and timing, 6 figures, 7 observations, 8 places to
update and open points, 9 reproduction.

## 1. Setting

### 1.1 Workcell and task

| item | value |
| --- | --- |
| arm | $m_i = 1$ kg, $l_i = 4$ m, $l_{c,i} = 2$ m, $I_i = 0.02$ kg m$^2$, horizontal plane, $g = 0$; output $y = h(q)$, branch $q_2 > 0$ |
| racks | columns $y_1 = \mp 3.2$, compartment heights $1.8, 2.8, 3.8, 4.8$ (pitch 1.0), compartments 1.5 deep |
| assignment | object $k$ to slot $\sigma(k)$, $\sigma = (3, 1, 4, 2)$: red $1.8 \to 3.8$, blue $2.8 \to 1.8$, green $3.8 \to 4.8$, orange $4.8 \to 2.8$; three crossings |
| fixture | box $\lvert y_1\rvert \le 0.6$, $1.1 \le y_2 \le 2.28$ between the racks |
| this task | red, forward: source $(-3.2, 1.8)$, slot $(3.2, 3.8)$, path length 7.0, rise 2.0 (largest span, tied with orange) |

### 1.2 Safe set and target

Reference path $y_2 = c(y_1)$ with $x = y_1/3.6$,
$c = m + d\,(3x - x^3)/2 + b\,(1 - x^2)^2$, bridge $b = 1$, $(m, d)$ such that
$c(\mp 3.2)$ equals the two compartment heights; the shape is flattest 0.4
inside the compartment, so the corridor stays level there.

$$\psi(y) = \frac{w(y_1)\,[h^2 - (y_2 - c(y_1))^2] - \varepsilon\,(y_1/3.2)^{10}}{\psi_{\max}},\qquad w = 1 + 0.02\,(y_1 + 4.5)^2,$$

$h = 0.38$, $\varepsilon = \kappa (X_0 + a) h^2 X_0 / 5$ (so that $\partial\psi/\partial y_1 = 0$ at the slot),
$\psi_{\max} = 0.3014$, degree 10. Target $\mathcal{X}_T = \{\lvert y - (3.2, 3.8)\rvert < 0.35\}$.

Properties (proved in one line and checked on the grid): $\partial\psi/\partial y_2 = 0$ only on the path, the remaining condition depends on $y_1$ alone,
so $\psi$ has exactly one critical point, its maximum, exactly at the slot. $\psi(\text{source}) = 0.45\,\psi_{\max}$.
$\lvert y\rvert \in [2.87, 5.97]$ on $\mathcal{X}_S$ (no singularity), corridor inside its compartments, 0.67 clear of the fixture,
0.59 clear of the other objects and slots.

### 1.3 Certificate, scale and bounds

| item | value |
| --- | --- |
| $k_1$ | (12) at the normalised scale, outputs $y/4$, degree 8, $\lvert k_{1,i}\rvert \le 1$ on $\{\psi \ge 0\}$, $\lambda$ fixed; $\lambda = 0.05$ is the largest of $\{0.02, 0.05, 0.1, 0.2\}$ with $\delta_s/\lambda = 0$ (0.1 gives 0.20, 0.2 gives 0.59) |
| base design check | ECGBF margin on the grid $\ge +1.0\times 10^{-3}$, $\max\lvert k_{1,i}\rvert = 0.99$ |
| $\mu_0$ | 0.1 (Section 7.3) |
| exact time scaling | with $g = 0$: $\tau_s(D_s x) = s^2\tau_1(x)$ and $T_\varepsilon(s) = T_\varepsilon(1)/s$; checked: $T_\varepsilon(1) = 33.5$ ms, $T_\varepsilon(s) = 12.50$ ms $= T_\varepsilon(1)/s$ |
| scale | $s = 0.8\,T_\varepsilon(1)/T = 2.68$ at $T = 10$ ms; $\lambda_s = 0.134$, $\mu_s = 0.719$ |
| torque bounds | $U = 1.25\,s^2 \max\lvert\tau_1\rvert$ on $\mathcal{X}_{RA}^\delta$: $\lvert\tau_1\rvert \le 35.5$, $\lvert\tau_2\rvert \le 8.5$ N m; $s_U = 3.00 > s$ |
| terminal controller on $\mathcal{X}_{RA}^\delta$ | $\max\lvert\tau\rvert = (28.3, 6.6)$ N m, inside the bounds |
| second-order bound | $\eta = 7.4\times 10^{-2}$, $M_2 = 43.4$, $T^\star = 2\eta/M_2 = 3.4$ ms (reported, not used, as for the Dubins car) |
| MPC | $T = 10$ ms, $N = 60$ (0.6 s), 5 RK4 sub-steps, state constraint at the sub-step nodes with back-off 0.005, $\delta = 0.5\,M_\Psi$, $\omega = 1.2\sup\ell/\varepsilon = 8.1\times 10^{5}$, $Q_y = 5$, $R_u = 10^{-3}$ |
| URA | the same $k_1$ at $s = 5.36$ (twice the terminal scale, four times the torque) |
| initial set | 150 states, positions uniform in $\mathcal{X}_S \setminus \mathcal{X}_T$, output velocities uniform in a disc of radius $0.3 \cdot 5.36$ around the URA's field $5.36\,k_1(y)$ (Dubins convention); 75 % of them in the URA's $\mathcal{X}_{RA}$ |

Horizon calibration on 30 states (`calib_horizon.out`): $N = 40$ vanilla 11 / RA-MPC 21 (3 states vanilla only),
$N = 60$ 16 / 28 (none vanilla only), $N = 100$ 29 / 29. $N = 60$ is the most informative.

## 2. Three-method comparison (150 states)

| | reach $\mathcal{X}_T$ | remark |
| --- | --- | --- |
| URA ($s = 5.36$) | 146 | 4 unsafe; exceeds the bounds in 142 of the successes, peaks $(118, 27)$ N m, above the bounds for a median 46 % of the run |
| vanilla MPC | 66 | feasible at $k = 0$ from 95; **29 runs lose feasibility later** (median after 31 steps) |
| RA-MPC | 139 | 11 infeasible at $k = 0$, **none afterwards** |
| terminal controller alone, $s = 2.68$ | 80 continuous, 79 sampled | |

Pairwise: RA-MPC and not vanilla 73, **vanilla and not RA-MPC 0**, RA-MPC and not URA 2, URA and not RA-MPC 9, all three 66, none 2.
Nine of the eleven RA-MPC failures start within 1.8 of the slot with the URA's speed and cannot stop inside the corridor under the bounds; two start at the far end.

| | states | time to $\mathcal{X}_T$, median (max) | cost, median |
| --- | --- | --- | --- |
| URA vs RA-MPC | 137 | 0.85 (3.87) vs 0.82 (1.56) s | 26.8 vs 19.7 |
| all three | 66 | URA 0.37, RA-MPC 0.38, vanilla 0.37 s | 3.0, 2.4, 2.2 |

RA-MPC uses the bounds actively: $\tau_1$ within 2 % of its bound 46 % of the time, $\tau_2$ 86 %.
Smallest $\psi$ along the continuous closed loops before entering $\mathcal{X}_T$: RA-MPC 0.0055, vanilla 0.0050.

## 3. Verification of Lemma 2 (held torque, terminal set)

| $T$ [ms] | arcs | smallest one-step increase of $\Psi$ | increasing | end in terminal set | second-order bound |
| --- | --- | --- | --- | --- | --- |
| 1 | 6611 | $+7.2\times 10^{-5}$ | 100 % | 100 % | $+5.3\times 10^{-5}$ |
| 5 | 6596 | $+2.3\times 10^{-4}$ | 100 % | 100 % | $-1.8\times 10^{-4}$ |
| **10** | 6578 | $+2.5\times 10^{-4}$ | 100 % | 100 % | $-1.5\times 10^{-3}$ |
| 12.5 | 6571 | $+1.2\times 10^{-4}$ | 100 % | 100 % | $-2.5\times 10^{-3}$ |
| 15 | 6559 | $-6.3\times 10^{-4}$ | 99.74 % | 99.77 % | |
| 20 | 6549 | $-2.9\times 10^{-3}$ | 99.14 % | 99.22 % | |
| 50 | 6429 | $-3.9\times 10^{-2}$ | 94.51 % | 95.86 % | |

Sampled-data closed loop (400 states each): from $\mathcal{X}_{RA}^\delta$ all reach $\mathcal{X}_T$ safely at every period up to 30 ms, with $\Psi$ increasing at every sample up to 12.5 ms (396/400 at 15 ms); from $\mathcal{X}_{RA}$ all reach $\mathcal{X}_T$ at every period, median 0.79 s.
The guarantee of the lemma at $T = 10$ ms rests on the sampled check (threshold 12.5 ms), the second-order bound (3.4 ms) is conservative by a factor 3.7.

## 4. Verification of Theorem 2 and Corollary 1

139 runs feasible at $k = 0$, all reach $\mathcal{X}_T$; 10 405 solves, no infeasible solve after $k = 0$; the shifted candidate was kept in 5254 steps;
$u^c$ at the predicted terminal states within the bounds in every step.
Cost decrease over 10 266 steps: smallest 0.78, $J^*_k - J^*_{k+1} \ge \ell(x_k, u_k)$ in every step (smallest ratio 0.999999), never below $\ell_{\min} = 0.6125$.
Horizon nonincreasing in all runs. Corollary 1 bound holds in 139/139 runs (entry median 0.82 s, bound median 1816 s: valid, very loose because $\omega$ is large).

## 5. Period sweep and timing

40 states, $NT = 0.6$ s, terminal design ($\delta$, $\omega$) fixed:

| $T$ [ms] | $N$ | vanilla feasible (reached) | $u^c$ held | RA-MPC | later infeasible | decrease violations | min $\Delta\Psi$ on terminal set |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 5 | 120 | 26 (19) | 25 | 38 | 0 | 0/5249 | $+2.7\times 10^{-4}$ |
| 7.5 | 80 | 26 (19) | 25 | 38 | 0 | 0/3494 | $+3.0\times 10^{-4}$ |
| 10 | 60 | 26 (19) | 25 | 38 | 0 | 0/2616 | $+2.8\times 10^{-4}$ |
| 12.5 | 48 | 26 (19) | 25 | 38 | 0 | 0/2087 | $+1.8\times 10^{-4}$ |
| 15 | 40 | 26 (19) | 25 | 38 | 0 | 0/1738 | $-4.8\times 10^{-4}$ |
| 20 | 30 | 26 (19) | 25 | 38 | 0 | 0/1298 | $-2.7\times 10^{-3}$ |

The success rate is flat; beyond 12.5 ms the one-step increase fails on part of the terminal set, so the guarantee lapses although the closed loop still succeeds.

Timing (one run at a time, Apple M5): RA-MPC 185 ms per solve (first 453 ms), vanilla 129 ms, $u^c$ 0.24 ms. The MPC is about 18 times slower than the 10 ms period here (Section 8).

## 6. Figures (`arm_top/figs/red_fwd/`)

| file | content |
| --- | --- |
| `fig_manip_vanilla_far`, `fig_manip_uncon_far`, `fig_manip_rampc_far` | paper panels in the workcell: corridor, slot, successes and failures over the initial positions, trajectory from the farthest common start (URA 1.34 s, peaks (56, 24) N m; RA-MPC 1.50 s) |
| `fig_manip_compare_far` | $\rho$ against $T$, torques of URA and RA-MPC against the bounds |
| `fig_manip_cost_far` | $J^*_k$, decrease against $\ell$, horizon |
| same names without `_far` | shared start where all three succeed (state 13) |
| `verify_outcomes`, `verify_lemma2`, `verify_theorem2`, `verify_sweep` | verification panels of Sections 2 to 5 |
| `../../arm/figs/arm_racks.png` | the workcell with all four tasks in both phases |

## 7. Observations

1. The comparison now has the intended shape: vanilla-feasible is a subset of RA-MPC-feasible, the URA and RA-MPC orders are mixed, and the URA succeeds only by exceeding the bounds.
2. Vanilla MPC loses feasibility in 29 of its 95 initially feasible runs (no terminal ingredient); RA-MPC never does. This is a direct illustration of recursive feasibility and worth a sentence in the paper.
3. The choice $\mu_0 = 0.1$ makes rest states mostly uncertified (7 % in $\mathcal{X}_{RA}$); $\mu_0 = 1$ would certify 75 % of rest states but makes $T_\varepsilon(1)$ 4.5 ms and the runs about eight times longer. With $\mu_0 = 0.1$ the initial set follows the Dubins convention (velocities around the URA's field). Rest starts remain possible for the demonstration cycle.

## 8. Places to update in the paper, and open points

- Example 2 statement: horizontal plane ($g = 0$, no gravity term), workcell and task description, $\psi$ of Section 1.2, target disc, bounds 35.5 / 8.5 N m, $T = 10$ ms, $N = 60$.
- Table I: re-measure on the paper's machine; here RA-MPC 185 ms per step against a 10 ms period, so the example is not real time on this hardware (as for the option-A arm).
- Open: the other seven certificates (same pipeline, `TASK=<colour>_<fwd|rev>`), and one RA-MPC run of the full forward and reverse cycle for the illustration and the later simulation; transit moves between tasks are outside the planar model.
- Open, carried over: second-order Lemma 2 in the paper (D7), notebooks not updated, `controllers/*.py` untouched.

## 9. Reproduction (`experiments/strict_xra/arm_top/`, environment `rab_mpc`)

```
python make_task.py red fwd                       # tasks/red_fwd.json from ../arm/screen_racks.py
matlab -batch "TASK='red_fwd'; DEG=8; DESIGN_LAMS=[0.02 0.05 0.1 0.2]; run('run_design_top.m')"
python scale_study_top.py                         # mu0 trade-off
python solve_scale_top.py                         # s, bounds, URA scale -> data/scale_red_fwd.json
python calib_horizon.py 40,60,100 30
python run_main_top.py && python run_lemma2_top.py && python run_sweep_top.py 40 && python run_timing_top.py
python summarize_top.py > data/red_fwd/summary.txt && python make_figs_top.py
```
