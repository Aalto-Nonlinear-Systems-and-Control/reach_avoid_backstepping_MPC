# TAC revision: experiment update plan (draft for confirmation)

Branch `tac-revision`, created from `lcss-revision` (the branch that produced
the figures currently in the TAC manuscript). This file lists what must
change in the code so that the experiments match the revised Section IV of
the paper (`TAC/method.tex`: problem (21), terminal-cost condition (22),
Algorithm 1, Corollary 1, Remark 4). Nothing below has been implemented yet.
Tick the boxes as items are confirmed.

Hard constraint carried over from the LCSS revision: do not overwrite
`controllers/*.py`. They are the certified controllers behind the reported
numbers. Any MATLAB re-run must use throwaway export names.

Note for editing: keep `|` out of inline math (write `\mid`) and avoid `\|`,
otherwise the Markdown table formatter of the IDE breaks the tables.

## 1. Paper changes and their code implications

### 1.1 Terminal set

- Paper: terminal set is the shrunk set $\mathcal{X}_{RA}^{\delta} = \{\Psi \ge \delta\}$
  with $\delta > \beta(T)/\lambda$, so that the sampling condition of Lemma 2
  holds at every terminal state.
- Code now: `cert_vec(x_N) >= 0`, the zero-superlevel set.
- Change: `cert_vec(x_N) >= delta`, with `delta` computed from the constants
  in Section 2.

### 1.2 Cost function

- Paper: $J = \sum_{i<N} \ell(\hat{x}_i, u_i) + V_f(\hat{x}_N)$ with
  $\ell \ge \underline{\ell} > 0$ outside $\mathcal{X}_T$ and $V_f$
  satisfying the one-step decrease (22) along ZOH $u^c$.
- Code now: stage cost $\lVert h(x) - y^* \rVert_Q^2 + \lVert u \rVert_R^2$,
  terminal cost $\lVert h(x_N) - y^* \rVert_{Q_f}^2$.
- Change: keep the stage cost (its lower bound outside $\mathcal{X}_T$ is
  automatic, see Section 2). Replace the quadratic terminal term by
  $V_f(x) = \omega\,(M_\Psi - \Psi(x))$ as in Remark 4. Adding $V_f$ on top
  of the $Q_f$ term is not allowed: the tracking term need not decrease
  along $u^c$, which would break (22).

### 1.3 Sampling period and the `dt` sweep

- Paper: fixed period $T$ with $T < T^*(\delta)$, where $T^*(\delta)$ is the
  largest $T$ with $\beta(T)/\lambda < \delta$. Algorithm 1 takes $\delta$
  as input and picks $T \in (0, \beta^{-1}(\lambda\delta))$ offline; the
  experiments do the converse (fix `dt`, then choose $\delta > \beta(T)/\lambda$).
  Both satisfy the same condition, so no code change follows from this.
- Code now: fixed `dt_mpc` for the base runs; the sweep varies `dt` with the
  same controller design.
- Change for the base runs: none, but $\delta$ and $\omega$ are computed at
  the base period $T_{\text{base}}$ (0.05 s Dubins, 0.01 s manipulator).
- Change for the sweep: keep the sweep, but with $\delta$ and $\omega$
  fixed at their base values, so that only the sampling period changes.
  Recomputing $\delta(T)$, $\omega(T)$ per point would shrink the terminal
  set with $T$ and force $\rho = 0$ by construction once
  $\beta(T)/\lambda \ge M_\Psi$; that curve would measure the
  conservativeness of $\beta(T)$, not the closed loop. Mark the certified
  range $T < T^*(\delta)$ on the $\rho$ versus $T$ plot. Inside it the
  theory predicts success from every initially feasible state (nominal
  model), so any failure there signals an implementation or integration
  error; outside it the curve is empirical, as in the current paper text.

### 1.4 Shifted candidate (Algorithm 1, line 7)

- Paper: $\mathbf{U}_{k+1} = \{u^*_{1 \mid k}, \dots, u^*_{N-1 \mid k}, u^c(\hat{x}^*(NT \mid k))\}$.
- Code now: the warm start repeats the last input (`U_warm[:, -1:]`).
- Change: build the candidate with $u^c$ evaluated at the predicted terminal
  state (numpy `lambdify` of `u_opt`). Verify that this input respects the
  bounds. At $k = 0$ there is no candidate ($\mathbf{U}_0 = \emptyset$ in
  Algorithm 1): the first solve uses the solver's default initial guess, as
  the notebooks already do.

### 1.5 Cost comparison rule (suboptimal MPC)

- Paper: the corollary only needs the applied sequence to cost no more than
  the candidate.
- Code now: the solver output is used as is.
- Change: warm start at the candidate and compare the returned cost with the
  candidate cost; keep the cheaper one.

### 1.6 Adaptive horizon (Algorithm 1, lines 9 and 10)

- Paper: after building $\mathbf{U}_{k+1}$, roll it out over $[0, NT]$. If the
  rollout enters $\mathcal{X}_T$ at any time, set $N \gets \max(N - 1, 1)$.
  There is no fallback and no open-loop execution: the next step solves (21)
  with the new horizon, warm-started at $\mathbf{U}_{k+1}$ (after a shrink the
  candidate has one entry more than the horizon; only the first $N$ entries
  are used). The shortened candidate $\{u^*_{1 \mid k}, \dots, u^*_{N-1 \mid k}\}$
  ends at $\hat{x}^*(NT \mid k) \in \mathcal{X}_{RA}^{\delta}$, so it is feasible
  and costs $J^*_k - \ell(x_k, u^*_{0 \mid k})$; no $u^c$ step is ever appended
  from a state at or inside the target. The horizon never grows back.
- Code now: `N_hor` is a fixed global; an infeasible solve marks the run as
  failed.
- Change: make the horizon a per-step variable carried through the closed
  loop. Dubins: `solve_mpc` builds a fresh `Opti` per call, so pass `N` as an
  argument. Manipulator: the pre-compiled solver is built for one `N`; cache
  one solver per horizon length (at most `N_hor` of them, built on demand).
  Check the candidate rollout against $\mathcal{X}_T$ at every RK4 node and,
  for the last ($u^c$) segment, at integrator sub-steps as well, since the
  paper's condition is on the continuous trajectory. Under the nominal model
  an infeasible solve can no longer occur; treat it as a violation and log it.

### 1.7 Time bound (Corollary 1)

- Paper: entry into $\mathcal{X}_T$ no later than
  $t = (\lfloor J_0^*/\underline{\ell} \rfloor + 1)\,T$; the proof uses
  $J^*_{k+1} \le J^*_k - \underline{\ell}$, which holds both for the full
  candidate (via (22)) and for the shortened one (same terminal state).
- Code now: $J^*_k$ is not logged.
- Change: log $J^*_k$ per closed-loop step and check both the monotone
  decrease and the bound. Validation only unless we decide to plot it.

## 2. Constants needed per system

None of these exist on the Python side yet. Notation follows the paper.
"Over $\mathcal{X}$" means over the compact sampling box used in MATLAB
(`bound_min` and `bound_max` in `example_*.m`) unless we decide to use
$\{\Psi \ge 0\}$ instead (decision D3).

| Symbol | Meaning | How to get it |
| --- | --- | --- |
| $\lambda$ | ECGBF rate | Numeric `k1_lambda` inside `solvesop_bounded_control.m`, not exported. Option (a): re-run `example_*.m` with throwaway export names and print it. Option (b): estimate as the minimum of $\dot{\Psi}/\Psi$ under $u^c$ over samples with $\Psi \ge \epsilon$. Prefer (a). |
| $M_{g_j}$ | max of $\lvert L_{g_j}\Psi \rvert$ over $\mathcal{X}$ | sympy gradient of `certificate_opt` times column $g_j$, sampled |
| $\ell_{u^c_j}$ | Lipschitz constant of $u^c_j$ | max of $\lVert \nabla u^c_j \rVert_2$ over $\mathcal{X}$, sampled |
| $\ell_f$, $\ell_{g_j}$ | Lipschitz constants of $f$ and $g_j$ | max induced 2-norm of the Jacobians, sampled. Dubins: $g$ constant, so $\ell_{g_j} = 0$. Manipulator: $g = M(q)^{-1}$ columns, nonzero. |
| $M_{u^c_j}$ | max of $\lvert u^c_j \rvert$ over $\mathcal{X}$ | sampled; should not exceed the input bounds on $\{\Psi \ge 0\}$ |
| $\ell_{sys}$ | $\ell_f + \sum_j \ell_{g_j} M_{u^c_j}$ | formula |
| $M_v$ | max of $\lVert f + \sum_j g_j u^c_j \rVert$ over $\mathcal{X}$ | sampled |
| $\gamma(T)$ | $(\sum_j M_{g_j}\ell_{u^c_j})\,(e^{\ell_{sys}T} - 1)/\ell_{sys}$ | formula |
| $\beta(T)$ | $\gamma(T)\,M_v$ | formula |
| $\delta$ | terminal margin, $\delta > \beta(T)/\lambda$ | choice, see D2 |
| $M_\Psi$ | upper bound of $\Psi$ on the closure of $\mathcal{X}_{RA}$ | sampled max of `certificate_opt` on $\{\Psi \ge 0\}$, slightly inflated |
| $\sup \ell$ | max stage cost on $\mathcal{X} \times \mathcal{U}$ | $Q$ times max of $\lVert h(x) - y^* \rVert^2$ plus $R\,u_{\max}^2$, sampled |
| $\underline{\ell}$ | stage-cost lower bound outside $\mathcal{X}_T$ | automatic from the tracking term. Dubins target $y_2^2 + 4(y_1 + 1.7)^2 \le 0.4$ gives $\lVert h(x) - y^* \rVert^2 > 0.1$ outside, so $\underline{\ell} = 0.5$ with $Q = 5I$. Manipulator ellipse gives $\lVert h(x) - y^* \rVert^2 > 0.32$, so $\underline{\ell} = 1.6$. Verify numerically. |
| $\omega$ | terminal-cost weight | $\omega \ge \sup\ell / [(e^{\lambda T} - 1)(\delta - \beta(T)/\lambda)]$ from Remark 4; use the minimal admissible value |

Caveat on the exported certificate. `certificate_opt` is already shifted by
$-\delta_{SOS}/\lambda$ in MATLAB (`reach_avoid_controller.m`), so it is the
function that satisfies $\dot{\Psi} \ge \lambda \Psi$ exactly and it is the
paper's $\Psi$. The SOS slack $\delta_{SOS}$ is unrelated to the terminal
margin $\delta$ above (the same symbol clash exists in the paper).

Risk. $\beta(T)$ is built from global Lipschitz constants and may be
conservative. If $\beta(T)/\lambda$ comes out close to or above $M_\Psi$ at
$T = 0.05$ s (Dubins) or $T = 0.01$ s (manipulator), the terminal set is
empty and the theory demands a smaller $T$. Compute the numbers first; this
is the gate for everything else.

Proposed home: a new module `python/ra_mpc_constants.py` with one function
per system that returns the constants together with `delta(T)` and
`omega(T)`, plus a notebook cell that prints them.

## 3. File-by-file changes

### `notebooks/example_dubins_car_reach_avoid_mpc.ipynb`

- [x] cell 4, `solve_mpc`: terminal constraint `cert_vec(X[:, N]) >= delta`.
- [x] cell 4, cost: drop the `Qf_y` term; add `omega * (M_psi - cert_vec(X[:, N]))`.
- [x] cell 4, warm start: shifted candidate whose last input is `u_c(x_N_pred)`
      (new numpy callable from `u_opt`); roll it out with the same RK4 map
      `F_disc` (the current rollout uses explicit Euler).
- [x] cell 4, return value: also return the solver cost; add a
      `candidate_cost(x0, U)` helper and apply the "keep the cheaper" rule.
- [x] cells 5 and 7, `_feasible_receding` and `simulate_ra_mpc_traj`: carry
      `N` through the loop; after each solve build the candidate, roll it
      out, and shrink `N` if the rollout touches $\mathcal{X}_T$ (1.6). Any
      infeasible solve is now a violation; log it separately from "failed".
- [x] log $J^*_k$ and check $J^*_{k+1} \le J^*_k - \underline{\ell}$.
- [x] header comment and printouts (they still say `V(x_N) >= 0`).

### `notebooks/example_manipulator_reach_avoid_mpc.ipynb`

- [ ] cell 3: the same four changes in `solve_mpc` and in the shared
      pre-compiled solver (`solve_mpc_fast`, `make_fast_solver`); add a
      per-horizon solver cache for the adaptive `N`.
- [ ] cell 4: the same closed-loop changes.

### `notebooks/make_paper_figs_dubins.ipynb` and `notebooks/make_paper_figs_manipulator.ipynb`

- [~] `feas_rampc` and `make_fast_solver(..., "cert")`: the same MPC changes.
- [ ] `dt` sweep: keep `delta` and `omega` at their base values for every
      entry of `DT_SWEEP` (see 1.3); `_sweep_sig()` must include the constants
      or a version tag so that the caches in `figures/*.npz` are invalidated.
- [ ] compute $T^*(\delta)$ and draw it on the sweep panel (vertical line or
      shaded certified range); check that $\rho = 100\%$ for the sweep points
      inside it.
- [ ] optional (D4): a small table with $\delta(T)$, $\omega(T)$ and whether
      $\{\Psi \ge \delta(T)\}$ is non-empty for each `dt`, to show how the
      certified design changes with $T$.

### `notebooks/noise_robustness_dubins.ipynb`

- [ ] builds its own RA-MPC solver: the same changes.

### Baseline notebooks

- `example_*_vanilla_mpc.ipynb` and `example_*_unconstrained_reach_avoid.ipynb`:
  no change.

### `python/`

- [x] add `python/ra_mpc_constants.py` (Section 2).
- [ ] optional (D5): factor the RA-MPC builder into `python/ra_mpc.py` so
      that the four notebooks share one implementation.

### `matlab/`

- [ ] `export_to_python.m`: also write `lambda_val` and `mu_val` into the
      exported file. Run only with throwaway export names, then copy the
      scalar by hand into `ra_mpc_constants.py` or a small JSON file next to
      the controller file. Do not regenerate `controllers/*.py`.

### `README.md`

- [ ] "Controller Settings Per Notebook": terminal constraint column becomes
      $V(x_N) \ge \delta$; add $\delta$ and $\omega$ columns; drop $Q_{f,y}$
      for the RA-MPC rows.
- [ ] "MPC Formulation": new cost and terminal set; mention the candidate
      warm start, the cost comparison rule, and the fallback.

### `data/` and `figures/`

- [ ] regenerate `traj_controls_*.npz`, the `*_cache.npz` files, and the
      PNG/PDF figures; copy the new figures to `TAC/figs/` in the paper repo.

## 4. Validation checklist (per system, before regenerating figures)

- [ ] constants printed and sanity-checked; $\beta(T)/\lambda < M_\Psi$ with
      room for $\delta$.
- [ ] $\{\Psi \ge \delta\}$ non-empty and still reachable within $N$ steps from
      a comparable share of initial states (compare feasibility counts before
      and after).
- [ ] candidate input $u^c(\hat{x}^*(NT \mid k))$ within the bounds at every step.
- [ ] $J^*_k$ decreases by at least $\underline{\ell}$ at every solved step.
- [ ] reaching time no later than
      $(\lfloor J_0^*/\underline{\ell} \rfloor + 1)\,T$ on every successful run.
- [ ] no infeasible solve at all; `N` nonincreasing and never below 1.
- [ ] success rates and closed-loop plots regenerated; numbers in the paper
      (Fig. 4 and Fig. 6 counts, sweep curves, computation-time table) updated.

## 5. Open decisions

- D1, resolved: the paper now uses the adaptive horizon of 1.6 (no
  fallback, no open-loop execution). Code follows the paper.
- D2, choice of $\delta$. Candidates: $\delta = 2\beta(T)/\lambda$, or
  $\delta = \beta(T)/\lambda + \theta M_\Psi$ with a small $\theta$. A larger
  $\delta$ shrinks the terminal set but enlarges the guaranteed increment and
  reduces $\omega$.
- D3, domain for the constants. Sampling box $\mathcal{X}$ from MATLAB versus
  the tighter $\{\Psi \ge 0\}$. The paper states the constants on $\mathcal{X}$.
- D4, sweep presentation. The main sweep uses the base design (1.3). Decide
  whether to add the per-`dt` redesign table ($\delta(T)$, $\omega(T)$,
  terminal set non-empty or not), and how the paper words the gap between
  the certified $T^*(\delta)$ and the empirically successful periods.
- D5, refactor. Shared `python/ra_mpc.py` versus editing the four notebooks in
  place.
- D6, magnitude of $\omega$. The theoretical minimum may dominate the tracking
  cost. Report the value; if the closed-loop behaviour degrades, discuss
  whether to present the tracking cost as illustrative in the paper.

## 6. Results: Dubins car (branch `tac-revision`, for confirmation)

Superseded for the paper: these results use the committed controller, whose
certificate is not rigorous. The final Dubins setting, data and figures are in
`experiments/strict_xra/PAPER_DUBINS.md`.

Implementation: `python/ra_mpc_dubins.py` (Algorithm 1: terminal set
$\{\Psi \ge \delta\}$, $V_f = \omega(M_\Psi - \Psi)$, shifted candidate with
cheaper-of rule, adaptive horizon), `python/ra_mpc_constants.py` (Lemma 2
constants), `scripts/eval_dubins_tac.py` (batch evaluation). Both Dubins
notebooks now call the module; `make_paper_figs_dubins.ipynb` feasibility and
`dt` sweep use Algorithm 1 for RA-MPC (vanilla and unconstrained unchanged).
The manipulator item for `make_fast_solver` is still open.

Settings: $T = 0.05$ s, $N = 25$, $\delta = 0.1$, $\omega = 1000$,
$M_\Psi = 11.63$, $\ell_{\min} = 0.5$, stage cost unchanged, safety
constraint tightened to $\psi_{\text{safe}} \ge 10^{-3}$ (see 6.3).

### 6.1 Theory constants are vacuous for this controller

- MATLAB `solve_vanilla_k1_controller` gives $\lambda = 1.09 \times 10^{-8}$
  (independent numerical estimate of $\min \dot\Psi/\Psi$ agrees, about $10^{-6}$ or below).
- On $\{\Psi \ge 0\}$ within the synthesis box ($\theta \in [2\pi/3, 4\pi/3]$,
  $v \in [0.1, 1]$): $M_{u^c} = (96, 11)$, so $u^c$ exceeds the input bounds on
  34 % ($\omega$) and 6 % ($a$) of that set; $\ell_{u^c} = (1246, 85)$, $M_v = 106$.
- Hence $\beta(0.05)/\lambda \approx 4.7 \times 10^{11} \gg M_\Psi$ and the certified
  period is $T^*(\delta) \approx 10^{-14}$ s. Remark 4's $\omega$ is infinite.
- Under one ZOH step of $u^c$ ($T = 0.05$ s) from $\{\Psi \ge \delta\} \setminus \mathcal{X}_T$,
  $\Psi$ decreases in 44 % of sampled states.
- Consequence: $\delta$ and $\omega$ are practical choices; (22) is not
  guaranteed a priori and is checked a posteriori instead.

### 6.2 Closed loop, nominal plant (MPC model = plant), 300 initial states, 10 s

| $\omega$ | reached | infeasible at $k=0$ | infeasible later | timeout | steps violating $J^*_{k+1} \le J^*_k - \ell_{\min}$ (reached runs) | $t_{\text{entry}} \le t_{\text{bound}}$ |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 285 | 12 | 0 | 3 | 1 / 5964 | 285 / 285 |
| 100 | 285 | 10 | 0 | 5 | 0 / 5935 | 285 / 285 |
| 1000 | 284 | 11 | 0 | 5 | 0 / 5900 | 284 / 284 |

- Recursive feasibility holds numerically: no infeasible solve after $k = 0$.
- The candidate $u^c$ step was within the input bounds at every step, and the
  horizon shrank in every successful run (the plan reaches $\mathcal{X}_T$ well
  within $N$ from about 0.5 s on).
- Time bound: all reached runs satisfy it, but it is very loose:
  $t_{\text{entry}}$ 0.15 to 1.85 s versus $t_{\text{bound}}$ 102 to 359 s
  (ratio below 0.01). The reason is that $J^*$ floors near
  $\omega(M_\Psi - \max\Psi) \approx 1000$, not near 0 (figure
  `figures/tac_dubins_cost.png`). Using $J^*_0 - \inf J^*$ in place of $J^*_0$
  would tighten the corollary; worth considering for the paper.
- The 5 timeouts all stall at $(1.41, 0)$ with $v = 0$, the point of the ring
  opposite the target. There the decoupling matrix is singular (Assumption 1
  fails), $\Psi$ does not increase under $u^c$, and (22) is violated at almost
  every step (J decreases by about 0.4 per step, less than $\ell_{\min}$).
- Solve time about 40 ms per step (fresh `Opti` per call, larger than the 22 ms
  of the pre-built solver below).

### 6.3 Closed loop, continuous ZOH plant (as in the paper figures), 20 s

- Without the safety back-off, 40 runs ended "unsafe" with violations of order
  $10^{-5}$ or smaller while travelling at 3 to 4 m/s along the boundary: the
  RK4 prediction error at $T = 0.05$ s pushes the state just outside. With the
  back-off $\psi_{\text{safe}} \ge 10^{-3}$: no safety violation.
- Result with back-off: RA-MPC reaches the target from 258 / 300 states
  (paper currently reports 255 / 300 with the LCSS formulation); vanilla MPC
  245 / 300, unconstrained 197 / 300 (unchanged).
- 26 runs become infeasible mid-run, all on the tightened boundary at about
  3.3 m/s: plant-model mismatch breaks recursive feasibility (no such case with
  the nominal plant). This is the robustness issue of the paper's model
  uncertainty remark, not a defect of Algorithm 1.
- Shared trajectory of the paper figure: reached at 1.40 s, $t_{\text{bound}} = 249$ s,
  no decrease violation in 27 steps (figure `figures/fig_dubins_cost.png`).
- `dt` sweep (150 states, $T N = 1.25$ s): $\rho$ = 96.7, 88.0, 77.3, 84.0 %
  for $T$ = 0.025, 0.05, 0.125, 0.25 s. Not monotone any more (was 96.7, 85.3,
  66.7, 58.7 % before the back-off, and monotone in the LCSS version). The
  certified range $T < T^*(\delta)$ is about $10^{-14}$ s and cannot be drawn.

### 6.4 Points to decide

1. How the paper treats 6.1: the Lemma 2 constants are existence-level only for
   this controller; the experiments cannot certify $T$, $\delta$, $\omega$.
2. Whether to keep $\omega = 1000$ (results are insensitive to $\omega$ in
   10 to 1000) and $\delta = 0.1$.
3. Whether to tighten Corollary 1 with $J^*_0 - \inf J^*$.
4. Whether to state the safety back-off in the paper (it is a standard
   constraint tightening for discretisation error).
5. The non-monotone `dt` sweep and the new counts (258 / 300, figure 4 and 5 of the paper).
6. Whether to update the paper figures now (`figures/fig_dubins_*.png`
   regenerated on this branch; not yet copied to `TAC/figs/`).

## 7. Time scaling of the backstepping design (derivation, for confirmation)

Purpose: replace the input-constraint part of the synthesis by a scalar
design variable $s$ (the time scale of the terminal controller), with the
admissible range of $s$ derived from the backstepping relations, and tie the
sampling condition of Lemma 2 to the same variable. Written for the
relative-degree-two case of the examples; the extension to higher relative
degree needs level-dependent exponents and is not written here.

### 7.1 Base design

System $\dot x = f(x) + g(x)u$, output $y = h(x) \in \mathbb{R}^p$, and in
output coordinates

$$
\ddot y = b(x) + A(x)\,u ,
$$

with $A$ nonsingular on $\overline{\mathcal{X}_S \setminus \mathcal{X}_T}$
(Assumption 1). The state is split as $x = (q, \dot q)$ into configuration and
velocity coordinates, with $y = h(q)$ and $\dot y = J(q)\dot q$.

Program (12) solved with $\lambda$ **fixed** and the normalisation
$\lvert k_{1,i} \rvert \le K$ on $\{\psi \ge 0\}$ returns $(k_1, \delta_s)$ with

$$
\nabla\psi^\top k_1 - \lambda\psi + \delta_s \;\ge\; 0
\quad \text{on } \mathcal{X}_S \setminus \mathcal{X}_T
\tag{7.1}
$$

(the SOS multipliers are omitted). With a chosen $\mu > 0$ the certificate and
the controller are

$$
\Psi(x) = \psi(y) - \frac{\delta_s}{\lambda} - \frac{\lVert e \rVert^2}{2\mu},
\qquad e = \dot y - k_1(y),
\tag{7.2}
$$

$$
u^c(x) = A(x)^{-1}\big(\ddot y_{des}(x) - b(x)\big), \qquad
\ddot y_{des} = \mu\nabla\psi + J_{k_1}(y)\,\dot y + \tfrac{\lambda}{2}\,e .
\tag{7.3}
$$

Along the closed loop $\ddot y = \ddot y_{des}$, so by (7.1)

$$
\dot\Psi = \nabla\psi^\top\dot y - \frac{e^\top(\ddot y - J_{k_1}\dot y)}{\mu}
= \nabla\psi^\top k_1 - \frac{\lambda\lVert e\rVert^2}{2\mu} \;\ge\; \lambda\Psi .
$$

This is the ECGBF condition of the paper; (7.3) is what `design/cert.py`
implements (`uc`).

### 7.2 The scaled family

For $s > 0$ define the dilation $D_s(q, \dot q) = (q, s\dot q)$ and the scaled
design

$$
k_1^s = s\,k_1, \qquad \lambda^s = s\lambda, \qquad \delta_s^s = s\,\delta_s,
\qquad \mu^s = s^2\mu .
\tag{7.4}
$$

Claims, each one line to check:

1. **(7.1) is preserved.** Multiplying (7.1) by $s$ gives
   $\nabla\psi^\top k_1^s - \lambda^s\psi + \delta_s^s \ge 0$, and the SOS
   multipliers, multiplied by $s$, stay SOS. The ratio
   $\delta_s^s/\lambda^s = \delta_s/\lambda$ and the position set
   $\{\psi \ge \delta_s/\lambda\}$ do not change. The bound becomes
   $\lvert k^s_{1,i}\rvert \le sK$.
2. **Certificate.** With $\dot y = J(q)\dot q$,

   $$
   \Psi^s(q, \dot q)
   = \psi(y) - \frac{\delta_s}{\lambda} - \frac{\lVert \dot y - s\,k_1(y)\rVert^2}{2 s^2\mu}
   = \Psi\big(q, \dot q/s\big),
   $$

   that is $\Psi^s = \Psi \circ D_s^{-1}$. Hence
   $\mathcal{X}_{RA}^s = D_s\mathcal{X}_{RA}$: the same positions with the
   velocities scaled by $s$, while $M_\Psi$, $\delta$ and the share of
   certified positions are unchanged.
3. **Closed loop.** If $y(t)$ solves $\ddot y = \ddot y_{des}$ for the base
   design, then $y^s(t) = y(st)$ solves it for the scaled design, because
   $\ddot y^s(t) = s^2\,\ddot y_{des}(x(st))$ and the desired acceleration of
   the scaled design at $D_s x$ is

   $$
   s^2\mu\nabla\psi + J_{s k_1}\,(s\dot y) + \frac{s\lambda}{2}\,(s\dot y - s k_1)
   = s^2\,\ddot y_{des}(x).
   $$

   Consequently $\Psi^s(x^s(t)) = \Psi(x(st))$: the same path, $s$ times
   faster, with the same certificate values.
4. **Inputs.** At corresponding states
   $u^{c,s}(D_s x) = A(D_s x)^{-1}\big(s^2\ddot y_{des}(x) - b(D_s x)\big)$.
   Expanding $b$ and $A$ by their homogeneity in the velocities gives a finite
   sum

   $$
   u^{c,s}(D_s x) = \sum_{k} s^{k}\,u_k(x),
   \tag{7.5}
   $$

   with the $u_k$ computable from the base design. Two cases of the paper:
   - Dubins car: $b = 0$, $A(D_s x) = A(x)\,\mathrm{diag}(s, 1)$, so
     $\omega^s = s\,\omega_1$, $a^s = s^2 a_1$.
   - Manipulator: $A = J M^{-1}$ depends on $q$ only;
     $b = J M^{-1}(-C\dot q - G) + \dot J\dot q$ with $C\dot q$ and
     $\dot J\dot q$ quadratic in $\dot q$ and $G$ independent of it, so
     $\tau^s = G(q) + s^2\,\tilde\tau_1(x)$, where $\tilde\tau_1$ is the base
     torque without the gravity term.

### 7.3 Admissible scales

**Input bounds.** Theorem 2 needs $u^c \in \mathcal{U}$ on the terminal set
$\mathcal{X}_{RA}^{\delta}$ only. With (7.5) this is, for every
$x \in \mathcal{X}_{RA}^{\delta}$ of the base design, a polynomial inequality
in $s$; the admissible range is $s \le s_U$ with

$$
s_U = \sup\Big\{ s > 0 : \sum_k s^k u_k(x) \in \mathcal{U}
\ \ \forall x \in \mathcal{X}_{RA}^{\delta} \Big\}.
\tag{7.6}
$$

For box bounds and the Dubins car this is explicit,

$$
s_U = \min\Big( \frac{\omega_{\max}}{\sup_{\mathcal{X}_{RA}^{\delta}} \lvert\omega_1\rvert},\
\sqrt{\frac{a_{\max}}{\sup_{\mathcal{X}_{RA}^{\delta}} \lvert a_1\rvert}} \Big);
\tag{7.7}
$$

for the manipulator it is the smallest root over the set of the quadratic
inequalities $\lvert G_i(q) + s^2\tilde\tau_{1,i}(x)\rvert \le \tau_{\max,i}$,
which requires the static condition $\lvert G_i(q)\rvert < \tau_{\max,i}$ on the
set. If the zero-order term alone violates the bounds, no $s$ helps and the
shape or the set has to change; that is the case left for (14).

**Sampling period.** In the second-order form of Lemma 2,
$\Psi(x^s(t_{k+1})) \ge \Psi(x^s(t_k)) + \eta T - \tfrac12 M_2 T^2$ with
$\eta = \min_{\mathcal{X}_{RA}^{\delta}}\dot\Psi$ and
$M_2 = \sup\lvert\ddot\Psi\rvert$ along ZOH arcs, so $T^\star = 2\eta/M_2$.
Under (7.4), $\eta^s = s\,\eta_1$. For a plant with the dilation symmetry
(Dubins: $(t, \dot q, u) \to (t/s, s\dot q, \mathrm{diag}(s^k)u)$ leaves the
dynamics invariant) a ZOH arc of period $T$ of the scaled design is the base
arc of period $sT$, hence $M_2^s = s^2 M_{2,1}$ and

$$
T^\star(s) = \frac{T^\star_1}{s}, \qquad s_T = \frac{T^\star_1}{T}.
\tag{7.8}
$$

For the manipulator the frozen gravity compensation $G(q_k)$ differs from
$G(q(t))$ along the arc by a term of first order in $s$, so
$M_2^s = s^2 c_2 + s c_1$ and $T^\star(s) = 2\eta_1/(s c_2 + c_1)$: still
decreasing in $s$, $s_T$ by bisection, but with the ceiling
$T^\star \to 2\eta_1/c_1$ as $s \to 0$. The chosen $T$ must lie below that
ceiling.

**Design rule.** $s = \min(s_U, s_T)$; then $\delta$ as a fraction of $M_\Psi$,
$\varepsilon^s(T) = \varepsilon_1(sT)$ for the one-step increase (because
$\Psi^s = \Psi\circ D_s^{-1}$), and $\omega \ge \sup\ell/\varepsilon^s(T)$ as in
the Remark on the terminal cost.

### 7.4 What this changes with respect to (12) and (14)

- In (12), $\lambda$ is fixed (the program is positively homogeneous in
  $(k_1, \lambda, \delta_s)$, so minimising $\delta$ with $\lambda$ free is
  ill-posed and returns solver-tolerance values); the bound
  $\lvert k_{1,i}\rvert \le K$ is added so that the unconstrained program is
  bounded. (12) then fixes the **shape** of $k_1$; $s$ is the only scale.
- The input constraints of (14) are replaced by $s \le s_U$, which is exact
  on $\mathcal{X}_{RA}^{\delta}$, needs no sampled states and no sign
  condition on $\det A$. (14) is kept for systems whose zero-order drift
  alone is near the bounds, or when reshaping $k_1$ is wanted for strongly
  asymmetric bounds.
- The sampling condition and the input bounds act on the same variable, so
  Algorithm 1 gets one offline line: solve (12) at the normalised scale,
  compute $s_U$ and $s_T$, set $s = \min(s_U, s_T)$.
- The comparison with the unconstrained controller becomes a statement about
  scales: URA is the same shape at a performance scale $s_{\text{URA}} > s_U$,
  the terminal controller is at $s \le \min(s_U, s_T)$.

### 7.5 Numbers for the Dubins car (sampled suprema, base scale $s = 1$)

From `experiments/strict_xra/final/data/lemma2.pkl` and the terminal-set
inputs of `PAPER_DUBINS.md` Section 1.2, rescaled to $s = 1$:

| quantity | value |
| --- | --- |
| $\sup_{\mathcal{X}_{RA}^{\delta}}\lvert\omega_1\rvert$, $\sup\lvert a_1\rvert$ | 0.957 rad/s, 0.233 m/s$^2$ |
| $s_U$ from (7.7) with $(\omega_{\max}, a_{\max}) = (3, 3)$ | $\min(3.13, 3.59) = 3.1$ |
| $\eta_1$, $M_{2,1}$ | $1.26 \times 10^{-2}$, 3.15 |
| $T^\star_1 = 2\eta_1/M_{2,1}$ | 8.0 ms |
| $s_T$ from (7.8) with $T = 0.05$ s | 0.16 (0.146 with the certified rate $\lambda\delta$ in place of $\eta_1$) |
| terminal controller | $s = 0.159 = \min(s_U, s_T)$; the input bounds are inactive by a factor 20 |
| URA of Section 10 of `PAPER_DUBINS.md` | $s = 5 > s_U$: inputs up to 5.4 rad/s and 6.7 m/s$^2$ |

For the manipulator, with the parameters of `python/systems.py`
($m_i = 1$, $l_i = 4$, $l_{c,i} = 2$): $\max\lvert G_1\rvert = 78.5$ N m and
$\max\lvert G_2\rvert = 19.6$ N m against the bounds 680 and 500, so the static
condition holds with room, and (7.6) reduces to a quadratic inequality in
$s^2$ per state. Its $s_U$, $T^\star(s)$ and the ceiling $2\eta_1/c_1$ are the
first things to compute in the manipulator audit.

### 7.6 Code

`design/cert.py`, `scaled(s)` implements (7.4) exactly; `uc` is (7.3). $s_U$
and $s_T$ are currently sampled (`design/scan.py`, `final/run_lemma2.py`);
certified values would need SOS or interval bounds for the two suprema and
for $\eta_1$, $M_{2,1}$. `final/run_direct.m` shows that solving (12) directly
at the final $\lambda$ and $K$ reproduces the scaled $k_1$ within 6.5 %, with
a slightly larger certified set.

### 7.7 Solution flow under the scaling design (offline and online)

Inputs: system $(f, g, h)$, sets $\psi$ and $\phi$, input constraint
$P(x)u \le Q(x)$, sampling period $T$, horizon $N$, stage cost $\ell$.

| step | what is solved or computed | output | status in the code |
| --- | --- | --- | --- |
| 0. shape of $\psi$ | check that $\psi$ has its only critical point inside $\mathcal{X}_T$ on the relevant part of $\mathcal{X}_S$; otherwise redesign the weight | $\psi$ | done for Dubins (`design/screen_w3.py`) |
| 1. shape of $k_1$ | (12) with $\lambda = \lambda_0$ fixed and $\lvert k_{1,i}\rvert \le K_0$ on $\{\psi \ge 0\}$, solved in scaled outputs, degree $d$, objective $\min \delta_s$; a short scan picks the largest $\lambda_0$ at which the $K_0$ bound is still inactive | $k_1$, $\delta_s$, certified position share | `design/run_design.m` |
| 2. certificate | choose $\mu_0$ and $\delta = \theta M_\Psi$; build $\Psi$; check that $\mathcal{X}_{RA}^{\delta}$ stays away from the singularities of $A$ (Dubins: smallest speed on the set positive) | $\Psi$, $\mathcal{X}_{RA}$, $\mathcal{X}_{RA}^{\delta}$ | `design/cert.py`, `design/scan.py` |
| 3. input decomposition | from (7.5): the functions $u_k(x)$ of the base design | $u_k$ | implicit in `cert.py` (`uc` at any scale); symbolic form to add |
| 4. $s_U$ | (a) sample $\mathcal{X}_{RA}^{\delta}$, evaluate $\rho(x)$, take the minimum with a safety factor; (b) certify the chosen value by one SOS feasibility test in $x$ at fixed $s$, bisection if it fails | $s_U$ | (a) done (sampled); (b) to write (MATLAB, multipliers only) |
| 5. $s_T$ | $\eta_1 \ge \lambda_0\delta$ (certified by (12)); $M_{2,1}$ along ZOH arcs (sampled, or bounded); $T^\star_1 = 2\eta_1/M_{2,1}$; $s_T = T^\star_1/T$, or bisection on $T^\star(s) \ge T$ when the plant has a zero-order drift | $s_T$ | `final/run_lemma2.py` (sampled) |
| 6. scale | $s = \min(s_U, s_T)$, optionally times a safety factor; final design $(s k_1, s\lambda_0, s\delta_s, s^2\mu_0)$ | terminal controller | `cert.scaled(s)`, `final/config.py` |
| 7. verification at the final design | one-step increase $\varepsilon(T) > 0$ on $\mathcal{X}_{RA}^{\delta}$ (sampled, or $\eta T - \tfrac12 M_2 T^2 > 0$); inputs within the bounds on $\mathcal{X}_{RA}^{\delta}$; $\dot\Psi \ge \lambda\Psi$ spot check | pass or back to 6 with a smaller $s$ | `final/run_lemma2.py`, `ra_mpc_cert.terminal_stats` |
| 8. terminal cost | $\sup \ell$ on $\mathcal{X}_{RA}^{\delta} \setminus \mathcal{X}_T$, $\omega = c\,\sup\ell/\varepsilon(T)$ with $c > 1$, $M_\Psi = \max\psi - \delta_s/\lambda$ | $V_f$ | `ra_mpc_cert.terminal_stats` |
| 9. MPC | model with RK4 sub-steps, state constraint at the sub-step nodes with a back-off, terminal constraint $\Psi \ge \delta$, terminal cost $V_f$, horizon $N$ | problem (22) | `design/ra_mpc_cert.py` |
| online | Algorithm 1: solve, apply $u_{0\mid k}$, shifted candidate with $u^c(\hat x(NT\mid k))$, keep the cheaper, shrink $N$ when the candidate touches $\mathcal{X}_T$ | closed loop | `python/ra_mpc_dubins.py`, `ra_mpc_cert.py` |
| baseline | URA: the same $k_1$ at a performance scale $s_{\text{URA}}$ (no step 4 to 6) | URA | `final/config.py` (`URA_SCALE`) |

Remarks on the flow:

- Steps 1 to 3 are solved once; steps 4 to 6 involve only scalar searches and
  verification problems without decision variables other than multipliers.
- Only two sampled quantities would need certification for a fully rigorous
  design, $s_U$ (step 4b) and $M_{2,1}$ (step 5); everything else is either an
  SOS output or a closed form.
- The order of steps 4 and 5 does not matter; for the Dubins car step 5 is the
  binding one by a factor 20.
- If step 7 fails, reduce $s$ (or $T$); if the zero-order term $u_0$ alone
  violates the bounds in step 4, no $s$ works and the shape (14) or a smaller
  set is needed.
