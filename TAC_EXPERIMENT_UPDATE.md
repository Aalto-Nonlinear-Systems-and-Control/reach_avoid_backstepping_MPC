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
  largest $T$ with $\beta(T)/\lambda < \delta$.
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
  bounds.

### 1.5 Cost comparison rule (suboptimal MPC)

- Paper: the corollary only needs the applied sequence to cost no more than
  the candidate.
- Code now: the solver output is used as is.
- Change: warm start at the candidate and compare the returned cost with the
  candidate cost; keep the cheaper one.

### 1.6 Fallback branch (Algorithm 1, lines 8 and 9)

- Paper: if the candidate's predicted trajectory enters $\mathcal{X}_T$,
  execute it (current text: open loop until entry).
- Code now: an infeasible solve marks the trajectory as failed.
- Change: add the fallback branch. See decision D1 for open loop versus a
  shrinking-horizon reach problem.

### 1.7 Time bound (Corollary 1)

- Paper: entry into $\mathcal{X}_T$ no later than
  $t = (\lfloor J_0^*/\underline{\ell} \rfloor + N + 1)\,T$; the proof uses
  $J^*_{k+1} \le J^*_k - \underline{\ell}$.
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

- [ ] cell 4, `solve_mpc`: terminal constraint `cert_vec(X[:, N]) >= delta`.
- [ ] cell 4, cost: drop the `Qf_y` term; add `omega * (M_psi - cert_vec(X[:, N]))`.
- [ ] cell 4, warm start: shifted candidate whose last input is `u_c(x_N_pred)`
      (new numpy callable from `u_opt`); roll it out with the same RK4 map
      `F_disc` (the current rollout uses explicit Euler).
- [ ] cell 4, return value: also return the solver cost; add a
      `candidate_cost(x0, U)` helper and apply the "keep the cheaper" rule.
- [ ] cells 5 and 7, `_feasible_receding` and `simulate_ra_mpc_traj`: check
      whether the candidate rollout enters $\mathcal{X}_T$ within the horizon
      and add the fallback branch. An infeasible solve outside the fallback is
      now a violation; log it separately from "failed".
- [ ] log $J^*_k$ and check $J^*_{k+1} \le J^*_k - \underline{\ell}$.
- [ ] header comment and printouts (they still say `V(x_N) >= 0`).

### `notebooks/example_manipulator_reach_avoid_mpc.ipynb`

- [ ] cell 3: the same four changes in `solve_mpc` and in the shared
      pre-compiled solver (`solve_mpc_fast`, `make_fast_solver`).
- [ ] cell 4: the same closed-loop changes.

### `notebooks/make_paper_figs_dubins.ipynb` and `notebooks/make_paper_figs_manipulator.ipynb`

- [ ] `feas_rampc` and `make_fast_solver(..., "cert")`: the same MPC changes.
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

- [ ] add `python/ra_mpc_constants.py` (Section 2).
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
      $(\lfloor J_0^*/\underline{\ell} \rfloor + N + 1)\,T$ on every successful run.
- [ ] no infeasible solve outside the fallback branch.
- [ ] success rates and closed-loop plots regenerated; numbers in the paper
      (Fig. 4 and Fig. 6 counts, sweep curves, computation-time table) updated.

## 5. Open decisions

- D1, fallback implementation. Paper Algorithm 1 currently executes the
  stored candidate open loop until entry. Alternative discussed: keep
  re-solving a shrinking-horizon reach problem with terminal constraint
  $\hat{x}(\tau_e - kT) \in \mathcal{X}_T$, warm-started at the candidate,
  which retains feedback. Code and paper must agree.
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
