# TAC revision — experiment update plan (DRAFT for confirmation)

Branch: `tac-revision` (created from `lcss-revision`, which produced the
figures currently in the TAC manuscript). This file lists what must change
in the code so that the experiments match the revised Section IV of the
paper (`TAC/method.tex`: problem (21), terminal-cost condition (22),
Algorithm 1, Corollary 1, Remark 4). Nothing below has been implemented
yet; tick the boxes as items are confirmed.

Hard constraint carried over from the LCSS revision: **do not overwrite
`controllers/*.py`** (they are the certified controllers behind the reported
numbers). Any MATLAB re-run must use throwaway export names.

---

## 1. Paper changes and their code implications

| Paper (TAC/method.tex) | Current code | Required change |
| --- | --- | --- |
| Terminal set is the shrunk set $\mathcal X_{RA}^{\delta}=\{\Psi\ge\delta\}$ with $\delta>\beta(T)/\lambda$ (sampling condition of Lemma 2 holds at every terminal state) | `cert_vec(x_N) >= 0` (zero-superlevel set) | `cert_vec(x_N) >= delta`, with `delta` computed from the constants in Sec. 2 |
| Cost $J=\sum_{i<N}\ell(\hat x_i,u_i)+V_f(\hat x_N)$ with $\ell\ge\underline{\ell}>0$ outside $\mathcal X_T$ and $V_f$ satisfying the one-step decrease (22) along ZOH $u^c$ | stage: $\|h(x)-y^*\|_{Q}^2+\|u\|_R^2$; terminal: $\|h(x_N)-y^*\|_{Q_f}^2$ | keep the stage cost (its lower bound outside $\mathcal X_T$ is automatic, see Sec. 2); **replace** the quadratic terminal term by $V_f(x)=\omega\,(M_\Psi-\Psi(x))$ (Remark 4). Adding $V_f$ on top of $Q_f$ is not allowed: the tracking term need not decrease along $u^c$, which would break (22) |
| Fixed sampling period $T$ (= `dt_mpc`), with $T<T^*(\delta)$ | fixed `dt_mpc` | unchanged, but $\delta$ and $\omega$ depend on $T$; every `dt` in the sweeps needs its own $\delta(T)$, $\omega(T)$ |
| Algorithm 1 line 7: shifted candidate $\mathbf U_{k+1}=\{u^*_{1|k},\dots,u^*_{N-1|k},u^c(\hat x^*(NT|k))\}$ | warm start repeats the last input (`U_warm[:, -1:]`) | build the candidate with $u^c$ evaluated at the predicted terminal state (numpy `lambdify` of `u_opt`); verify it respects the input bounds |
| Implementation rule behind the corollary (suboptimal MPC): the applied sequence must not cost more than the candidate | solver output used as is | warm start at the candidate **and** compare the returned cost with the candidate cost; keep the cheaper one |
| Algorithm 1 lines 8–9: if the candidate's predicted trajectory enters $\mathcal X_T$, execute it (paper text: open loop until entry) | infeasible solve ⇒ trajectory marked failed | add the fallback branch; see open decision D1 for open-loop vs. shrinking-horizon reach problem |
| Corollary 1: entry into $\mathcal X_T$ no later than $t=(\lfloor J_0^*/\underline{\ell}\rfloor+N+1)T$; proof uses $J^*_{k+1}\le J^*_k-\underline{\ell}$ | not logged | log $J_k^*$ per closed-loop step; check monotone decrease and the bound (validation only, not a figure unless we decide so) |

---

## 2. Constants needed per system (all currently missing on the Python side)

Notation as in the paper; "over $\mathcal X$" means over the compact
sampling box used in MATLAB (`bound_min`/`bound_max` in `example_*.m`) unless
we decide to use $\{\Psi\ge 0\}$ instead (see D3).

| Symbol | Meaning | How to get it |
| --- | --- | --- |
| $\lambda$ | ECGBF rate | Numeric value `k1_lambda` inside `solvesop_bounded_control.m`; **not exported**. Options: (a) re-run `example_*.m` with throwaway export names and print `k1_lambda`; (b) estimate numerically as $\min$ of $\dot\Psi/\Psi$ under $u^c$ over samples with $\Psi\ge\epsilon$ (verification, conservative). Prefer (a) |
| $M_{g_j}$ | $\max_{\mathcal X}\lvert L_{g_j}\Psi\rvert$ | sympy gradient of `certificate_opt` times column $g_j$, sampled |
| $\ell_{u^c_j}$ | Lipschitz constant of $u^c_j$ | $\max_{\mathcal X}\lVert\nabla u^c_j\rVert_2$, sampled |
| $\ell_f,\ \ell_{g_j}$ | Lipschitz constants of $f$, $g_j$ | $\max_{\mathcal X}\lVert\partial f/\partial x\rVert_2$ etc. Dubins: $g$ constant so $\ell_{g_j}=0$; manipulator: $g=M(q)^{-1}$ columns, nonzero |
| $M_{u^c_j}$ | $\max_{\mathcal X}\lvert u^c_j\rvert$ | sampled (should be $\le$ input bounds on $\{\Psi\ge0\}$) |
| $\ell_{sys}$ | $\ell_f+\sum_j\ell_{g_j}M_{u^c_j}$ | formula |
| $M_v$ | $\max_{\mathcal X}\lVert f+\sum_j g_ju^c_j\rVert$ | sampled |
| $\gamma(T)$ | $(\sum_j M_{g_j}\ell_{u^c_j})\,(e^{\ell_{sys}T}-1)/\ell_{sys}$ | formula |
| $\beta(T)$ | $\gamma(T)\,M_v$ | formula |
| $\delta$ | terminal margin, $\delta>\beta(T)/\lambda$ | choice, see D2 |
| $M_\Psi$ | upper bound of $\Psi$ on $\overline{\mathcal X_{RA}}$ | sampled max of `certificate_opt` on $\{\Psi\ge0\}$, inflated slightly |
| $\sup\ell$ | max stage cost on $\mathcal X\times\mathcal U$ | $Q\max\lVert h(x)-y^*\rVert^2+R\,u_{\max}^2$, sampled |
| $\underline{\ell}$ | stage-cost lower bound outside $\mathcal X_T$ | automatic from the tracking term: Dubins target $y_2^2+4(y_1+1.7)^2\le0.4$ gives $\lVert h(x)-y^*\rVert^2>0.1$ outside, so $\underline{\ell}=0.5$ with $Q=5I$; manipulator ellipse gives $\lVert h(x)-y^*\rVert^2>0.32$, so $\underline{\ell}=1.6$. Verify numerically |
| $\omega$ | terminal-cost weight | $\omega\ge\sup\ell/[(e^{\lambda T}-1)(\delta-\beta(T)/\lambda)]$ (Remark 4); use the minimal admissible value |

Caveat on the exported certificate: `certificate_opt` is already shifted by
$-\delta_{SOS}/\lambda$ in MATLAB (`reach_avoid_controller.m`), so it is the
function that satisfies $\dot\Psi\ge\lambda\Psi$ exactly and is the paper's
$\Psi$. The SOS slack $\delta_{SOS}$ is unrelated to the terminal margin
$\delta$ above (same symbol clash as in the paper).

Risk: $\beta(T)$ uses global Lipschitz constants and may be conservative. If
$\beta(T)/\lambda$ comes out close to or above $M_\Psi$ at $T=0.05$ s
(Dubins) or $0.01$ s (manipulator), the terminal set is empty and the
theory demands a smaller $T$. Compute the numbers first; this is the gate
for everything else.

Proposed home: new module `python/ra_mpc_constants.py` with one function
per system returning a dict of the constants plus `delta(T)`, `omega(T)`,
and a notebook cell that prints them.

---

## 3. File-by-file changes

### `notebooks/example_dubins_car_reach_avoid_mpc.ipynb`
- [ ] cell 4 `solve_mpc`: terminal constraint `cert_vec(X[:, N]) >= delta`.
- [ ] cell 4 cost: drop `Qf_y` term; add `omega * (M_psi - cert_vec(X[:, N]))`.
- [ ] cell 4 warm start: shifted candidate with last input `u_c(x_N_pred)`
      (new numpy callable from `u_opt`); rollout with the same RK4 map
      `F_disc` (currently uses explicit Euler for the rollout).
- [ ] cell 4 return value: also return the solver cost and expose a
      `candidate_cost(x0, U)` helper; apply the "keep the cheaper" rule.
- [ ] cells 5 and 7 (`_feasible_receding`, `simulate_ra_mpc_traj`):
      candidate check (does the candidate rollout enter $\mathcal X_T$ within
      the horizon?) and fallback branch; infeasibility outside the fallback
      is now a violation, log it separately from "failed".
- [ ] log $J_k^*$ and check $J^*_{k+1}\le J^*_k-\underline{\ell}$.
- [ ] header comment and printouts (currently say "V(x_N) >= 0").

### `notebooks/example_manipulator_reach_avoid_mpc.ipynb`
- [ ] cell 3: same four changes in `solve_mpc` and in the shared
      pre-compiled solver (`solve_mpc_fast` / `make_fast_solver`).
- [ ] cell 4: same closed-loop changes.

### `notebooks/make_paper_figs_dubins.ipynb`, `notebooks/make_paper_figs_manipulator.ipynb`
- [ ] `feas_rampc` / `make_fast_solver(..., "cert")`: same MPC changes.
- [ ] `dt` sweep: recompute `delta(dt)` and `omega(dt)` for every `dt` in
      `DT_SWEEP`; `_sweep_sig()` must include the constants (or a version
      tag) so the caches in `figures/*.npz` are invalidated.
- [ ] expect the sweep to change qualitatively: for the larger `dt` values
      the sampling condition may be violated; decide how the paper reports
      this (D4).

### `notebooks/noise_robustness_dubins.ipynb`
- [ ] builds its own RA-MPC solver: same changes.

### `notebooks/example_*_vanilla_mpc.ipynb`, `notebooks/example_*_unconstrained_reach_avoid.ipynb`
- No change (baselines).

### `python/`
- [ ] add `python/ra_mpc_constants.py` (Sec. 2).
- [ ] optional (D5): factor the RA-MPC builder into `python/ra_mpc.py` so
      the four notebooks share one implementation.

### `matlab/`
- [ ] `export_to_python.m`: also write `lambda_val` (and `mu_val`) into the
      exported file. Only run with throwaway export names; then copy the
      scalar by hand into `ra_mpc_constants.py` or a small JSON next to the
      controller file. Do not regenerate `controllers/*.py`.

### `README.md`
- [ ] "Controller Settings Per Notebook": terminal constraint column
      $V(x_N)\ge\delta$, add $\delta$, $\omega$ columns, drop $Q_{f,y}$ for
      RA-MPC rows.
- [ ] "MPC Formulation": new cost and terminal set; mention the candidate
      warm start, the cost comparison rule, and the fallback.

### `data/`, `figures/`
- [ ] regenerate `traj_controls_*.npz`, `*_cache.npz`, and the PNG/PDF
      figures; copy the new figures to `TAC/figs/` in the paper repo.

---

## 4. Validation checklist (per system, before regenerating figures)

- [ ] constants printed and sanity-checked; $\beta(T)/\lambda<M_\Psi$ with
      room for $\delta$.
- [ ] $\{\Psi\ge\delta\}$ non-empty and still reachable within $N$ steps from
      a comparable share of initial states (compare feasibility counts
      before/after).
- [ ] candidate $u^c(\hat x^*(NT|k))$ within the input bounds at every step.
- [ ] $J_k^*$ decreases by at least $\underline{\ell}$ at every solved step.
- [ ] reaching time $\le(\lfloor J_0^*/\underline{\ell}\rfloor+N+1)T$ on every
      successful run.
- [ ] no infeasible solve outside the fallback branch.
- [ ] success rates and closed-loop plots regenerated; numbers in the paper
      (Fig. 4/6 counts, sweep curves, Table of computation times) updated.

---

## 5. Open decisions (need confirmation)

- **D1 Fallback implementation.** Paper Algorithm 1 currently executes the
  stored candidate open loop until entry. Alternative discussed: keep
  re-solving a shrinking-horizon reach problem (terminal constraint
  $\hat x(\tau_e-kT)\in\mathcal X_T$) warm-started at the candidate, which
  retains feedback. Code and paper must agree.
- **D2 Choice of $\delta$.** Candidates: $\delta=2\beta(T)/\lambda$, or
  $\delta=\beta(T)/\lambda+\theta M_\Psi$ with a small $\theta$. Larger
  $\delta$ shrinks the terminal set but enlarges the guaranteed increment
  and reduces $\omega$.
- **D3 Domain for the constants.** Sampling box $\mathcal X$ from MATLAB vs
  the tighter $\{\Psi\ge0\}$. The paper's constants are stated on $\mathcal X$.
- **D4 Sweep semantics.** If $\delta(dt)$ exceeds $M_\Psi$ for the coarse
  `dt` values, report $\rho=0$ there, or drop those points and say so.
- **D5 Refactor.** Shared `python/ra_mpc.py` vs editing four notebooks in
  place.
- **D6 $\omega$ magnitude.** The theoretical minimum may dominate the
  tracking cost. Report the value; if the closed-loop behaviour degrades,
  discuss whether to present the tracking cost as illustrative in the paper.
