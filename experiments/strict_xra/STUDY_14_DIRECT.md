# Study: program (14) solved directly against (12) plus time scaling (top-view arm, 2026-10-10)

Question: is it better to obtain $k_1$ by solving (14) directly, with the input bounds enforced at sampled states, than to
solve (12) once and choose the time scale $s$ afterwards (the method behind all reported data)?

## Setting

- **A, current method.** (12) at $\lambda_0 = 0.1$ with $\lvert k_{1,i}\rvert \le 1$ on $\{\psi \ge 0\}$, scaled by $s$: $k_1 \to s k_1$, $\lambda = s\lambda_0$, $\mu = s^2\mu_0$ ($\mu_0 = 0.1$). The deployed $s = \min(s_{\rm bnd}, 0.8\, s_\varepsilon)$.
- **B, (14) directly.** The same $\lambda = s\lambda_0$ and $\mu = s^2\mu_0$ are fixed.
  - The ECGBF condition of (12) is kept as an SOS constraint, with $\lvert k_{1,i}\rvert \le s$ on $\{\psi \ge 0\}$; this is the same box as A.
  - The input bounds $\lvert u^c_i(x_j)\rvert \le U_i$ are imposed at states $x_j$ of the terminal set $\mathcal{X}_{RA}^{\delta}$. $u^c(x_j)$ is affine in the coefficients of $k_1$ (Theorem 1), so these are linear constraints in the same SOS program (`arm_top/solve_k1_sampled14.m`).
  - The terminal set depends on $k_1$, so B runs an exchange loop. Each round solves the program, samples about 20 000 states of the new design's own terminal set, and adds the 250 states closest to or beyond the bounds. The loop stops when a fresh sample is within the bounds (`arm_top/study14.py`).
  - Same degree 8, scaled outputs $y/4$, same solver.
- **Acceptance of a design.** Each design must pass three checks, all of which A also meets at its deployed $s$:
  - within the bounds on a second fresh terminal sample;
  - ECGBF margin $\ge 0$ on a $601^2$ grid;
  - sampling-period threshold $T_\varepsilon \ge T/0.8 = 12.5$ ms, i.e. A's 0.8 safety factor.
- **Common conditions.** Common torque bounds (28.5, 8.5) N m, $T = 10$ ms, $N = 60$. RA-MPC is run on each task's 150-state common initial set; everything except the terminal design is unchanged (`arm_top/eval14_mpc.py`).

## Designs (max $\lvert u_i\rvert/U_i$ on the terminal set, sampling-period threshold)

| task (A's limit) | $s$ | A: ratio, $T_\varepsilon$ [ms] | B: ratio, $T_\varepsilon$ [ms], grid margin | B rounds, time | B accepted |
| --- | --- | --- | --- | --- | --- |
| blue_fwd (bounds, $s^\ast = 1.747$) | 1.747 | 0.80, 16.4 | 0.69, 16.1, +1.7e-3 | 1, 17 s | yes |
| | 2.0 | 1.04, 14.4 | 0.81, 14.6, +2.2e-3 | 1, 17 s | yes |
| | 2.3 | 1.38, 12.5 | 0.93, 14.2, +3.0e-3 | 1, 17 s | yes |
| | 2.6 | 1.76, 11.1 | 0.97, 16.1, +3.4e-3 | 4, 229 s | yes |
| | 2.9 | 2.19, 9.9 | 0.98, 12.9, +4.9e-3 | 4, 231 s | yes |
| | 3.2 | 2.67, 9.0 | infeasible after 4 rounds | 5, 394 s | no |
| green_fwd (bounds, $s^\ast = 1.623$) | 1.623 | 0.80, 16.0 | 0.66, 17.6, +1.2e-3 | 1 | yes |
| | 1.85 | 1.05, 14.0 | 0.74, 14.5, +6.4e-4 | 1 | yes |
| | 2.074 | 1.31, 12.5 | 0.80, 13.3, +3.2e-3 | 1 | yes |
| | 2.3 | 1.61, 11.3 | 0.94, 11.7, **−8.1e-4** | 1 | no (margin, period) |
| orange_fwd (bounds, $s^\ast = 1.690$) | 1.69 | 0.80, 14.1 | 0.71, 14.2, +7.7e-4 | 1 | yes |
| | 1.8 | 0.91, 13.2 | 0.78, 13.5, +8.5e-4 | 1 | yes |
| | 1.901 | 1.02, 12.5 | 0.83, 12.9, +9.1e-4 | 1 | yes |
| | 2.1 | 1.24, 11.3 | 0.94, 12.3, +1.1e-3 | 1 | no (period) |
| red_fwd (period, $s^\ast = 2.386$) | 2.386 | 0.80, 12.5 | 0.75, 11.9, +4.2e-3 | 1 | no (period) |
| | 2.669 | 1.01, 11.2 | 0.91, 10.4, +4.7e-3 | 1 | no (period) |
| | 3.0 | 1.27, 9.9 | 0.99, 9.0, +1.9e-3 | 4 | no (period) |

**Without the bound on $k_1$** (blue, $s$ = 2.6 and 3.0) the loop diverges. $k_1$ grows where there are no samples, and the terminal set moves there. Fresh samples then reach 89 and 650 times the bounds, and after two rounds the terminal set is empty. (14) as printed, with $\lambda$ free and no bound on $k_1$, cannot be solved this way. The sampled input constraint also contains $\lambda\,k_1$, which is bilinear unless $\lambda$ is fixed.

Green's design A has the same numerical difficulty. Its SOS solve reports a Mosek feasratio of 3.2, and the B solves report about 1.5. Its grid margin is positive (+2.4e-3, the same on $601^2$, $1201^2$ and $2401^2$ grids), so it remains a valid certificate.

## RA-MPC on the common initial sets (150 states each)

| task | A ($s^\ast$) | B ($s$) | reached A / B | same states | median time A / B [s] | median cost A / B | Theorem 2 checks B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| blue_fwd | 1.747 | 1.747 | 137 / 137 | yes | 0.94 / 0.94 | 17.8 / 17.8 | all pass |
| blue_fwd | 1.747 | **2.9** | 137 / 137 | yes | 0.94 / 0.87 | 17.8 / 18.3 | all pass |
| green_fwd | 1.623 | **2.074** | 144 / 144 | yes | 0.77 / 0.74 | 12.2 / 11.9 | all pass |
| orange_fwd | 1.690 | **1.901** | 138 / 138 | yes | 0.92 / 0.90 | 19.8 / 19.9 | all pass |
| red_fwd | 2.386 | 2.386 | 137 / 137 | yes | 0.65 / 0.65 | 9.8 / 9.8 | all pass |

"All pass" means:
- no infeasible solve after $k = 0$;
- no decrease violation;
- horizon nonincreasing;
- $u^c$ within the bounds at every candidate;
- smallest $\psi$ along the closed loops 0.0045 to 0.0054.

The states that fail are the same for A and B: 13, 6, 12 and 13 states that are infeasible at $k = 0$. They fail because of the bounds and the corridor, not the terminal set; they are the states that the URA reaches only by exceeding the bounds.

## Conclusions

1. **Where the bounds bind, (14) solved directly gives a faster terminal controller under the same bounds and the same period rule.**
   - Blue: 1.66 times faster (s up to 2.9 against 1.747).
   - Green: 1.28 times faster.
   - Orange: 1.12 times faster.
   - It reshapes $k_1$ so that $u^c$ uses about 20 % less torque at equal speed. The sampling period then becomes the limit.
2. **Where the period binds, it gives nothing.** That is the case for red, the paper's task, and for the Dubins example, whose bounds have a slack of 2.3 to 4 at the deployed $s$. Program (14) does not address the period, and B's period threshold is equal or slightly lower.
3. **RA-MPC is unchanged in every case.**
   - The same states are reached.
   - The median time to the target is up to 7 % shorter.
   - The cost is within ±3 %.
   - A faster terminal controller does not enlarge what the MPC can reach within its horizon under the bounds.
4. **Cost and side effects of B.**
   - B costs one SOS solve per round and 1 to 4 rounds per speed (17 s to 4 min), plus a search over speeds. A needs one solve (about 1 s) and a one-dimensional choice of $s$.
   - B loses the exact scaling structure. That structure gives the closed-form relation $T_\varepsilon(s) = T_\varepsilon(1)/s$ and defines the URA baseline as the same $k_1$ at $2s$.
   - Both methods certify the bounds by sampling.
5. **For the paper.**
   - Keep (12) with $\lambda$ fixed and the bound on $k_1$, followed by the choice of $s$, as the method.
   - (14) can be kept as a remark. The input bounds can be added at sampled terminal states, which is convex because $u^c$ is affine in $\Theta$ for fixed $\lambda$ and $\mu$, and this helps only when the bounds rather than the period limit the terminal controller.
   - In that form (14) needs $\lambda$ fixed and the bound on $k_1$, and its sampled states must come from the new design's own terminal set, through the exchange loop.

## Files

- `arm_top/solve_k1_sampled14.m`, `arm_top/study14.py`, `arm_top/eval14_mpc.py`, `arm_top/run_study14.sh`, `arm_top/run_eval14.sh`.
- Designs: `arm_top/k1/<task>_s14_s<s>_k1.txt`.
- Logs and results: `arm_top/data/study14/` (`<task>_s<s>_k1.json` per design, `mpc_<task>_s<s>_k1.pkl` with all closed loops, `sweep.log`, `sweep2.log`, `mpc.log`, `mpc2.log`).
- No existing design or result file was changed.
