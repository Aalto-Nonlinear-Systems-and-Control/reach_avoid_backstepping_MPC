# Strict X_RA test for the Dubins car (TAC revision)

Question: can the vanilla k1 SOS give a reach-avoid certificate that holds
rigorously (not only up to solver tolerance), and what does it cost?
Drivers run from this folder in MATLAB (it shadows the two SOS helpers of
`matlab/`); exports go to `controllers/tmp_*.py` and must be moved out
afterwards. `controllers/*.py` of the paper are never touched.

Helpers (copies of `matlab/` with test switches, all via globals):
`LMIN` lower bound on lambda, `KMAX` bound on |k1_i| on {psi >= 0}
(S-procedure), `SCALE` solve in ys = y / SCALE and map back.

## Findings

1. The original setting (lambda > 1e-8, degree 4) returns lambda = 1.1e-8,
   delta = 1.7e-10, i.e. a solution at the solver-tolerance scale. The SOS is
   homogeneous in (k1, lambda, delta, s1, s2), so delta/lambda should not
   depend on the scale; away from 1e-8 it converges to 10.8 = max psi, and at
   degree 8 the tiny-scale solve even returns delta < 0. The large X_RA of the
   paper is therefore a numerical artefact (consistent with Psidot < 0 on
   0.85 % of sampled states for the committed controller).
2. psi = alpha(y) h_raw(y) has a local maximum outside the target at
   (1.774, 0) with psi = 9.0742 (the Euclidean weight alpha is locally
   maximal there on the flat side of the superellipse ring). There
   grad psi = 0, so every k1 needs delta/lambda >= 9.0742: the largest
   certifiable set is {psi > 9.0742}, about 22 % of the safe-set area.
   (On an annulus some critical point besides the maximum is unavoidable;
   it must lie inside X_T to be harmless.)
3. With lambda fixed, |k1| is unbounded (solver returned |k1| ~ 2e4), so
   {Psi >= 0} is empty in the physical state box. Fix: |k1_i| <= K = 1.
4. Unscaled outputs (y up to 2.2, multipliers of degree 16) give violations
   of about 1e-3 even at O(1) scale. Solving in ys = y / 2 removes them.
5. Rigorous setting: SCALE = 2, K = 1, lambda = 0.02, degree 8:
   delta/lambda = 9.0743 (the limit of 2), k1 condition margin > 0 on a
   1501 x 1501 grid, |k1| <= 0.944, certified set 1.63 of 7.41 (area).

## Full-state certificate with the rigorous k1 (mu = 0.1)

| quantity | value |
| --- | --- |
| ECGBF in full state, min (Psidot - lambda Psi) on 9e3 samples of X_RA \ X_T | +1.9e-3 (no violation) |
| median Psidot / Psi | 0.85 (committed controller: 2e-5) |
| X_RA share of the synthesis box / of the 300 MPC candidates | 2.6 % / 11 of 300 |
| u^c within the bounds 5 on X_RA | 95-96 % (max abs omega 18.8) |
| continuous u^c from 40 states | 40 reached, Psi nondecreasing in all, median 3.1 s |
| ZOH one step, share with Psi decreasing | T = 0.05: 41 %, 0.01: 9.5 %, 0.005: 2.7 %, 0.002: 0.4 %, 0.001: 0.1 % |
| Lemma 2 bound T*(delta = 0.1 max Psi) | 1.1e-7 s (was 1e-14 s) |

mu = 1 enlarges X_RA (54 of 300 candidates) but u^c meets the bounds on only
16 % of it. The SOP stage (`solve_k1_controller_sop`, unscaled) reintroduces
small violations (0.2 % of samples); it needs the same scaling.
