# All eight kitting tasks (square layout, common torque bounds (28.5, 8.5) N m, T = 10 ms, N = 60, 150 initial states each)

## Design and scale

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

## Comparison

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

## Theorem 2 / Corollary 1 and Lemma 2

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
