# Two-link arm: audit of the paper's setting and proposal (2026-10-07)

Scripts in `arm/`: `screen_psi_arm.py` (safe-set functions), `run_design_arm.m`
(k1 synthesis), `cert_arm.py` (dynamics, controller, certificate, scaling,
sampling, ZOH), `scan_arm.py`, `m2_split.py`, `scan_arm2.py`. The MPC and the
three-method comparison are not built yet; they follow the Dubins pipeline
once the points of Section 5 are decided.

## 1. What the paper's example uses

| item | value | source |
| --- | --- | --- |
| arm | $m_i = 1$ kg, $l_i = 4$ m, $l_{c,i} = 2$ m, $I_i = 0.02$, $g = 9.81$ | `python/systems.py`, `matlab/example_manipulator.m` |
| output | end-effector position $y = h(q)$, elbow-down branch ($q_2 > 0$) | notebooks |
| safe set | $b(y) = -(4(y_1 - 2) - 2y_2^3)^2 + 0.8\,y_2^3 + 10 > 0$ | paper, code |
| target | $\phi(y) = (y_1 - 5.5)^2/1.2^2 + (y_2 - 1.8)^2/0.4^2 - 2 < 0$ | paper, code |
| bounds | $\lvert\tau_1\rvert \le 680$, $\lvert\tau_2\rvert \le 500$ N m | paper |
| certificate | $\mu = 15$, $k_1$ of degree 2, $\lambda$ free in the SOS (lower bound $10^{-8}$) | `example_manipulator.m` |
| MPC | $T = 0.01$ s, $N = 20$ (horizon 0.2 s) | paper |
| initial states | 100 to 150 positions in the band, at rest, elbow-down | notebooks |
| "URA success" in the figures | membership of the state in the URA's certified set, not a closed-loop run | `make_paper_figs_manipulator.ipynb` |

## 2. Problems found

1. **$b$ has no maximum and the band leaves the workspace.** Along its ridge
   $y_1 = 2 + y_2^3/2$ the function is $0.8 y_2^3 + 10$, increasing without
   bound; its only critical point is a saddle at $(2, 0)$ with
   $b = 10$. Inside the reach $\lvert y\rvert < 8$ the largest value, 19.1, is at
   $(7.67, 2.25)$ outside the target (largest value inside the target 17.8),
   and the band crosses the singular circle $\lvert y\rvert = 8$ ($\det J = 0$) at
   $y_2 \approx 2.25$. The assumption of a unique maximum in
   $\mathcal{X}_T \cap \mathcal{X}_S$ fails, every certified field has to push
   states beyond the target towards the singularity, and the SOS condition is
   imposed on the unreachable part of the band as well.
2. **The certificate is a solver-tolerance artefact**, as for the Dubins car:
   $\lambda$ is free with lower bound $10^{-8}$ and the objective drives it there.
3. **The torque bounds are far above what the arm needs.** Gravity needs at
   most 78.5 and 19.6 N m; moving the end effector at 1 m/s$^2$ needs about
   10 N m. With 680 and 500 N m the MPC may accelerate at 60 m/s$^2$ and no
   reach-avoid controller of reasonable speed comes near the bounds.
4. **URA success was set membership.** The closed loop of the old controller
   was not run.
5. **Sampling period.** With a held torque the frozen gravity compensation
   gives $\Psi$ a second-derivative term that does not shrink with the time
   scale (Section 4): no scale makes $T = 0.01$ s certifiable, and the sampled
   one-step increase already fails at 5 ms.

## 3. Redesign of the safe set

Screened $\psi = w(y)\,b(y) - c\lvert y\rvert^4$ over $c \in \{0.005, 0.01, 0.02\}$
and $w \in \{1,\ 1 + \kappa(y_2 + 2.5)\}$, $\kappa \in \{0.3, 0.6, 1\}$
(`arm/screen_psi_arm.py`). The $\lvert y\rvert^4$ term closes the band inside the
reach; the weight removes the saddle by making $\psi$ increase along the
ridge up to the target. Only one candidate satisfies everything:

$$
\psi(y) = \big[(1 + 0.6\,(y_2 + 2.5))\,b(y) - 0.02\,\lvert y\rvert^4\big] / 39.04 ,
$$

normalised to maximum 1. Properties: safe set with $\lvert y\rvert \in [1.21, 7.74]$
($q_2 \in [30^\circ, 163^\circ]$, $\lvert\sin q_2\rvert \ge 0.3$, so the decoupling
matrix stays away from singularity), unique maximum at $(4.21, 1.66)$ inside
the target, no other critical point, 28 % of the target ellipse inside the
safe set. The old boundary and the new set are in `arm/figs/arm_sets.png`.

Synthesis of $k_1$ with (12) as for the Dubins car (scaled outputs $y/4$,
$\lambda$ fixed, $\lvert k_{1,i}\rvert \le 1$ m/s, degree 8, `arm/run_design_arm.m`):

| $\lambda_0$ | $\delta_s/\lambda_0$ | certified share of the band | $\max\lvert k_{1,i}\rvert$ | note |
| --- | --- | --- | --- | --- |
| 0.02 | 0.000 | 100 % | 0.76 m/s | bound inactive; ECGBF margin on the grid $+4.5 \times 10^{-4}$ |
| 0.03 | 0.000 | 100 % | | |
| 0.04 | 0.123 | 88 % | 1.00 m/s | bound active |
| 0.05 | 0.298 | 61 % | 1.00 m/s | |
| 0.10 | 0.646 | 18 % | 1.00 m/s | |

The output-space flow $\dot y = k_1(y)$ of the 0.02 and 0.04 designs reaches
the target from 300 of 300 band positions. $\lambda_0 = 0.04$ is used below
(larger rate, 88 % of the band).

## 4. Sampling period and gravity

For the arm the input of the scaled design is $\tau_s = G(q) + s^2\tilde\tau_1(x)$
(Section 7 of `TAC_EXPERIMENT_UPDATE.md`). Under a zero-order hold of the
torque the compensation $G(q_k)$ is frozen while $q$ moves, which adds to
$\ddot\Psi$ a term that does not scale: $M_2(s) = c_0 + c_2 s^2$ with
$c_0 = 8$ to $24$ (`arm/m2_split.py`). The second-order bound
$T^\star = 2 s\eta_1/(c_0 + c_2 s^2)$ is then largest at $s^\star = \sqrt{c_0/c_2}$
and never exceeds $\eta_1/\sqrt{c_0 c_2}$:

| $\lambda_0$, $\mu_0$, $\delta$ | $c_0$ | $c_2$ | best $T^\star$ |
| --- | --- | --- | --- |
| 0.02, 0.1, $0.2 M_\Psi$ | 19.4 | 5.2 | 0.5 ms |
| 0.04, 0.1, $0.5 M_\Psi$ | 17.9 | 3.1 | 2.8 ms |
| 0.04, 0.3, $0.5 M_\Psi$ | 10.9 | 10.4 | 1.9 ms |

The sampled one-step increase confirms it: negative already at $T = 5$ ms for
every $\mu_0$ at $s = 1$ (`arm/scan_arm.py`). Two ways out:

- **A. Keep the held torque and use $T$ of 1 to 2 ms** (horizon of 0.2 to
  0.5 s needs $N$ of 100 to 500).
- **B. Compensate gravity continuously.** The input held by the sampled-data
  loop is $v = \tau - G(q)$, the plant receives $v + G(q(t))$; the bound becomes
  the state-dependent affine constraint $\lvert G(q) + v\rvert \le \tau_{\max}$,
  which is the form (13) of the paper allows. Then $c_0 = 0$ exactly and the
  scaling rule of the Dubins car applies (`arm/scan_arm2.py`, $T = 0.01$ s):

| $\lambda_0 = 0.04$, $\mu_0$, $\delta$ | terminal-set positions, share of the band | $T^\star_1$ | $s_\varepsilon(T = 0.01)$ sampled |
| --- | --- | --- | --- |
| 0.1, $0.3 M_\Psi$ | 57 % | | 0.56 |
| 0.1, $0.5 M_\Psi$ | 32 % | 11.0 ms | 1.10 |
| 0.3, $0.3 M_\Psi$ | 57 % | | 0.19 |
| 0.3, $0.5 M_\Psi$ | 32 % | 3.7 ms | 0.48 |

With B and $\mu_0 = 0.1$, $\delta = 0.5 M_\Psi$ the terminal controller can run
at $s \approx 0.8$ for $T = 0.01$ s, and even the second-order bound certifies
$s \le 1.1$.

## 5. Points to decide

1. **Gravity handling: A or B.** B is recommended; it needs one sentence in
   the example and makes the state-dependent input constraint of the paper
   meaningful.
2. **Torque bounds.** Candidates well above the gravity load, e.g.
   $(\tau_{1,\max}, \tau_{2,\max}) = (120, 40)$ or $(150, 50)$ N m. The terminal
   controller is admissible up to $s_U = 5.6$ or $6.7$ ($\mu_0 = 0.1$, $\delta = 0.5 M_\Psi$).
   URA from rest at $s = 3$ has peaks of 66 to 98 N m in $\tau_1$, below all
   candidates; a URA that exceeds the bounds needs $s \ge 5$ (end-effector
   speeds above 4 m/s) or bounds close to the gravity load, where the MPC
   has little torque left for motion. The contrast of the Dubins figure is
   not available here at reasonable speeds; the comparison would rest on
   time, cost and the domain of attraction.
3. **$\mu_0$.** 0.1 gives the faster terminal controller ($s_\varepsilon = 1.1$)
   but no rest state inside $\mathcal{X}_{RA}$ at band positions (URA from
   rest still reaches 88 of 150 positions); 0.3 puts 47 of 150 rest states
   inside $\mathcal{X}_{RA}$ (URA reaches 128) at the price of $s_\varepsilon = 0.48$.
4. **Initial states.** At rest at positions of the band (as in the paper) is
   natural for the arm, because its decoupling matrix is regular at rest.
5. **Horizon.** With torques of this size the end effector accelerates at
   1 to 5 m/s$^2$; a horizon of 0.2 s ($N = 20$) is short, 0.5 s ($N = 50$ at
   $T = 0.01$ s) is suggested.
6. **$\delta$.** $0.5 M_\Psi$ for the sampling margin (32 % of the band in the
   terminal set) or $0.3 M_\Psi$ (57 %, $s$ halved).
