# Rigour of the reach-avoid certificate and of Lemma 2: Dubins car

Branch `tac-revision`, folder `experiments/strict_xra/`. All numbers below
come from the scripts in this folder and from `python/ra_mpc_constants.py`.
The committed controllers in `controllers/*.py` were not modified.

Notation follows the paper. Output $y = (x_1, x_2) \in \mathbb{R}^2$, state
$x = (x_1, x_2, \theta, v)$, input $u = (\omega, a)$ with
$\lvert \omega \rvert, \lvert a \rvert \le 5$. Safe-set function
$\psi(y) = \alpha(y)\, h(y)$ with

$$
h(y) = -(y_1^4 + y_2^4 - 16)(y_1^4 + y_2^4 - 4), \qquad
\alpha(y) = 10^{-3}\,\big(300 - \phi(y)\big), \qquad
\phi(y) = y_2^2 + \Big(\frac{y_1 + 1.7}{0.5}\Big)^2 - 0.4 ,
$$

safe set $\{\psi > 0\}$ (an annulus), target set $\mathcal{X}_T = \{\phi \le 0\}$.
Area of the safe set in the output plane: $7.41$.

---

## 1. The k1 SOS program and what it certifies

The first backstepping layer treats the output as a single integrator
$\dot y = k_1(y)$. The program in `matlab/solve_vanilla_k1_controller.m` is

$$
\begin{aligned}
\min_{k_1,\,\lambda,\,\delta,\,s_1,\,s_2}\quad & \delta \\
\text{s.t.}\quad & \nabla\psi(y)^\top k_1(y) - \lambda\,\psi(y) + \delta
  - s_1(y)\,\psi(y) - s_2(y)\,\phi(y) \in \Sigma[y], \\
& s_1, s_2 \in \Sigma[y], \qquad \delta \ge 0, \qquad \lambda \ge \xi_0 ,
\end{aligned}
\tag{1}
$$

with $\xi_0 = 10^{-8}$, $\deg k_1 = d_v$, $\deg s_i = d_s$ and $\Sigma[y]$ the
SOS cone. By the S-procedure, (1) implies

$$
\nabla\psi^\top k_1 \;\ge\; \lambda\,\psi - \delta
\qquad \text{on } \{\psi \ge 0\} \cap \{\phi \ge 0\},
\tag{2}
$$

so $V \doteq \psi - \delta/\lambda$ satisfies $\dot V = \nabla\psi^\top k_1 \ge \lambda V$
outside the target. The certified set of this layer is

$$
\{\, V \ge 0 \,\} = \{\, \psi \ge \delta/\lambda \,\}.
\tag{3}
$$

Hence the only quantity that matters for the size of the set is the ratio
$\delta/\lambda$: the smaller, the larger the set. Since
$\max \psi = 10.81$, any $\delta/\lambda$ close to $10.8$ leaves nothing.

The backstepping certificate of Lemma 1 for the Dubins car (relative degree 2
in each output channel, $\dot y = v(\cos\theta, \sin\theta)$) is

$$
\Psi(x) = \psi(y) - \frac{\delta}{\lambda}
  - \frac{1}{2\mu}\,\big\lVert \dot y - k_1(y) \big\rVert^2 ,
\tag{4}
$$

so $\{\Psi \ge 0\}$ requires both $\psi \ge \delta/\lambda$ and a heading and speed
close to the field $k_1$:
$\lVert \dot y - k_1 \rVert^2 \le 2\mu\,(\psi - \delta/\lambda)$.

---

## 2. The committed controller: a certificate at solver tolerance

### 2.1 Homogeneity of (1)

All constraints of (1) are linear in $(k_1, \lambda, \delta, s_1, s_2)$ and the
cones are invariant under positive scaling. If a tuple is feasible, so is
$c\,(k_1, \lambda, \delta, s_1, s_2)$ for every $c > 0$. Therefore

$$
r^\star \doteq \inf \frac{\delta}{\lambda}
\quad\text{is scale invariant, and the optimum of (1) is}\quad
\delta^\star = r^\star \xi_0 ,\ \ \lambda^\star = \xi_0 .
\tag{5}
$$

The solver drives the whole solution to the lower bound $\xi_0 = 10^{-8}$,
which is the order of the interior-point tolerances.

### 2.2 Evidence

Sweep of the lower bound on $\lambda$ (degree $d_v = d_s = 4$; MOSEK status
`pinf = dinf = numerr = 0` in every solve, each solve below 0.1 s):

| lower bound on $\lambda$ | $\lambda$ | $\delta$ | $\delta/\lambda$ | certified area outside $\mathcal{X}_T$ | SOP samples within bounds (of 1000) |
| --- | --- | --- | --- | --- | --- |
| $10^{-8}$ (committed) | $1.09 \times 10^{-8}$ | $1.67 \times 10^{-10}$ | 0.015 | 6.79 | 249 |
| $10^{-4}$ | $1.0 \times 10^{-4}$ | $9.9 \times 10^{-4}$ | 9.87 | 0.88 | 9 |
| $10^{-2}$ | $10^{-2}$ | 0.103 | 10.31 | 0.49 | 4 |
| $0.1$ | 0.1 | 1.08 | 10.76 | 0.03 | 0 |
| $1$ | 1 | 10.80 | 10.80 | 0.001 | 0 |
| $5$ | 5 | 54.0 | 10.80 | 0.001 | 0 |

By (5) the ratio should not depend on the row. It converges to about $10.8$
once the solution leaves the tolerance scale. Raising the degree at the
committed bound makes it worse: degree 6 gives $\delta/\lambda = 3.5 \times 10^{-5}$
and degree 8 returns $\delta = -1.2 \times 10^{-11} < 0$, which violates a
constraint of (1). The committed value $0.015$ is a numerical artefact:
with $\lambda = 10^{-8}$, condition (2) is violated by about
$(r^\star - 0.015)\,\xi_0 \approx 10^{-7}$, which is within tolerance.

Consistent observations on the committed full-state certificate (sampled on
$\{\Psi \ge 0\} \setminus \mathcal{X}_T$, $\theta \in [2\pi/3, 4\pi/3]$, $v \in [0.1, 1]$):

| quantity | value |
| --- | --- |
| share of states with $\dot\Psi < 0$ (continuous time) | 0.85 % |
| median of $\dot\Psi/\Psi$ | $2.2 \times 10^{-5}$ |
| median of $\lvert \dot\Psi \rvert$, of $\lvert \ddot\Psi \rvert$ (frozen input) | $1.2 \times 10^{-4}$, $46$ |
| $u^c$ within $\lvert u \rvert \le 5$ | 63 % (max $\lvert \omega \rvert = 96$) |

### 2.3 Why the continuous controller still works

From 40 states in $\{\Psi \ge 0.1\}$ the continuous loop under $u^c$ reaches the
target 37 times, median 7.4 s, while $\Psi$ only moves from 4.80 to 4.86.
The exponential bound would predict a reaching time of order

$$
t \approx \frac{1}{\lambda}\,\ln\frac{M_\Psi}{\Psi_0} \approx 8 \times 10^{7}\ \text{s}.
$$

The motion is driven by the $k_1$ field; $\Psi$ is almost flat along it and
does not explain the behaviour.

---

## 3. Topological limit of the current $\psi$

### 3.1 A critical point outside the target bounds $\delta/\lambda$

If $\nabla\psi(y^\star) = 0$, $\psi(y^\star) > 0$ and $\phi(y^\star) > 0$, then (2) at
$y^\star$ reads $0 \ge \lambda\psi(y^\star) - \delta$, hence for every $k_1$

$$
\frac{\delta}{\lambda} \;\ge\; \max\big\{\psi(y^\star) : \nabla\psi(y^\star) = 0,\ \phi(y^\star) > 0\big\}.
\tag{6}
$$

Critical points of $\psi$ (Newton from symmetric guesses, Hessian
classification):

| point | $\psi$ | type | location |
| --- | --- | --- | --- |
| $(-1.7782,\, 0)$ | 10.8135 | global maximum | inside $\mathcal{X}_T$ |
| $(1.7743,\, 0)$ | 9.0742 | local maximum | outside $\mathcal{X}_T$ |

So $\delta/\lambda \ge 9.0742$ for any $k_1$ and any degree. The largest set
the current $\psi$ can certify is $\{\psi > 9.0742\}$, a crescent around the
target with area $1.66$ outside $\mathcal{X}_T$ (22 % of the safe set).

### 3.2 Origin of the local maximum

On the ridge $y_1^4 + y_2^4 = 10$ near $(10^{1/4}, 0)$ we have
$y_1 = 10^{1/4} - y_2^4/(4 \cdot 10^{3/4}) + O(y_2^8)$, so

$$
\phi = \phi_0 + y_2^2 + O(y_2^4)
\quad\Longrightarrow\quad
\alpha = \alpha_0 - 10^{-3}\,y_2^2 + O(y_2^4).
$$

The flat side of the superellipse makes the far point a local minimum of the
distance to the target centre, so the weight $\alpha$ and hence $\psi$ have a
local maximum there.

### 3.3 Some second critical point is unavoidable

$\psi$ vanishes on both boundary curves of the annulus and is positive inside.
For such a function the critical points satisfy
$\#\max - \#\text{saddle} + \#\min = \chi(\text{annulus}) = 0$. A single
maximum therefore forces at least one saddle. A redesigned $\psi$ can only
make it harmless by placing every critical point other than the global
maximum inside $\mathcal{X}_T$, which needs a weight that is monotone along the
ring outside the target.

---

## 4. Making the certificate rigorous

### 4.1 Fix the scale and bound the field

With $\lambda$ fixed, $k_1$ is unbounded in (1). For $\lambda \in \{0.1, 0.5, 1\}$
(degree 6) the solver returned $\delta/\lambda \approx 9.065$ but
$\lVert k_1 \rVert \approx 2 \times 10^4$. Then $\{\Psi \ge 0\}$ in (4) needs
$\dot y \approx k_1$, which no state with $v \le 1$ satisfies: $\{\Psi \ge 0\}$ was
empty on $4 \times 10^5$ samples and the SOP stage had no samples.

Bound added to (1) for $i = 1, 2$:

$$
K \mp k_{1,i}(y) - s_{\pm,i}(y)\,\psi(y) \in \Sigma[y], \quad s_{\pm,i} \in \Sigma[y]
\quad\Longrightarrow\quad \lvert k_{1,i} \rvert \le K \ \text{on } \{\psi \ge 0\},
\tag{7}
$$

with $K = 1$ (the speed is at most about 1). The ratio $K/\lambda$ is now the
real design parameter.

### 4.2 Numerical accuracy needs output scaling

Without scaling (degree 6, $K = 1$), the a-posteriori check of (2) on a
$1501 \times 1501$ grid of $\{\psi \ge 0, \phi \ge 0\}$ finds violations of order
$10^{-3}$. The margin is $m(y) \doteq \nabla\psi^\top k_1 - \lambda\psi + \delta$, and
the relative value is $\min m / (\lambda \max\psi)$.

| $\lambda$ | $\delta/\lambda$ | $\min m$ (relative) | $m$ at the far maximum | $\max \lvert k_1 \rvert$ | certified area |
| --- | --- | --- | --- | --- | --- |
| 0.005 | 8.040 | $-5.3 \times 10^{-3}$ ($-9.8 \times 10^{-2}$) | $-5.2 \times 10^{-3}$ | 1.004 | 2.79 |
| 0.01 | 8.755 | $-3.3 \times 10^{-3}$ ($-3.1 \times 10^{-2}$) | $-3.2 \times 10^{-3}$ | 1.005 | 2.13 |
| 0.02 | 9.413 | $-1.5 \times 10^{-3}$ ($-7.1 \times 10^{-3}$) | $+6.8 \times 10^{-3}$ | 1.004 | 1.28 |
| 0.05 | 10.174 | $-4.4 \times 10^{-3}$ ($-8.1 \times 10^{-3}$) | $+5.5 \times 10^{-2}$ | 1.030 | 0.61 |
| 0.1 | 10.508 | $+2.2 \times 10^{-3}$ | $+0.144$ | 1.018 | 0.30 |
| 0.2 | 10.658 | $+2.9 \times 10^{-3}$ | $+0.317$ | 1.022 | 0.14 |

The first two rows "beat" bound (6) only by violating (2) at the far maximum.
The cause is conditioning: $y \in [-2.2, 2.2]$ and the products $s_1\psi$ reach
degree 16, with $2.2^{16} \approx 3 \times 10^5$.

Scaling. Solve (1) and (7) in $\tilde y = y/\sigma$ with

$$
\tilde\psi(\tilde y) = \psi(\sigma\tilde y), \quad
\tilde\phi(\tilde y) = \phi(\sigma\tilde y), \quad
\tilde k_1 = k_1/\sigma, \quad \tilde K = K/\sigma ,
$$

then map back with $k_1(y) = \sigma\,\tilde k_1(y/\sigma)$. Condition (2) is
unchanged because
$\nabla_{\tilde y}\tilde\psi^\top \tilde k_1 = \sigma\,\nabla\psi^\top (k_1/\sigma)$,
and $\lambda$, $\delta$ are the same. Results with $\sigma = 2$, $K = 1$:

| degree | $\lambda$ | $\delta/\lambda$ | $\min m$ (relative) | $m$ at the far maximum | $\max \lvert k_1 \rvert$ | certified area |
| --- | --- | --- | --- | --- | --- | --- |
| 6 | 0.01 | 9.128 | $-2.8 \times 10^{-4}$ ($-2.6 \times 10^{-3}$) | $+5.4 \times 10^{-4}$ | 0.998 | 1.56 |
| 6 | 0.02 | 9.919 | $+1.0 \times 10^{-4}$ | $+1.7 \times 10^{-2}$ | 0.999 | 0.83 |
| 6 | 0.05 | 10.445 | $+1.2 \times 10^{-4}$ | $+6.9 \times 10^{-2}$ | 0.999 | 0.36 |
| 6 | 0.1 | 10.620 | $+8.6 \times 10^{-5}$ | $+0.155$ | 0.999 | 0.18 |
| 8 | 0.01 | 9.0731 | $-3.9 \times 10^{-6}$ ($-3.6 \times 10^{-5}$) | $-1.1 \times 10^{-5}$ | 0.989 | 1.63 |
| **8** | **0.02** | **9.0743** | $+3.0 \times 10^{-5}$ | $+1.2 \times 10^{-6}$ | **0.944** | **1.63** |
| 8 | 0.05 | 9.966 | $+3.5 \times 10^{-4}$ | $+4.5 \times 10^{-2}$ | 0.999 | 0.79 |
| 8 | 0.1 | 10.382 | $+3.5 \times 10^{-4}$ | $+0.131$ | 0.999 | 0.42 |

The bold row attains the topological limit (6), satisfies (2) with a positive
margin on the whole grid and respects (7). It is the rigorous choice:
$\sigma = 2$, $K = 1$, $\lambda = 0.02$, $d_v = d_s = 8$. Solve time below 0.6 s.

Trade-off visible in the table: for fixed $K$ a larger $\lambda$ needs a faster
field than (7) allows near the far side, so $\delta/\lambda$ grows and the set
shrinks. The limit (6) is reached only for $\lambda \lesssim 0.02$ at $K = 1$.

---

## 5. Full-state certificate with the rigorous $k_1$

The full synthesis (backstepping + bounded-input SOP, 3000 samples) was run
with the bold setting for $\mu \in \{0.1, 1\}$. "vanilla $k_1$" uses the
rigorous $k_1$ of Section 4; "SOP $k_1$" is re-solved by
`solve_k1_controller_sop` with sampled input constraints (still unscaled).
Checks on $4 \times 10^5$ samples of the synthesis box
($\theta \in [2\pi/3, 4\pi/3]$, $v \in [0.1, 1]$):

| quantity | $\mu = 0.1$, vanilla $k_1$ | $\mu = 0.1$, SOP $k_1$ | $\mu = 1$, vanilla $k_1$ | $\mu = 1$, SOP $k_1$ |
| --- | --- | --- | --- | --- |
| $\{\Psi \ge 0\}$ share of the box | 2.6 % | 3.3 % | 8.7 % | 12.1 % |
| MPC initial states inside $\{\Psi \ge 0\}$ (of 300) | 11 | 13 | 54 | 83 |
| $\min(\dot\Psi - \lambda\Psi)$ on $\{\Psi \ge 0\} \setminus \mathcal{X}_T$ | $+1.9 \times 10^{-3}$ | $-4.6 \times 10^{-3}$ | $+6.2 \times 10^{-4}$ | $-9.4 \times 10^{-3}$ |
| share violating $\dot\Psi \ge \lambda\Psi$ | 0 | 0.2 % | 0 | 0.6 % |
| median $\dot\Psi/\Psi$ | 0.85 | 1.89 | 0.78 | 1.33 |
| $u^c$ within $\lvert u \rvert \le 5$ on $\{\Psi \ge 0\}$ | 95.4 % | 93.9 % | 16.2 % | 16.3 % |
| $\max \lvert \omega \rvert$, $\max \lvert a \rvert$ | 18.8, 2.9 | 23.4, 3.1 | 278, 35 | 329, 43 |
| continuous loop from 40 states: reached / $\Psi$ nondecreasing | 40 / 40 | 40 / 40 | 40 / 40 | 40 / 37 |
| median (max) reaching time | 3.1 s (5.8 s) | 3.2 s (6.8 s) | 2.5 s (8.1 s) | 2.7 s (12.9 s) |

Observations:

- With the vanilla $k_1$ the full-state ECGBF condition holds on every sample,
  and $\Psi$ now grows at a meaningful rate (median $\dot\Psi/\Psi$ of about
  0.8 instead of $2 \times 10^{-5}$). The certificate explains the motion.
- The SOP stage reintroduces small violations because it is solved unscaled.
  It needs the same scaling (and bound (7)) to be rigorous.
- A larger $\mu$ enlarges $\{\Psi \ge 0\}$ through the velocity term of (4), but
  raises the gains: with $\mu = 1$ the input bounds hold on only 16 % of the set.

---

## 6. Sampled data: Lemma 2 versus direct ZOH simulation

### 6.1 Lemma 2 constants

With $\Delta_k$ the drift bound (18) of the paper,

$$
\gamma(T) = \Big(\sum_j M_{g_j}\,\ell_{u^c_j}\Big)\frac{e^{\ell_{sys}T} - 1}{\ell_{sys}},
\qquad \beta(T) = \gamma(T)\,M_v ,
$$

and the certified period for the terminal margin $\delta$ is the solution of
$\beta(T) = \lambda\delta$:

$$
T^\star(\delta) = \frac{1}{\ell_{sys}}
\ln\!\Big(1 + \frac{\ell_{sys}\,\lambda\,\delta}{\big(\sum_j M_{g_j}\ell_{u^c_j}\big) M_v}\Big)
\approx \frac{\lambda\,\delta}{\big(\sum_j M_{g_j}\ell_{u^c_j}\big) M_v}.
\tag{8}
$$

Sampled constants (inflated by 10 %) on $\{\Psi \ge 0\}$ of the synthesis box,
$\ell_f = 1.1$, $\ell_{g_j} = 0$:

| certificate | $\lambda$ | $M_{g_1}$, $M_{g_2}$ | $\ell_{u^c_1}$, $\ell_{u^c_2}$ | $M_v$ | $T^\star(0.1 \max\Psi)$ |
| --- | --- | --- | --- | --- | --- |
| committed | $1.09 \times 10^{-8}$ | $4.9 \times 10^{-4}$, 11.0 | 1246, 85 | 106 | $1.1 \times 10^{-14}$ s |
| rigorous, $\mu = 0.1$, vanilla $k_1$ | 0.02 | 5.3, 6.2 | 299, 27 | 17.1 | $1.1 \times 10^{-7}$ s |
| rigorous, $\mu = 1$, vanilla $k_1$ | 0.02 | 1.0, 1.9 | 5057, 481 | 306 | $1.9 \times 10^{-9}$ s |

The rigorous certificate gains seven orders of magnitude, but (8) is still
far from practical because of the global Lipschitz constant $\ell_{u^c_1}$ and the
small $\lambda$.

### 6.2 Direct ZOH test

One ZOH step of length $T$ with the frozen input $u^c(x)$ from states in
$\{\Psi \ge 0\} \setminus \mathcal{X}_T$, integrated with 50 RK4 sub-steps. To
second order

$$
\Psi(x^+) - \Psi(x) = \dot\Psi(x)\,T + \tfrac{1}{2}\ddot\Psi(x)\,T^2 + O(T^3),
\qquad
T_c \doteq \frac{2\,\dot\Psi}{\lvert \ddot\Psi \rvert},
\tag{9}
$$

so the increment is guaranteed positive only for $T \ll T_c$. For the
committed certificate the median $T_c$ is $5 \times 10^{-6}$ s.

Share of states where $\Psi$ decreases over one step:

| $T$ [s] | 0.05 | 0.02 | 0.01 | 0.005 | 0.002 | 0.001 | $10^{-4}$ | $10^{-6}$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| committed | 44 % | | 43 % | | | 43 % | 41 % | 13 % |
| rigorous, $\mu = 0.1$ | 41 % | 21 % | 9.5 % | 2.7 % | 0.4 % | 0.1 % | | |

The committed entry at $T = 10^{-6}$ s is the second-order estimate (9); all other entries are simulated.
For the rigorous certificate the smallest one-step change is
$-0.16$, $-4.4 \times 10^{-3}$, $-1.2 \times 10^{-4}$ at $T = 0.05$, $0.01$, $0.002$ s.
Restricting to states where $u^c$ meets the bounds changes the shares by at most
1.5 points.

---

## 7. Conclusions

1. Re-solving with a different bound on $\lambda$ does not make the SOS harder.
   Every solve converged with zero infeasibility and numerical error flags,
   in under a second.
2. The reach-avoid set reported with the committed controller is not
   certified: the solution sits at the solver tolerance (Section 2).
3. A rigorous certificate is obtained with output scaling, a fixed $\lambda$ and
   the field bound (7): $\sigma = 2$, $K = 1$, $\lambda = 0.02$, degree 8.
   It satisfies the ECGBF condition with positive margin in both the output
   and the full state, and the continuous loop reaches the target from every
   tested state with $\Psi$ nondecreasing.
4. Its size is limited by the far-side local maximum of $\psi$ (bound (6)):
   22 % of the safe set in the output plane, 2.6 % of the state box at
   $\mu = 0.1$, 11 of the 300 MPC initial states.
5. Sampled data: Lemma 2 still certifies only $T^\star \approx 10^{-7}$ s. In
   simulation, monotone growth of $\Psi$ per step holds from about $T = 2$ ms on.
   At $T = 0.05$ s, $\Psi$ still decreases from 41 % of the states.

## 8. Open points

- Scale the SOP stage (`solve_k1_controller_sop`) and add (7) there, so that
  the input-bounded certificate is also rigorous; then check the 4 % of
  $\{\Psi \ge 0\}$ where $u^c$ exceeds the bounds.
- Decide whether to redesign $\psi$ so that its non-global critical points lie
  inside $\mathcal{X}_T$ (Section 3.3). This is the only way to certify more than
  the crescent.
- Sampling period: either use $T \approx 2$ to 5 ms in the experiments (MPC
  horizon $N T$ of about 1 s then needs $N$ of 200 to 500), or keep $T = 0.05$ s
  and state that Lemma 2 is an existence result whose bound is not met.
- Tighten (8): the global constant $\ell_{u^c}$ dominates. A bound using
  $\lvert \ddot\Psi \rvert$ along ZOH arcs, as in (9), would be far less conservative.
- Re-run the MPC experiments with the rigorous certificate. The smaller
  $\mathcal{X}_{RA}$ makes the enlargement by MPC more visible.

## 9. Reproduction

```
cd experiments/strict_xra
matlab -batch "run('run_lmin.m')"    # Section 2.2 (lambda sweep, degree 4)
matlab -batch "run('run_deg.m')"     # degree test
matlab -batch "run('run_kmax.m')"    # Section 4.2, unscaled, K = 1
matlab -batch "run('run_scaled.m')"  # Section 4.2, scaled
matlab -batch "run('run_full.m')"    # Section 5, full synthesis
# exports land in controllers/tmp_*.py: move them to out/ before anything else
LAM=0.02 python check_full.py out/rigorous_mu0.1_vanilla_k1.py
```
