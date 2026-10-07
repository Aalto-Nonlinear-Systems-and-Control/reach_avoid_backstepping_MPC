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

---

## 9. Redesign of $\psi$: a steeper weight

The safe set is the zero-superlevel set of $h$ alone, so any weight
$\alpha > 0$ on the safe set gives the same safe set for $\psi = \alpha h$. The
committed weight is almost constant ($\alpha$ at the target over $\alpha$ at the
far side is 1.2), which makes the far critical value high (Section 3) and
$\psi$ flat. A weight that grows monotonically along the ring towards the
target turns the far critical point into a saddle with a small value.

### 9.1 Critical points for candidate weights

By symmetry in $y_2$, the critical points of $\psi = \alpha(y_1)\,h$ lie on the
axis $y_2 = 0$: off the axis $\partial\psi/\partial y_2 = \alpha\,\partial h/\partial y_2$
vanishes only on the ridge $y_1^4 + y_2^4 = 10$, where
$\nabla\psi = h\,\nabla\alpha \ne 0$. On the axis we solve
$\partial\psi/\partial y_1 = 0$ and classify with the second derivatives (the
fourth derivative in $y_2$ when the second vanishes).

| $\alpha(y)$ | $\min\alpha$ on safe set | far critical point | its $\psi$ / $\max\psi$ | largest certifiable area outside $\mathcal{X}_T$ |
| --- | --- | --- | --- | --- |
| $10^{-3}(300 - \phi)$ (committed) | 0.245 | local maximum at $(1.774, 0)$ | 83.9 % | 24 % |
| $(2.5 - y_1)^2$ | 0.252 | saddle at $(1.679, 0)$ | 3.2 % | 74 % |
| $(2.2 - y_1)^3$ | 0.008 | saddle at $(1.593, 0)$ | 0.23 % | 78 % |
| $(2.2 - y_1)^4$ | 0.000 | saddle at $(1.562, 0)$ | 0.04 % | 79 % |

$(2.5 - y_1)^2$ keeps the degree of $\psi$ at 10 and stays well away from zero,
so it is the choice tested below, normalised as
$\psi = (2.5 - y_1)^2\,h/600$ ($\max\psi = 1.10$, saddle value $0.0357$).
The remaining quarter of the ring around the far saddle cannot be certified by
any weight (Section 3.3).

### 9.2 Rigorous $k_1$ (scaled, $K = 1$, degree 8)

| $\lambda$ | $\delta/\lambda$ | $\min m$ (relative) | $\max \lvert k_1 \rvert$ | certified area outside $\mathcal{X}_T$ |
| --- | --- | --- | --- | --- |
| 0.0005 | 0.0318 | $-1.0 \times 10^{-6}$ ($-1.8 \times 10^{-3}$), poor solve (feasratio 0.45) | 0.770 | 79 % |
| 0.001 | 0.0355 | $+5.4 \times 10^{-7}$ | 0.760 | 75 % |
| 0.002 | 0.0353 | $+1.1 \times 10^{-7}$ | 0.827 | 75 % |
| **0.005** | **0.0356** | $+2.9 \times 10^{-7}$ | **0.764** | **75 %** |
| 0.01 | 0.0355 | $-1.3 \times 10^{-6}$ ($-1.2 \times 10^{-4}$) | 0.823 | 75 % |
| 0.02 | 0.0958 | $+9.9 \times 10^{-5}$ | 1.000 | 56 % |
| 0.05 | 0.698 | $+9.6 \times 10^{-5}$ | 1.000 | 14 % |
| 0.1 | 0.899 | $+9.8 \times 10^{-5}$ | 1.000 | 7 % |

For $\lambda \le 0.01$ the solution reaches the saddle bound $0.0357$. Above that
the field bound (7) becomes active and the set shrinks. Note that $\psi$ is
normalised differently from Section 4 ($\max\psi = 1.1$ instead of 10.8), so
the values of $\lambda$ are not directly comparable; the relevant rate is
$\dot\Psi/\Psi$ below.

### 9.3 Full-state certificate ($\lambda = 0.005$, rigorous vanilla $k_1$)

| quantity | committed | Section 5, $\mu = 0.1$ | $\alpha = (2.5 - y_1)^2$, $\mu = 0.1$ | $\alpha = (2.5 - y_1)^2$, $\mu = 0.5$ |
| --- | --- | --- | --- | --- |
| violations of $\dot\Psi \ge \lambda\Psi$ | 0.85 % | 0 | 0 (min margin $+6.4 \times 10^{-3}$) | 0 (min margin $+6.0 \times 10^{-3}$) |
| median $\dot\Psi/\Psi$ | $2 \times 10^{-5}$ | 0.85 | 1.26 | 1.09 |
| $\{\Psi \ge 0\}$ share of the box | (not certified) | 2.6 % | 4.1 % | 15.9 % |
| MPC initial states inside (of 300) | (not certified) | 11 | 17 | 79 |
| $u^c$ within $\lvert u \rvert \le 5$ | 63 % | 95.4 % | 99.8 % (max $\lvert\omega\rvert$ 7.1) | 91.7 % (max $\lvert\omega\rvert$ 46) |
| continuous loop, 40 states | 37 reached | 40, monotone | 40, monotone, median 2.7 s | 40, monotone, median 3.2 s |
| ZOH, $\Psi$ decreases, $T = 0.05$ s | 44 % | 41 % | **0.9 %** | 10.4 % |
| ZOH, $\Psi$ decreases, $T = 0.01$ s | 43 % | 9.5 % | 0.1 % | 1.2 % |
| Lemma 2, $T^\star(0.1 \max\Psi)$ | $10^{-14}$ s | $1.1 \times 10^{-7}$ s | $2.0 \times 10^{-7}$ s | $1.1 \times 10^{-8}$ s |

For $\alpha = (2.5 - y_1)^2$, $\mu = 0.1$, the remaining decreases at
$T = 0.05$ s are spread over all levels of $\Psi$ (0.98 % for $\Psi \ge 0.01$,
0.99 % for $\Psi \ge 0.1$), are at most $7.5 \times 10^{-3}$ in size, and every
next state stays in $\{\Psi \ge 0\}$. At $T = 0.1$ and $0.2$ s the shares are
2.1 % and 6.1 %.

### 9.4 What the redesign fixes and what it does not

Fixed or largely improved:

- size of the certified set in the output plane: 24 % to 75 % of the ring;
- flatness of $\Psi$ and hence the per-step behaviour under ZOH: at the paper's
  $T = 0.05$ s, $\Psi$ decreases from 0.9 % of the states instead of 41 to 44 %;
- input bounds: $u^c$ within the bounds on 99.8 % of $\{\Psi \ge 0\}$ for $\mu = 0.1$.

Not fixed:

- Lemma 2's certified $T^\star$ remains of order $10^{-7}$ s. The global
  Lipschitz constant $\ell_{u^c}$ (91 here) and the factor $M_v$ dominate (8).
  This is a property of the bound, not of $\Psi$; a sharper lemma (for example
  bounding $\ddot\Psi$ along ZOH arcs, as in (9)) is needed for a practical $T$.
- The full-state set is still small for $\mu = 0.1$ (17 of 300 initial states)
  because of the velocity term in (4). $\mu = 0.5$ quadruples it but weakens
  the input bounds and the ZOH behaviour. MPC is what enlarges the region.
- The SOP stage is still unscaled and reintroduces tiny violations
  ($2 \times 10^{-4}$ of the samples for $\mu = 0.5$).

Scripts: `run_alpha.m` (Section 9.2) and `run_full_alpha.m` (Section 9.3);
exports in `out/rigorous_alpha2_mu0.1_vanilla_k1.py` and
`out/rigorous_alpha2_mu0.5_vanilla_k1.py`.

---

## 10. Where Lemma 2 loses (certificate of Section 9, $\mu = 0.1$, $\lambda = 0.005$)

### 10.1 The domain of the constants is not the cause

All constants were already taken on $\{\Psi \ge 0\}$ (minus $\mathcal{X}_T$ where
relevant), intersected with the sampling box, not on the whole safe set or
state space. The box limits $\theta \in [2\pi/3, 4\pi/3]$ and $v \in [0.1, 1]$, so
it is if anything smaller than $\mathcal{X}_{RA}$. Changing the domain moves
$T^\star$ by less than a factor 2:

| domain | $M_{g_1}$, $M_{g_2}$ | $\ell_{u^c_1}$, $\ell_{u^c_2}$ | $M_v$ | $T^\star(0.1 \max\Psi)$ |
| --- | --- | --- | --- | --- |
| $\{\Psi \ge 0\}$, all headings, $v \in [0.1, 1.5]$ | 3.03, 4.58 | 110, 8.9 | 9.4 | $1.5 \times 10^{-7}$ s |
| $\{\Psi \ge 0\}$, sampling box (used so far) | 3.03, 4.59 | 90, 9.1 | 7.7 | $2.1 \times 10^{-7}$ s |
| $\{\Psi \ge 0.1 \max\Psi\}$, sampling box | 2.83, 4.34 | 90, 8.7 | 7.1 | $2.5 \times 10^{-7}$ s |

### 10.2 Factor-by-factor comparison along real ZOH arcs

From 13 500 states in $\mathcal{X}_{RA} \setminus \mathcal{X}_T$, arcs of length
$T = 0.05$ s with the input frozen at $u^c(x_k)$:

| step of the proof | bound | actual (median / max) | overestimate |
| --- | --- | --- | --- |
| drift $\Delta_k(t)$, eq. (18) | $\frac{e^{\ell_{sys}T} - 1}{\ell_{sys}}\lVert \dot x_k \rVert$ | | 1.02 (tight) |
| $\lvert u^c(x(t)) - u^c(x_k) \rvert / \Delta_k(t)$ | $\ell_{u^c} = 99$ | 0.67 / 15 | 6.5 to 150 |
| mismatch $e(t) = \sum_j L_{g_j}\Psi\,(u^c_j(x(t)) - u^c_j(x_k))$ | $\beta_k(T)$ | 0.016 / 0.78 | median 730 |
| same, state-independent | $\beta(T) = 135$ | 0.016 / 0.78 | median 8 400 |
| available margin | $\lambda\Psi_k$, median $10^{-3}$ | true $\dot\Psi$: median $\dot\Psi/\Psi = 1.26$ | about 250 |

Two factors dominate: the product of separate suprema $M_g\,\ell_{u^c}\,M_v$
(three to four orders of magnitude), and the use of the uniform worst-case
rate $\lambda$ as the only margin (more than two orders). Even the true mismatch
exceeds $\lambda\Psi_k$ in 94 % of the states, yet $\Psi$ increases in 99 % of
the ZOH steps, because the actual $\dot\Psi$ is much larger than $\lambda\Psi$.

### 10.3 A second-order alternative

With a frozen input, $\dot\Psi$ at $t_k$ equals its continuous-time value, so

$$
\Psi(x(t_{k+1})) \;\ge\; \Psi_k + \dot\Psi(x_k)\,T - \tfrac{1}{2} M_2\,T^2,
\qquad
M_2 \doteq \sup \lvert \ddot\Psi \rvert \ \text{along ZOH arcs in } \mathcal{X}_{RA},
\tag{10}
$$

with $\ddot\Psi = (f + g u)^\top \nabla^2\Psi\,(f + g u) + \nabla\Psi^\top \tfrac{\partial f}{\partial x}(f + g u)$
for constant $u$. This bound needs no Lipschitz constant of $u^c$. Increase over
one step holds if $\dot\Psi(x_k) > \tfrac{1}{2} M_2 T$, i.e.

$$
T^\star = \frac{2\,\eta}{M_2}, \qquad
\eta = \lambda\delta \ \ \text{(certified rate)} \quad\text{or}\quad
\eta = \min_{\{\Psi \ge \delta\} \setminus \mathcal{X}_T} \dot\Psi \ \ \text{(additive margin)}.
\tag{11}
$$

| quantity | value |
| --- | --- |
| $M_2$ over $\mathcal{X}_{RA} \times$ (box of $u^c$ values) | 576 (too loose: pairs a state with inputs it never sees) |
| $M_2$ along arcs, $T = 0.01$ to $0.1$ s | 19 to 23 (99 % quantile about 7) |
| $\eta = \lambda\delta$, $\delta = 0.1 \max\Psi$ | $5.2 \times 10^{-4}$ |
| $\eta = \min \dot\Psi$ on $\{\Psi \ge \delta\} \setminus \mathcal{X}_T$ | 0.0072 to 0.0083 |

| bound | $T^\star$ |
| --- | --- |
| Lemma 2 as in the paper | $2.1 \times 10^{-7}$ s |
| (11) with $\eta = \lambda\delta$ | $5.4 \times 10^{-5}$ s |
| (11) with the additive margin | $7.5 \times 10^{-4}$ s |
| simulation: $\Psi$ increases in 99 % of steps | $0.05$ s |

The additive margin can be certified by SOS (require
$\dot\Psi \ge \lambda\Psi + \eta$ on $\{\Psi \ge \delta\} \setminus \mathcal{X}_T$), and
$M_2$ by an SOS or interval bound over $\mathcal{X}_{RA}$ and the drift
neighbourhood of each state; the numbers above are sampled. The remaining gap
to 0.05 s comes from $\eta$ (the slowest region, near the far saddle) and from
$M_2$ being a supremum; for recursive feasibility only states near the level
$\delta$ need an increase, which would allow $\eta$ to be taken on a band around
$\{\Psi = \delta\}$ only.

Script: `lemma2_tightness.py`.

### 10.4 Uniform period versus terminal level (certificate of Section 9)

Largest uniform period for which $\Psi$ increases from every sampled state of
$\{\Psi \ge \delta\} \setminus \mathcal{X}_T$ (script `t_vs_delta.py`, box
sampling): 2 ms for $\delta \le 0.2 \max\Psi$, 5 ms at $0.3$, 10 ms at $0.5$,
100 ms at $0.7$. The slowest states sit at the edge of the target,
$(-1.78, \pm 0.56)$, where the ring is flat and $\psi$ of Section 9 barely
changes along it. This motivates the cross term of Section 11.

---

## 11. Parameter design for the Dubins car

Folder `design/`. Goal: one consistent set ($\psi$, $\lambda$, $\mu$, $\delta$, $T$)
with a rigorous certificate, inputs within bounds, and $T = 0.05$ s covered by
the second-order bound (11).

### 11.1 Sampling without a state box

The earlier analysis sampled a box in $(\theta, v)$ chosen by hand. It is not
needed: by (4), in the coordinates $(y, \dot y)$ with
$\dot y = v(\cos\theta, \sin\theta)$,

$$
\{\Psi \ge \delta\} = \Big\{ (y, \dot y) :
\lVert \dot y - k_1(y) \rVert^2 \le 2\mu\big(\psi(y) - \tfrac{\delta_s}{\lambda} - \delta\big) \Big\},
\tag{12}
$$

a velocity disc of radius $r(y) = \sqrt{2\mu(\psi - \delta_s/\lambda - \delta)}$
around the field $k_1(y)$. `design/cert.py` samples $y$ in the ring and
$\dot y$ in the disc (half of the samples on its boundary, where $\Psi = \delta$),
so all constants are taken on the set itself, for every heading. The speed
singularity $v = 0$ of the decoupling matrix is inside the set exactly when
$\lVert k_1(y) \rVert \le r(y)$ somewhere; the smallest speed
$\min_y(\lVert k_1 \rVert - r)$ is reported below. `cert.py` builds $\Psi$ and
$u^c$ directly from $k_1$ and agrees with the MATLAB export to $10^{-14}$
(`design/validate.py`).

### 11.2 Weight with a cross term

For a weight depending on $y_1$ only, $\lvert\nabla\psi\rvert \ge 0.02$ is the best
lower bound on any level set outside the target (flat near side, Section
10.4), and $\dot\Psi \le \sqrt{2}\,K\,\lvert\nabla\psi\rvert$ caps the margin $\eta$.
A term $c\,y_1 y_2^2\,(\cdot)$ has positive $y_2$-curvature on the far side
(keeps the saddle) and negative curvature on the target side (pulls towards
$y_2 = 0$). Screening (`design/screen_weights.py`, $\psi$ normalised to
$\max\psi = 1$; "share" is the part of the ring outside $\mathcal{X}_T$ with
$\psi \ge c$):

| weight $\alpha$ | critical points outside $\mathcal{X}_T$ | $\min\lvert\nabla\psi\rvert$ on $\{\psi \ge 0.1\}$ (share) | on $\{\psi \ge 0.3\}$ (share) |
| --- | --- | --- | --- |
| $(2.5 - y_1)^2$ | far saddle, 0.032 | 0.020 (54 %) | 0.020 (32 %) |
| $(2.5 - y_1)^2 + 0.5\,y_1 y_2^2$ | far saddle, 0.032 | 0.079 (63 %) | 0.079 (31 %) |
| $(2.5 - y_1)^2 + 0.2\,y_1 y_2^2 (2.5 - y_1)$ | far saddle, 0.032 | 0.119 (58 %) | 0.119 (31 %) |
| $\mathbf{(2.5 - y_1)^2 + 0.3\,y_1 y_2^2 (2.5 - y_1)}$ | far saddle, 0.032 | **0.131 (59 %)** | **0.141 (30 %)** |
| $(2.5 - y_1)^2 + 0.5\,y_1 y_2^2 (2.5 - y_1)$ | far saddle, 0.032 | 0.034 (62 %) | 0.034 (27 %) |
| $(2.5 - y_1)^2 + 0.7\,y_1 y_2^2 (2.5 - y_1)$ | two new maxima, 0.344 | 0.002 | 0.002 |

Chosen:

$$
\psi(y) = \frac{\alpha(y)\,h(y)}{661.409}, \qquad
\alpha(y) = (2.5 - y_1)^2 + 0.3\,y_1 y_2^2\,(2.5 - y_1),
\tag{13}
$$

of degree 12, with $\alpha \ge 0.25$ on the safe set.

### 11.3 $\lambda$, the field speed, and time scaling

The SOS minimises $\delta_s$ only, so the field it returns is as slow as the
constraints allow. With (13), scaled coordinates, $K = 1$, degree 8:

| $\lambda$ | $\delta_s/\lambda$ | $\max\lvert k_{1,i}\rvert$ | $\min m$ on the grid | certified ring area |
| --- | --- | --- | --- | --- |
| 0.005 | 0.0321 | 0.099 | $-1.3 \times 10^{-6}$ | 82 % |
| 0.02 | 0.0324 | 0.233 | $-5.6 \times 10^{-7}$ | 82 % |
| 0.05 | 0.0333 | 0.523 | $+4.7 \times 10^{-5}$ | 81 % |
| **0.06** | **0.0410** | **0.626** | $+4.4 \times 10^{-4}$ | **76 %** |
| 0.08 | 0.0484 | 0.804 | $+9.0 \times 10^{-4}$ | 73 % |
| 0.1 | 0.1428 | 0.927 | $+1.0 \times 10^{-3}$ | 51 % |

(The $y_1$-only weight loses the saddle bound already at $\lambda = 0.02$.) A
slow field with a fixed $\mu$ puts $v = 0$ inside the velocity discs:
for $\lambda = 0.005$, $\mu = 0.1$ the smallest speed is $-0.37$ and
$\lvert\omega\rvert$ reaches 35. $\lambda = 0.06$ is the base solution used below.

Time scaling. If $(k_1, \lambda, \delta_s)$ satisfies (2), so does
$s\,(k_1, \lambda, \delta_s)$ for $s > 0$. Together with $\mu \to s^2\mu$,

$$
\Psi_s(y, \dot y) = \Psi(y, \dot y / s), \qquad
\text{closed loop } y_s(t) = y(s t), \qquad
T^\star_s = T^\star / s .
\tag{14}
$$

Hence the certified period depends on the product of $T$ and the field
speed, and a slower terminal controller buys a larger certified $T$. Only
the terminal controller is slowed; the MPC closed loop is not.

### 11.4 Scan over $\mu$ and $\delta$ at $s = 1$ ($\lambda = 0.06$)

$\delta$ in units of $M_\Psi = 1 - \delta_s/\lambda = 0.959$. "share": 4D volume of
$\{\Psi \ge \delta\} \setminus \mathcal{X}_T$ relative to $\delta = 0$. $T^\star = 2\eta/M_2$
with sampled $\eta = \min\dot\Psi$ and $M_2$ along arcs. "sim": largest $T$ of the
grid for which $\Psi$ increases from every sampled state.

| $\mu$ | $\delta/M_\Psi$ | share | min speed | $\max\lvert\omega\rvert$ | $\eta$ | $M_2$ | $T^\star$ | sim |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.02 | 0.2 | 46 % | 0.16 | 1.0 | 0.0134 | 3.2 | 8.4 ms | 50 ms |
| 0.02 | 0.3 | 30 % | 0.17 | 0.8 | 0.0193 | 2.8 | 14.0 ms | 50 ms |
| 0.02 | 0.4 | 19 % | 0.19 | 0.8 | 0.0254 | 2.3 | 21.8 ms | 100 ms |
| 0.02 | 0.5 | 11 % | 0.20 | 0.7 | 0.0307 | 1.9 | 32.6 ms | 100 ms |
| 0.05 | 0.2 | 46 % | 0.06 | 2.8 | 0.0134 | 4.6 | 5.8 ms | 10 ms |
| 0.05 | 0.3 | 30 % | 0.08 | 2.2 | 0.0193 | 3.9 | 9.9 ms | 25 ms |
| 0.05 | 0.5 | 11 % | 0.13 | 1.6 | 0.0307 | 2.2 | 27.6 ms | 50 ms |

Smaller $\mu$ is better on every period metric (thinner discs, less velocity
mismatch); it only thins the set in the velocity directions. On the whole
certified set $\{\Psi \ge 0\}$ ($\mu = 0.02$): smallest speed 0.04, inputs within
the bounds on 100 % of the samples, $\max\lvert\omega\rvert = 1.3$,
$\max\lvert a\rvert = 0.25$. The bounded-input SOP stage is therefore not needed
for this example.

### 11.5 Certified settings and MPC closed loop

$s$ is chosen as $0.95\,T^\star(\delta)/T$ with $T = 0.05$ s, so that (11) covers
the sampling period. Verification on $4 \times 10^5$ samples, two seeds:

| setting | $s$ | $\lambda$ | $\mu$ | $\delta$ | $T^\star = 2\eta/M_2$ | speeds on the terminal set | $\max\lvert\omega\rvert$, $\max\lvert a\rvert$ | min one-step increase at 0.05 s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | 0.266 | 0.0160 | $1.42 \times 10^{-3}$ | $0.3 M_\Psi = 0.288$ | 52 ms | $\ge 0.046$ | 0.22, 0.016 | $+2.1 \times 10^{-4}$ |
| R2 | 0.159 | 0.0096 | $5.06 \times 10^{-4}$ | $0.2 M_\Psi = 0.192$ | 52 ms | $\ge 0.026$ | 0.16, 0.006 | $+8.6 \times 10^{-5}$ |

In both, $\Psi$ still increases from every sampled state at $T = 0.1$ s.

MPC (Algorithm 1, `design/ra_mpc_cert.py`): stage cost as in the paper,
$V_f = \omega(M_\Psi - \Psi)$ with $\omega = 1.2\,\sup\ell/\varepsilon_T$
($\varepsilon_T$ the sampled minimal one-step increase), safety back-off
$10^{-3}$, $T = 0.05$ s, the 300 initial states of the paper, continuous plant
(20 s). "Vanilla" is first-step feasibility of the terminal-target MPC.

| | $N = 25$ | $N = 40$ |
| --- | --- | --- |
| vanilla MPC | 245 | 288 |
| RA-MPC, R1 | 252 | 292 |
| RA-MPC, R2 | **280** | 292 |

| | infeasible after $k = 0$ | steps violating $J^*_{k+1} \le J^*_k - \underline{\ell}$ | $t_{\text{entry}}$ median (max) | $t_{\text{bound}}$ median | $\omega$ |
| --- | --- | --- | --- | --- | --- |
| R1, $N = 25$ | 0 | 0 / 4950 | 1.00 s (1.65 s) | 571 s | $1.7 \times 10^5$ |
| R1, $N = 40$ | 0 | 0 / 5927 | 1.10 s (2.00 s) | 43 s | $1.7 \times 10^5$ |
| R2, $N = 25$ | 0 | 0 / 5867 | 1.15 s (1.95 s) | 6030 s | $5.9 \times 10^5$ |
| R2, $N = 40$ | 0 | 0 / 6160 | 1.10 s (2.30 s) | 46 s | $5.9 \times 10^5$ |

All failures are infeasibility at $k = 0$ (terminal set not reachable within
$N$ steps); none occurs later. Nominal and continuous plants give the same
counts (R1, $N = 25$: 252 both).

Two implementation points were needed for this:

- The MPC model uses 4 RK4 sub-steps per sampling interval. With one RK4 step
  the prediction differs from the continuous plant by about $10^{-5}$, which
  is enough to break the terminal constraint of the candidate (29 to 42
  mid-run infeasibilities) because $\partial\Psi/\partial\dot y \sim 1/\mu$ is large.
- The cost is evaluated on the rolled-out trajectory of the returned inputs,
  not on the solver iterate, and the constraint tolerance is $10^{-8}$. With
  $\omega \sim 10^5$ a residual of $10^{-4}$ in the dynamics constraints shows
  up as spurious violations of the cost decrease.

### 11.6 Recommended parameters

- $\psi$ as in (13); $k_1$ from the scaled SOS ($\sigma = 2$, degree 8,
  $\lvert k_{1,i}\rvert \le K$), file `design/k1/W3_l0.06.txt` scaled by $s$.
- R2 (larger terminal set, clear gain over vanilla at $N = 25$):
  $K = 0.159$, $\lambda = 0.0096$, $\delta_s = 3.9 \times 10^{-4}$,
  $\mu = 5.06 \times 10^{-4}$, $\delta = 0.192$, $T = 0.05$ s, $N = 25$.
- R1 (more margin in the one-step increase, smaller $\omega$):
  $K = 0.266$, $\lambda = 0.0160$, $\delta_s = 6.6 \times 10^{-4}$,
  $\mu = 1.42 \times 10^{-3}$, $\delta = 0.288$, $T = 0.05$ s.

### 11.7 Consequences and limits

- $\eta$ and $M_2$ are sampled. A rigorous version needs an SOS or interval
  bound for both; the $k_1$ margin itself is positive on a $1601^2$ grid.
- The paper's Lemma 2 must be replaced by the second-order statement (11) for
  $T^\star \approx 52$ ms to be the certified bound; with the current proof the
  bound stays near $10^{-7}$ s.
- The certified period has only 4 % margin over $T$. A smaller $s$ gives more.
- The terminal controller is slow (0.03 to 0.17 m/s), so the guaranteed
  increase per step is small, $\omega$ is of order $10^5$ and the bound of
  Corollary 1 is 40 s to 6000 s against an actual entry time near 1 s.
- With $N = 40$ vanilla MPC and RA-MPC are almost equal (288 versus 292). The
  advantage of the terminal set shows at shorter horizons.
- The input-constrained SOP stage is inactive for this example with the new
  design. Whether it matters for the manipulator has to be checked.
- Example 1 of the paper changes: $\psi$, the certificate parameters, all
  Dubins figures and counts.

### 11.8 Input bounds (setting R2 fixed)

What the theory needs: the terminal controller must respect the bounds on the
terminal set. For R2, $\lvert\omega\rvert \le 0.16$ and $\lvert a\rvert \le 0.006$ on
$\{\Psi \ge \delta\}$ ($0.2$ and $0.01$ on all of $\{\Psi \ge 0\}$), so any
$\omega_{\max} \ge 0.2$, $a_{\max} \ge 0.01$ is admissible. Under the time scaling
(14), $\omega \propto s$ and $a \propto s^2$: input bounds and the sampling period
both limit the speed scale $s$ of the terminal controller, and at
$T = 0.05$ s the sampling period is the binding one.

The bounds are active in the MPC closed loop, not in the terminal
controller. Study with the certificate, $\delta$, cost and $T = 0.05$ s fixed
(`design/eval_bounds.py`; continuous plant, 300 initial states; "saturated"
is the share of applied inputs within 2 % of the bound; vanilla is first-step
feasibility of the terminal-target MPC):

| $(\omega_{\max}, a_{\max})$ | vanilla, $N=25$ | RA-MPC, $N=25$ | saturated $\omega$ / $a$ | median $t_{\text{entry}}$ | vanilla, $N=40$ | RA-MPC, $N=40$ | saturated $\omega$ / $a$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| (5, 5) | 245 | 280 | 18 % / 79 % | 1.15 s | 288 | 292 | 8 % / 19 % |
| **(3, 3)** | **147** | **233** | 31 % / 93 % | 1.40 s | 275 | 286 | 22 % / 62 % |
| (2, 2) | 95 | 195 | 38 % / 95 % | 1.70 s | 209 | 254 | 33 % / 81 % |
| (2, 1) | 51 | 167 | 27 % / 83 % | 2.50 s | 118 | 206 | 28 % / 96 % |
| (1, 2) | 61 | 183 | 80 % / 85 % | 2.05 s | 166 | 225 | 73 % / 84 % |
| (1, 1) | 27 | 153 | 51 % / 81 % | 2.70 s | 97 | 183 | 54 % / 96 % |
| (1, 0.5) | 12 | 85 | 24 % / 61 % | 3.75 s | 52 | 148 | 36 % / 90 % |
| (0.5, 0.5) | 7 | 50 | 73 % / 59 % | 4.50 s | 25 | 128 | 80 % / 83 % |

In all 16 cases: no infeasible solve after $k = 0$, no violation of the cost
decrease, no safety violation. RA-MPC reaches the target from every
vanilla-feasible state except 0 to 4 per case.

Reading:

- Tighter bounds widen the gap. Vanilla MPC must hit the target within the
  horizon; RA-MPC only needs the terminal set and then keeps going (median
  entry time exceeds the horizon $NT = 1.25$ s for the tighter bounds).
- At (5, 5) with $N = 40$ the two methods are equal and the bounds are
  rarely active; this setting does not show the method.
- Suggested: $(\omega_{\max}, a_{\max}) = (3, 3)$ with $N = 25$ (233 versus 147,
  78 % versus 49 %), or $(2, 2)$ with $N = 40$ (254 versus 209). Both keep the
  bounds active in more than half of the runs and are plausible for a small
  ground vehicle moving at about 1 m/s.

### 11.9 The continuous controller alone as a baseline

`design/eval_ura.py`: $u^c$ applied in continuous time without MPC.

| start | R2 ($s = 0.159$) | unscaled ($s = 1$) |
| --- | --- | --- |
| the 300 initial states: inside $\{\Psi \ge 0\}$ | 0 | 4 |
| the 300 initial states: reached / left the safe set | 5 / 292 | 23 / 274 |
| 60 states of $\{\Psi \ge 0.2 M_\Psi\}$: reached | 60 (median 19.8 s, max 39.5 s) | 60 (median 3.1 s, max 6.3 s) |
| largest inputs on those runs | 0.11, 0.006 | 0.69, 0.23 |

The certified set is thin in the velocity directions (radius
$\sqrt{2\mu(\cdot)}$), and the backstepping law does not regulate the velocity
mismatch ($\dot e = \mu\nabla\psi + \tfrac{\lambda}{2} e$), so from the paper's
initial states (speeds up to 1 m/s in arbitrary direction) the continuous
controller alone fails. The 197/300 of the committed controller came from its
large uncertified gains. In the rigorous setting the three-way comparison of
the paper becomes 5 (controller alone), 147 (vanilla), 233 (RA-MPC) at
(3, 3), $N = 25$, and the baseline no longer illustrates input-bound
violations: it is slow and gentle, with a small domain.

### 11.10 Initial states for the three-way comparison (R2)

With the rigorous certificate the box used so far ($\theta \in [2\pi/3, 4\pi/3]$,
$v \in [-1, 1]$) has no state inside $\mathcal{X}_{RA} = \{\Psi \ge 0\}$, so the
comparison "controller alone, vanilla MPC, RA-MPC" has no group on which the
controller is guaranteed. Scripts: `design/eval_init.py`, `design/eval_init2.py`,
`design/plot_init.py`. All runs use the 300 positions of the previous sections,
the continuous plant and $T = 0.05$ s.

**States at rest are never inside $\mathcal{X}_{RA}$.** At $v = 0$,
$\Psi = \psi(y) - \delta_s/\lambda - \lVert k_1(y)\rVert^2/(2\mu)$. On a
$1601^2$ grid of the safe positions outside the target the largest value is
$-0.009$. The ratio $\lVert k_1\rVert^2/(2\mu)$ is invariant under the time
scaling, so the speed scale $s$ does not change this. It is also structural:
$u^c$ divides by $v$ ($\omega = \ddot y \cdot n / v$), so a certified set on
which $u^c$ is bounded must exclude $v = 0$. For R2 the speeds inside
$\mathcal{X}_{RA}$ lie in $[0.007, 0.113]$ m/s ($[0.026, 0.107]$ on the
terminal set).

**Settings tested** at $(\omega_{\max}, a_{\max}) = (3, 3)$, $N = 25$:

| initial velocity | inside $\mathcal{X}_{RA}$ | $u^c$ alone reaches | vanilla | RA-MPC |
| --- | --- | --- | --- | --- |
| at rest, heading in $[2\pi/3, 4\pi/3]$ | 0 | not defined | 161 | 233 |
| at rest, any heading | 0 | not defined | 157 | 235 |
| speed 0.02, any heading | 1 | 43 (256 leave the safe set) | 154 | 235 |
| speed 0.05, any heading | 6 | 42 (258 leave the safe set) | 150 | 236 |
| $\dot y(0) = k_1(y) + w$, $\lVert w\rVert \le 0.02$ | 136 | 259 | 165 | 240 |
| $\dot y(0) = k_1(y) + w$, $\lVert w\rVert \le 0.05$ | 25 | 155 | 164 | 244 |
| $\dot y(0) = k_1(y) + w$, $\lVert w\rVert \le 0.1$ | 4 | 63 | 165 | 234 |

In every row all states inside $\mathcal{X}_{RA}$ reach the target under
$u^c$ (1/1, 6/6, 136/136, 25/25, 4/4), and RA-MPC has no infeasible solve
after $k = 0$. A small initial speed does not help the controller alone: what
matters is that the velocity is close to $k_1(y)$, not that it is small. The
single family around the field mixes certified and uncertified successes of
$u^c$ (259 with only 136 guaranteed), which makes the count hard to interpret.

**Two initial sets on the same 300 positions** (superseded by Section 11.11:
the comparison has to use one initial set for all methods).

- $S_{RA}$: for the 234 positions with $\psi(y) > \delta_s/\lambda$, the velocity
  is drawn uniformly from the disc
  $\lVert \dot y - k_1(y)\rVert^2 \le 2\mu(\psi(y) - \delta_s/\lambda)$. All 234
  states are in $\mathcal{X}_{RA}$, 65 of them in the terminal set
  $\mathcal{X}_{RA}^{\delta}$.
- $S_0$: all 300 positions at rest with heading uniform in $[-\pi, \pi]$. No
  state is in $\mathcal{X}_{RA}$ and $u^c$ is not defined.

Controller alone on $S_{RA}$: 234/234 in continuous time and 234/234 under
zero-order hold with $T = 0.05$ s, median time to the target 26.1 s (max
99.6 s), largest inputs 0.14 rad/s and 0.006 m/s$^2$. MPC results:

| $(\omega_{\max}, a_{\max})$, $N$ | $S_{RA}$: $u^c$ alone | $S_{RA}$: vanilla | $S_{RA}$: RA-MPC | $S_0$: $u^c$ alone | $S_0$: vanilla | $S_0$: RA-MPC |
| --- | --- | --- | --- | --- | --- | --- |
| (3, 3), 25 | 234 | 163 | 226 | n/a | 159 | 231 |
| **(3, 3), 30** | **234** | **206** | **234** | **n/a** | **206** | **268** |
| (3, 3), 40 | 234 | 234 | 234 | n/a | 286 | 298 |
| (2, 2), 25 | 234 | 125 | 204 | n/a | 102 | 208 |
| (2, 2), 30 | 234 | 158 | 224 | n/a | 144 | 226 |
| (2, 2), 40 | 234 | 234 | 234 | n/a | 246 | 285 |

In all 12 MPC runs: no infeasible solve after $k = 0$; RA-MPC succeeds from all
65 states of the terminal set and from all 133 positions with
$\psi(y) - \delta_s/\lambda \ge \delta$ in both sets, at rest included; vanilla
succeeds where RA-MPC fails in one state only ($S_0$, (2, 2), $N = 40$). Median
time to the target under RA-MPC is 1.15 to 1.95 s.

Reading:

- $S_{RA}$ validates the certificate: every state reaches under $u^c$, in
  continuous time and under sampling, with inputs far below the bounds, but
  slowly (26 s against 1.2 s for RA-MPC at (3, 3), $N = 30$).
- RA-MPC is guaranteed from $\mathcal{X}_{RA}^{\delta}$, not from all of
  $\mathcal{X}_{RA}$. At $N = 25$ it misses 8 states of $S_{RA}$, all at the far
  tips of the certified region with $\psi - \delta_s/\lambda \le 0.013$, from
  which $\mathcal{X}_{RA}^{\delta}$ is out of reach in 1.25 s. At $N = 30$ it
  covers all of $S_{RA}$.
- $S_0$ shows the enlargement: none of these states is in $\mathcal{X}_{RA}$,
  yet RA-MPC succeeds from 268 of 300 against 206 for vanilla MPC at
  (3, 3), $N = 30$. The remaining 32 lie on the far side where no certificate
  exists (`design/figs/init_sets_3_3_N30.png`).
- Setting suggested at this stage (see Section 11.11 for the one adopted): both sets, $(\omega_{\max}, a_{\max}) = (3, 3)$,
  $N = 30$. Then RA-MPC contains both baselines on both sets. $N = 25$ gives a
  larger gap to vanilla MPC (226 against 163, 231 against 159) at the price of
  the 8 states above.

### 11.11 One initial set for all three methods (R2)

Requirement: the three methods are compared on the same initial states; the
positions come from $\mathcal{X}_S \setminus \mathcal{X}_T$ and need not lie in
$\mathcal{X}_{RA}$; the initial velocity is drawn uniformly from a range that
leaves the controller alone enough admissible starts. Rule used
(`design/eval_single.py`, `build(r0)`):

- positions: the 300 points drawn uniformly from the safe ring outside the
  target (unchanged, `eval_mpc.py`);
- velocity: $\dot y(0) = k_1(y(0)) + w$ with $w$ uniform in the disc of radius
  $r_0$, then $\theta(0)$ and $v(0)$ from $\dot y(0)$. The same unit-disc draws
  are used for every $r_0$, so the sets are nested.

The velocity disc of $\mathcal{X}_{RA}$ at a position has radius
$\sqrt{2\mu(\psi(y) - \delta_s/\lambda)}$, at most 0.030 over the 300 positions
(median 0.015), so $r_0 = 0.03$ is the smallest radius whose tube contains
every velocity that $\mathcal{X}_{RA}$ admits.

**Radius scan, controller alone (URA, continuous time).**

| $r_0$ | inside $\mathcal{X}_{RA}$ | inside $\mathcal{X}_{RA}^{\delta}$ | URA reaches | URA leaves the safe set | URA median time |
| --- | --- | --- | --- | --- | --- |
| 0.01 | 191 | 107 | 291 | 9 | 31.7 s |
| 0.02 | 124 | 67 | 267 | 33 | 29.9 s |
| 0.03 | 69 | 32 | 229 | 71 | 25.7 s |
| 0.04 | 36 | 20 | 192 | 108 | 24.3 s |
| 0.05 | 23 | 11 | 149 | 151 | 21.2 s |

Every state inside $\mathcal{X}_{RA}$ reaches the target under URA in all five
rows.

**Three methods on the same 300 states** (continuous plant, $T = 0.05$ s;
vanilla is first-step feasibility as before; A minus B counts the states where
A succeeds and B does not):

| $r_0$ | $(\omega_{\max}, a_{\max})$, $N$ | URA | vanilla | RA-MPC | RA-MPC minus vanilla | vanilla minus RA-MPC | RA-MPC minus URA | URA minus RA-MPC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.02 | (3, 3), 25 | 267 | 166 | 238 | 72 | 0 | 4 | 33 |
| 0.02 | (3, 3), 30 | 267 | 216 | 281 | 65 | 0 | 21 | 7 |
| 0.02 | (2, 2), 40 | 267 | 256 | 282 | 26 | 0 | 21 | 6 |
| **0.03** | **(3, 3), 25** | **229** | **166** | **236** | **70** | **0** | **22** | **15** |
| 0.03 | (3, 3), 30 | 229 | 214 | 272 | 58 | 0 | 45 | 2 |
| 0.03 | (2, 2), 40 | 229 | 252 | 290 | 38 | 0 | 62 | 1 |

In all six runs every vanilla-feasible state is RA-MPC feasible, RA-MPC
succeeds from every state inside $\mathcal{X}_{RA}$ (124/124 and 69/69), and no
solve is infeasible after $k = 0$.

**URA against RA-MPC on the common successes.** Cost is the accumulated stage
cost $\int (5\lVert y - y^\star\rVert^2 + 0.05\lVert u\rVert^2)\,dt$ up to the
entry into the target.

| $r_0$ | $(\omega_{\max}, a_{\max})$, $N$ | common | median time URA / RA-MPC | median cost URA / RA-MPC | peak inputs URA | peak inputs RA-MPC |
| --- | --- | --- | --- | --- | --- | --- |
| 0.02 | (3, 3), 25 | 234 | 25.2 s / 1.30 s | 261 / 15.5 | 0.16, 0.006 | 3.0, 3.0 |
| 0.02 | (3, 3), 30 | 260 | 29.1 s / 1.33 s | 309 / 17.1 | 0.31, 0.006 | 3.0, 3.0 |
| 0.02 | (2, 2), 40 | 261 | 29.1 s / 1.40 s | 310 / 19.8 | 0.16, 0.006 | 2.0, 2.0 |
| 0.03 | (3, 3), 25 | 214 | 24.2 s / 1.25 s | 219 / 14.2 | 0.15, 0.007 | 3.0, 3.0 |
| 0.03 | (3, 3), 30 | 227 | 25.6 s / 1.15 s | 263 / 14.6 | 0.15, 0.007 | 3.0, 3.0 |
| 0.03 | (2, 2), 40 | 228 | 25.7 s / 1.30 s | 264 / 16.8 | 0.15, 0.007 | 2.0, 2.0 |

Reading:

- RA-MPC against vanilla MPC: on the same states RA-MPC is feasible wherever
  vanilla MPC is, plus 26 to 72 more. The gap is largest for the short
  horizon.
- RA-MPC against URA: neither set of successes contains the other. URA has no
  horizon, so it also succeeds from far positions that RA-MPC cannot connect
  to the terminal set within $NT$; RA-MPC succeeds from states whose velocity
  is too far from $k_1$ for URA, which then leaves the safe set. The order of
  the two counts depends on $r_0$ and $N$.
- What the URA comparison shows is the optimisation under input bounds: about
  20 times faster and 15 to 18 times cheaper on the common successes, with
  inputs at the bounds, whereas URA uses at most 10 % of $\omega_{\max}$ and
  0.2 % of $a_{\max}$.
- URA never exceeds the bounds in the rigorous design (it must be admissible on
  the terminal set, and the time scaling that makes $T = 0.05$ s certifiable
  makes it slow). The sentence of the paper that the unconstrained controller
  "severely violates actuator limits" does not hold any more and has to become
  a statement about speed and cost.
- Suggested: $r_0 = 0.03$ with $(3, 3)$, $N = 25$ (URA 229, vanilla 166,
  RA-MPC 236; `design/figs/single_r0.03_3_3_N25.png`). $r_0 = 0.02$ gives URA
  more starts (267, of which 124 inside $\mathcal{X}_{RA}$) with the same MPC
  picture (166 against 238).

### 11.12 Final setting and data for the paper

The setting adopted for the paper, the verification of Lemma 2, Theorem 2 and
Corollary 1, the sampling-period sweep and the paper figures are in
`PAPER_DUBINS.md`, produced by `final/`. Since Section 11.11 the design
changed in three ways: the terminal controller's time scale is chosen from
the input bounds and the sampled one-step increase at $T = 0.05$ s
($s = 0.764$ instead of 0.159, derivation in Section 7 of
`TAC_EXPERIMENT_UPDATE.md`); the URA baseline is the same synthesis at
$s = 5$ and the initial velocities are drawn around its field; the state
constraint is imposed at the RK4 sub-step nodes with a small back-off and
both MPCs share model, constraints and initial-guess fallbacks. Counts:
URA 229 (158 runs over the bounds), vanilla 213, RA-MPC 263.

## 12. Reproduction

```
cd experiments/strict_xra
matlab -batch "run('run_lmin.m')"    # Section 2.2 (lambda sweep, degree 4)
matlab -batch "run('run_deg.m')"     # degree test
matlab -batch "run('run_kmax.m')"    # Section 4.2, unscaled, K = 1
matlab -batch "run('run_scaled.m')"  # Section 4.2, scaled
matlab -batch "run('run_full.m')"    # Section 5, full synthesis
matlab -batch "run('run_alpha.m')"   # Section 9.2, weight (2.5 - y1)^2
matlab -batch "run('run_full_alpha.m')" # Section 9.3
python lemma2_tightness.py out/rigorous_alpha2_mu0.1_vanilla_k1.py 0.005   # Section 10
python t_vs_delta.py                       # Section 10.4
cd design
python screen_weights.py && python screen_w3.py   # Section 11.2
matlab -batch "run('run_design.m')"        # Section 11.3 (writes k1/*.txt)
python validate.py                         # cert.py versus the MATLAB export
python scan.py W3 0.06 0.02,0.05 0.2,0.3,0.4,0.5   # Section 11.4
python eval_mpc.py 25,40                   # Section 11.5 (certified settings, both horizons)
python eval_bounds.py 25,40                # Section 11.8 (input bounds, R2)
python eval_ura.py                         # Section 11.9 (controller alone)
python eval_init.py                        # Section 11.10 (rest, near rest, field-relative)
python eval_init2.py && python eval_init2.py sweep   # Section 11.10 (two initial sets)
python plot_init.py 3_3 30                 # outcome map, figs/init_sets_3_3_N30.png
python eval_single.py                      # Section 11.11 (one initial set, three methods)
python plot_single.py 0.03 3_3 25          # figs/single_r0.03_3_3_N25.png
# exports land in controllers/tmp_*.py: move them to out/ before anything else
LAM=0.02 python check_full.py out/rigorous_mu0.1_vanilla_k1.py
```
