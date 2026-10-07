# Dubins car: final setting, data and figures for the paper

Setting of 2026-10-07: shape of $k_1$ from (12) at the normalised scale, terminal
controller at the time scale $s = 0.764$ chosen from the input bounds and the
sampled one-step increase at $T = 0.05$ s (`TAC_EXPERIMENT_UPDATE.md`
Section 7), URA = the same synthesis at the performance scale $s = 5$, one
initial set for all methods with velocities around the URA's field,
$(\omega_{\max}, a_{\max}) = (3, 3)$, $N = 25$. Everything below is produced by
the scripts in `final/` (Section 9) on the continuous plant. The paper and the
notebooks are not yet changed. Earlier settings are listed in Section 10.

Contents: 1 setting, 2 three-method comparison, 3 Lemma 2, 4 Theorem 2 and
Corollary 1, 5 sampling-period sweep, 6 figures for the paper, 7 places to
update in the paper, 8 open points, 9 reproduction, 10 earlier settings.

## 1. Setting (Example 1)

### 1.1 Problem data

| item | value | remark |
| --- | --- | --- |
| output | $y = (x_1, x_2)$ | unchanged |
| safe set | $\mathcal{X}_S = \{\psi(y) > 0\}$, the ring $4 < y_1^4 + y_2^4 < 16$ | set unchanged, $\psi$ changed |
| $\psi$ | $\psi(y) = \alpha(y)\,h(y)/661.409$ | replaces $10^{-3}(300 - \phi)\,h$ |
| $h$ | $h(y) = -(y_1^4 + y_2^4 - 16)(y_1^4 + y_2^4 - 4)$ | unchanged |
| $\alpha$ | $\alpha(y) = (2.5 - y_1)^2 + 0.3\,y_1 y_2^2\,(2.5 - y_1)$ | positive on the ring, minimum 0.25 |
| maximum of $\psi$ | 1, at $y = (-1.794, 0) \in \mathcal{X}_T$ | unique maximum inside the target, as the assumption on the sets asks |
| target set | $\mathcal{X}_T = \{\phi(y) < 0\}$, $\phi(y) = 4(y_1 + 1.7)^2 + y_2^2 - 0.4$ | the paper prints $y_1^2 + 4(y_2 + 1.7)^2 - 0.4$, indices swapped with respect to the code and the figures |
| input set | $\lvert\omega\rvert \le 3$, $\lvert a\rvert \le 3$ | was 5 |

### 1.2 Shape of $k_1$ and the scaled family

Program (12) is solved once at the normalised scale ($\lambda_0 = 0.06$,
$\lvert k_{1,i}\rvert \le 1$ on $\{\psi \ge 0\}$, scaled outputs $y/2$, degree 8;
`design/run_design.m`). The base design is then scaled in time by $s$
(`TAC_EXPERIMENT_UPDATE.md` Section 7): $k_1 \to s k_1$, $\lambda \to s\lambda_0$,
$\delta_s \to s\delta_s$, $\mu \to s^2\mu_0$ with $\mu_0 = 0.02$. The certificate
$\Psi$, $\delta_s/\lambda = 0.0410$, $M_\Psi = 0.959$, the position part of
$\mathcal{X}_{RA}$ (76.5 % of $\mathcal{X}_S \setminus \mathcal{X}_T$) and of
$\mathcal{X}_{RA}^{\delta}$ (37.9 %) do not depend on $s$; the velocities,
the inputs ($\omega \propto s$, $a \propto s^2$) and the speed of the closed
loop do.

### 1.3 Choice of the time scale of the terminal controller

`final/solve_scale.py`, `final/verify_inputs_grid.py`, $T = 0.05$ s,
$\delta = 0.2\,M_\Psi = 0.192$:

| bound on $s$ | how | value |
| --- | --- | --- |
| input bounds on $\mathcal{X}_{RA}^{\delta}$, sampled suprema (61 322 states) | $s_U = \min(3/\sup\lvert\omega_1\rvert,\ \sqrt{3/\sup\lvert a_1\rvert}) = \min(3/0.994, \sqrt{3/0.237})$ | 3.02 |
| input bounds, analytic bound over each velocity fibre on a $2001^2$ position grid (no velocity sampling) | $\lvert\omega_1\rvert \le 1.73$, $\lvert a_1\rvert \le 0.40$ | 1.74 |
| one-step increase of $\Psi$ under ZOH at the period $T$ (sampled, 61 322 states) | positive for base arcs up to $\tau_\varepsilon = 47.7$ ms, $s_\varepsilon = \tau_\varepsilon/T$ | 0.95 |
| second-order bound of Lemma 2, for reference | $s_T = T^\star_1/T$, $T^\star_1 = 2\eta_1/M_{2,1} = 8.1$ ms | 0.16 |
| closed-form bound of the current proof, for reference | $T^\star = 0.35$ ms at $s = 0.764$ | about 0.005 |

The period, not the input bounds, limits the terminal controller. Chosen
$s = 0.8\,s_\varepsilon = 0.764$ (safety factor 0.8 on the sampled threshold;
within both input bounds by a factor 2.3 and 4). The conservative bounds are
not used for the choice; they are reported because Lemma 2 and Theorem 2 are
stated with them.

### 1.4 Terminal controller (the design used by RA-MPC)

| item | value |
| --- | --- |
| scale | $s = 0.764$ |
| $\lambda$, $\delta_s$, $\mu$ | $0.0459$, $1.88 \times 10^{-3}$, $0.0117$ |
| largest field speed | $\max\lvert k_{1,i}\rvert = 0.478$ m/s |
| $\mathcal{X}_{RA}$ | speeds 0.03 to 0.54 m/s; $u^c$ up to 1.02 rad/s, 0.145 m/s$^2$ |
| $\mathcal{X}_{RA}^{\delta}$ | speeds 0.12 to 0.51 m/s; $u^c$ up to 0.72 rad/s, 0.138 m/s$^2$ |
| one-step increase at $T = 0.05$ s on $\mathcal{X}_{RA}^{\delta}$ | $\varepsilon = 1.07 \times 10^{-4}$ (76 809 arcs), all increasing |
| offline time | one SOS solve, 1.9 s, plus the sampled quantities above |

### 1.5 URA baseline

The same $k_1$ at the performance scale $s = 5$: $\lambda = 0.30$, $\mu = 0.50$,
field speeds up to 3.13 m/s. It is the synthesis without regard to the input
bounds or the period; $s = 5$ exceeds $s_U$, so its inputs leave the bounds on
part of its runs (Section 2).

### 1.6 MPC

| item | value |
| --- | --- |
| sampling period, horizon | $T = 0.05$ s, $N = 25$ |
| stage cost | $\ell(x, u) = 5\lVert y - y^\star\rVert^2 + 0.05\lVert u\rVert^2$, $y^\star = (-1.7, 0)$, $\underline{\ell} = 0.5$ |
| terminal cost | $V_f = \omega\,(M_\Psi - \Psi)$, $\omega = 1.2\,\sup\ell/\varepsilon$ with $\varepsilon$ the sampled one-step increase (seed of `terminal_stats`) and $\sup\ell = 46.1$ |
| vanilla MPC | same model, stage cost, constraints and solver; terminal constraint $\hat x(NT) \in \mathcal{X}_T$, terminal cost $80\lVert y_N - y^\star\rVert^2$ |
| prediction model | RK4 with 4 sub-steps per period (12.5 ms) |
| state constraint | at every sub-step node with a back-off of about 1.2 mm |
| initial guesses at $k = 0$ (both MPCs) | cold start; on failure a steer-to-target rollout; RA-MPC also the vanilla solution and a continuation over the terminal set's scale ($s/4$, $s/2$, $s$) |
| plant | continuous time, input held over $T$ (RK4 at 5 ms) |

### 1.7 Initial states

One set of 300 states for all methods (`final/config.py`, `initial_set()`
with `SET_SCALE = 5`).

- Positions: uniform in $\mathcal{X}_S \setminus \mathcal{X}_T$ (unchanged).
- Velocities: $\dot y(0) = k_1^{(5)}(y(0)) + w$ with $k_1^{(5)}$ the URA's field
  and $w$ uniform in the disc of radius $r_0 = 0.943$ m/s, the largest velocity
  radius of the URA's $\mathcal{X}_{RA}$. Speeds range from 0.15 to 3.6 m/s.
- 69 states lie in the URA's $\mathcal{X}_{RA}$; none lies in the terminal
  controller's $\mathcal{X}_{RA}$, which is 6.5 times slower.
- Success: the closed loop enters $\mathcal{X}_T$ and stays in $\mathcal{X}_S$
  until then (URA run for 600 s, MPC for 20 s; vanilla counted by first-step
  feasibility, its closed loop then run from all those states).

## 2. Three-method comparison

### 2.1 Success counts (subfigure captions)

| | URA ($s = 5$) | vanilla MPC | RA-MPC |
| --- | --- | --- | --- |
| **reach $\mathcal{X}_T$ safely, of 300** | **229** | **213** | **263** |
| runs that exceed the input bounds | 170 (158 among the successes) | 0 | 0 |
| peak $\lvert\omega\rvert$, $\lvert a\rvert$ over all runs | 5.4 rad/s, 6.7 m/s$^2$ | 3, 3 | 3, 3 |
| inside the URA's $\mathcal{X}_{RA}$ (69 states) | 69 | | |

| pair | count |
| --- | --- |
| RA-MPC succeeds, vanilla MPC does not | 50 |
| vanilla MPC succeeds, RA-MPC does not | 0 |
| RA-MPC succeeds, URA does not | 38 |
| URA succeeds, RA-MPC does not | 4 |
| all three succeed | 199 |
| none succeeds | 33 |

- URA: 229 reach, 71 leave the safe set; on 158 of the 229 the inputs exceed
  the bounds (mostly $a$, up to 6.7 m/s$^2$; $\omega$ up to 5.4 rad/s).
- RA-MPC: the 37 failures are all infeasible at $k = 0$; no solve is
  infeasible afterwards (Section 4). The terminal controller alone reaches
  from 46 of the 300 states (the set is 6.5 times faster than its field).

![outcomes](final/figs/verify_outcomes.png)

### 2.2 URA against RA-MPC: time, cost, inputs

On the 225 states where both succeed (cost: accumulated stage cost up to the
entry into $\mathcal{X}_T$):

| quantity | URA | RA-MPC |
| --- | --- | --- |
| time to the target, median (maximum) | 0.80 s (4.3 s) | 0.70 s (2.4 s) |
| accumulated cost, median | 8.7 | 7.0 |
| largest $\lvert\omega\rvert$ | 4.6 rad/s | 3.0 rad/s (bound) |
| largest $\lvert a\rvert$ | 6.7 m/s$^2$ | 3.0 m/s$^2$ (bound) |
| share of applied inputs within 2 % of the bound | | $\omega$ 28 %, $a$ 77 % |

On the 199 states where all three succeed the median times are 0.80 s (URA),
0.60 s (vanilla MPC) and 0.65 s (RA-MPC), costs 6.1, 5.8 and 6.2. The point of
the comparison: URA is as fast only by exceeding the actuator limits, on
two thirds of its successful runs; RA-MPC matches its speed within the limits
and succeeds from 34 more states.

Shared initial state of the trajectory panels and of the input plot: state 6,
$x_0 = (1.330, 1.776, -2.820, 1.733)$, 3.5 from the target, the farthest state
from which all three methods succeed and URA exceeds the bounds. Times to
the target: vanilla MPC 1.20 s, URA 1.47 s with peaks $\omega$ 3.5 rad/s and
$a$ 5.7 m/s$^2$ (0.66 s above the bounds), RA-MPC 1.55 s within the bounds.

### 2.3 Computation time (Table I)

Measured on this machine (Apple M5, one process, 12 runs per method), not on
the i5-13400 named in the paper; to be measured again there.

| method | offline | online per step, mean (median) |
| --- | --- | --- |
| vanilla MPC | none | 47.0 ms (45.0 ms) |
| URA | 1.9 s | 0.23 ms |
| RA-MPC | 1.9 s | 69.7 ms (68.9 ms) |

- Both MPC schemes use the 4-sub-step model with the state constraint at the
  sub-step nodes and tight solver tolerances, hence slower than the 28.9 ms
  and 30.8 ms of the current table. Under the parallel evaluation the means
  are 97 ms and 141 ms.
- The first solve of RA-MPC takes 114 ms on average because of the fallback
  initial guesses (1.6); later solves 66.5 ms. The RA-MPC step exceeds
  $T = 0.05$ s on this machine.
- The offline time is one SOS solve (shared with URA) plus the sampled
  quantities of 1.3 and 1.6, a few seconds.

## 3. Verification of Lemma 2

Terminal controller at $s = 0.764$; `final/run_lemma2.py`, states sampled
inside the sets (half on the velocity boundary).

### 3.1 At $T = 0.05$ s

| check | result |
| --- | --- |
| one ZOH step from 76 896 states of $\mathcal{X}_{RA}^{\delta}$ (76 809 arcs outside $\mathcal{X}_T$): $\Psi$ increases | all, smallest increase $1.07 \times 10^{-4}$ |
| same arcs: end state in $\mathcal{X}_{RA}^{\delta}$, arc inside $\mathcal{X}_{RA}$ | all |
| closed loop from 600 states of $\mathcal{X}_{RA}^{\delta}$: reach $\mathcal{X}_T$ inside $\mathcal{X}_S$ | 600 of 600 |
| same runs: $\Psi$ increasing at every sampling instant | 600 of 600, smallest step $2.0 \times 10^{-4}$ |
| closed loop from 600 states of $\mathcal{X}_{RA}$ | 600 of 600 reach; $\Psi$ increasing in 597 |

### 3.2 Against the sampling period

One ZOH step from 9 295 states of $\mathcal{X}_{RA}^{\delta}$ (arcs entering
$\mathcal{X}_T$ left out) and the sampled-data closed loop from 600 states:

| $T$ [s] | smallest increase | arcs with increasing $\Psi$ | end in $\mathcal{X}_{RA}^{\delta}$ | closed loop from $\mathcal{X}_{RA}^{\delta}$: reach | $\Psi$ increasing at every sample | from $\mathcal{X}_{RA}$: reach |
| --- | --- | --- | --- | --- | --- | --- |
| 0.01 | $+1.2 \times 10^{-4}$ | 100 % | 100 % | 600 | 600 | 600 |
| 0.02 | $+2.1 \times 10^{-4}$ | 100 % | 100 % | 600 | 600 | 600 |
| 0.05 | $+2.9 \times 10^{-4}$ | 100 % | 100 % | 600 | 600 | 600 |
| 0.075 | $+2.1 \times 10^{-4}$ | 100 % | 100 % | | | |
| 0.1 | $-4.3 \times 10^{-4}$ | 99.97 % | 99.97 % | 600 | 597 | 600 |
| 0.2 | $-1.1 \times 10^{-2}$ | 98.6 % | 98.8 % | 600 | 570 | 600 |
| 0.5 | $-0.15$ | 86.8 % | 90.1 % | 600 | 405 | 600 |
| 1 | $-0.95$ | 66.1 % | 73.4 % | 600 | 233 | 582 |
| 2 | $-6.8$ | 30.5 % | 33.4 % | 410 | 262 | 323 |

### 3.3 Sampling-period bounds at $s = 0.764$

| bound on $T$ | value |
| --- | --- |
| closed form of the current proof ($M_g = 3.6,\ 10.9$; $\ell_{u^c} = 7.55,\ 0.59$; $M_v = 0.75$) | 0.35 ms |
| second-order bound $2\eta/M_2$ ($\eta = 9.6 \times 10^{-3}$, certified rate $\lambda\delta = 8.8 \times 10^{-3}$; $M_2 = 1.85$) | 10.4 ms (9.5 ms with $\lambda\delta$) |
| observed: $\Psi$ increases over one step from every sampled state | up to 75 ms |
| observed: the closed loop reaches $\mathcal{X}_T$ safely from every sampled state | up to 1 s |

Reading: with $T = 0.05$ s the conclusion of Lemma 2 holds with a factor 1.5
to the first one-step violation and a factor 20 to the first closed-loop
failure; both bounds are below $T$ (5 times for the second-order one), so the
paper's hypothesis $\delta > \beta(T)/\lambda$ is not met by the bounds and is
met by the sampled check, as decided for the experiments.

## 4. Verification of Theorem 2 and Corollary 1

Data: the 300 RA-MPC runs of Section 2 (`final/run_main.py`).

![theorem 2](final/figs/verify_theorem2.png)

### 4.1 Checks

| claim | check | result |
| --- | --- | --- |
| hypothesis (22) on the terminal cost | $\omega\,[\Psi(x^s(t_{k+1})) - \Psi(x^s(t_k))] \ge \ell(x^s(t_k), u^c)$ on 76 490 sampled states of $\mathcal{X}_{RA}^{\delta}$ | all; $\omega$ times the smallest increase is 74.8 against $\sup \ell = 46.1$ |
| terminal controller admissible | $u^c$ of the appended candidate step within the bounds | every step of every run |
| recursive feasibility | solves after $k = 0$ in the 263 runs feasible at $k = 0$ | 5 217 solves, none infeasible |
| reaches $\mathcal{X}_T$ | runs feasible at $k = 0$ that enter the target | 263 of 263 |
| strictly inside $\mathcal{X}_S$ | smallest $h$ along the continuous closed loop before the entry, 2.5 ms resolution | 0.165 > 0, no run with a negative value |
| cost decrease | $J^*_k - J^*_{k+1}$ over 4 954 steps | smallest 1.99; never below $\underline{\ell} = 0.5$; never below $\ell(x_k, u_k)$ beyond solver tolerance (ratio at least 0.999 999) |
| adaptive horizon | $N$ along each run | nonincreasing in 263 of 263, shortened in 256, smallest value 1 |
| Corollary 1 | entry time against $(\lfloor J^*_0/\underline{\ell}\rfloor + 1)\,T$ | holds in 263 of 263 |

- Entry times: median 0.75 s, maximum 2.40 s. Bound of Corollary 1: between
  0.40 s and 53 541 s, median 404 s. It is loose by three orders of magnitude
  for most states because $V_f$ with $\omega = 7.0 \times 10^5$ dominates
  $J^*_0$, and tight only near the target.
- The candidate was kept instead of the solver output in 2 959 of the 5 217
  solves. Once the horizon is being shortened the tail of the previous plan is
  optimal, the decrease equals $\ell(x_k, u_k)$, and the points of panel (c)
  lie on the diagonal.

### 4.2 Safety between sampling instants

With the state constraint imposed at the sampling instants only, as in all
earlier runs and in the L-CSS code, about half of the successful RA-MPC runs
of the earlier setting left $\mathcal{X}_S$ between two sampling instants, by
at most 1.2 mm. The MPC problem of Theorem 2 constrains
$\hat x(\tau \mid k)$ for all $\tau$, so this was an implementation gap. The
final runs impose the constraint at the four RK4 sub-step nodes of every
period with a back-off of 0.05 in the units of the safe-set function (about
1.2 mm at the inner boundary), for both MPC schemes.

## 5. Sampling-period sweep

`final/run_sweep.py`: $T$ varies with $NT = 1.25$ s, the terminal design is kept
at its base values ($s = 0.764$, $\delta = 0.192$, $\omega = 7.0 \times 10^5$),
model sub-steps and constraint nodes stay at 12.5 ms. Same 300 states.

| $T$ [s] | $N$ | vanilla MPC | $u^c$ alone, sampled-data | RA-MPC | $\rho$ of RA-MPC | infeasible after $k=0$ | steps violating the decrease | $\min \Delta\Psi$ on $\mathcal{X}_{RA}^{\delta}$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.025 | 50 | 213 | 46 | 265 | 88.3 % | 0 | 0 of 10 373 | $+1.5 \times 10^{-4}$ |
| 0.03125 | 40 | 213 | 46 | 264 | 88.0 % | 0 | 0 of 8 141 | $+1.6 \times 10^{-4}$ |
| **0.05** | **25** | **213** | **46** | **263** | **87.7 %** | **0** | **0 of 4 954** | $+7.9 \times 10^{-5}$ |
| 0.0625 | 20 | 213 | 46 | 264 | 88.0 % | 0 | 0 of 3 941 | $-4.4 \times 10^{-5}$ |
| 0.125 | 10 | 213 | 47 | 264 | 88.0 % | 0 | 0 of 1 829 | $-1.6 \times 10^{-3}$ |
| 0.25 | 5 | 213 | 47 | 264 | 88.0 % | 0 | 0 of 835 | $-2.4 \times 10^{-2}$ |
| 0.3125 | 4 | 213 | 47 | 261 | 87.0 % | 0 | 0 of 649 | $-4.4 \times 10^{-2}$ |
| 0.625 | 2 | 209 | 48 | 248 | 82.7 % | 0 | 0 of 364 | $-0.27$ |

![sweep](final/figs/verify_sweep.png)

- For every $T$, every run feasible at $k = 0$ reaches the target, no run is
  unsafe or times out, and vanilla-feasible states are RA-MPC feasible
  without exception.
- Up to $T = 0.25$ s the success rate is constant (263 to 265); it falls to
  82.7 % at $T = 0.625$ s, through the set of states feasible at $k = 0$
  (coarser input parametrisation), not through recursive feasibility.
- The one-step increase of $\Psi$ on the terminal set is lost from
  $T = 0.0625$ s on (last column), so hypothesis (22) fails on part of the
  terminal set there, yet the closed loop never needs it: the candidate is
  only a fallback and the solver always found a feasible plan.
- The current sentence of the paper, that the drop of $\rho$ with $T$
  "directly corroborates" Lemma 2, is not supported: the drop sets in an
  order of magnitude beyond the period at which the lemma's mechanism (the
  one-step increase) stops holding.

## 6. Figures for the paper

Files in `final/figs/`, PNG and PDF, same names, size and style as the current
figures in `TAC/figs/`, so they can be copied over once confirmed.

### 6.1 Three method panels (figure* in `method.tex`)

| vanilla MPC, 213/300 | URA, 229/300 | RA-MPC, 263/300 |
| --- | --- | --- |
| ![vanilla](final/figs/fig_dubins_vanilla.png) | ![ura](final/figs/fig_dubins_uncon.png) | ![rampc](final/figs/fig_dubins_rampc.png) |

Blue dots succeed, red crosses fail, the dashed line is the closed loop from
the shared initial state, the grey streamlines are $k_1$.

### 6.2 Sweep and inputs (`fig_dubins_compare`)

![compare](final/figs/fig_dubins_compare.png)

Top: $\rho$ of RA-MPC against $T$ (Section 5), $\rho$ axis zoomed to the data,
dashed line at $T = 0.05$ s (the shaded certified range of the second-order
bound, $T < 10$ ms, lies left of the plotted range). Below: $\omega$ and $a$
from the shared initial state, URA in blue (1.47 s, outside the bounds between
0.24 s and 0.65 s and again at the end) and RA-MPC in purple (1.55 s, inside),
bounds dashed, linear time axis.

### 6.3 Variant with a farther start (`*_far` files)

URA and RA-MPC from state 297, $x_0 = (1.878, -1.306, -3.070, 0.389)$, 3.8 from
the target (URA 2.61 s with peaks 2.9 and 4.1, RA-MPC 2.10 s); vanilla from its
own farthest feasible state 65 at $(1.77, 0.73)$ (1.45 s).

| vanilla MPC, own start | URA, far start | RA-MPC, far start |
| --- | --- | --- |
| ![vanilla](final/figs/fig_dubins_vanilla_far.png) | ![ura](final/figs/fig_dubins_uncon_far.png) | ![rampc](final/figs/fig_dubins_rampc_far.png) |

![compare far](final/figs/fig_dubins_compare_far.png)

### 6.4 Optional: cost and horizon (`fig_dubins_cost`)

![cost](final/figs/fig_dubins_cost.png)

Three rows for the shared initial state: $J^*_k$, the decrease
$J^*_k - J^*_{k+1}$ against $\ell(x_k, u_k)$ and $\underline{\ell}$, and the
horizon used at step $k$. The plateau of $J^*_k$ is the terminal cost of the
shortened plan; the middle row shows the decrease directly.

## 7. Places to update in the paper

| file | place | change |
| --- | --- | --- |
| `Example_dubins.tex` | definition of $\psi$ and $\phi$ | new $\psi$ (1.1); fix the swapped indices in $\phi$ |
| `Example_dubins.tex` | "300 safe initial states randomly sampled outside the target set" | describe the initial set (1.7) |
| `Example_dubins.tex` | "actuator limits $\lvert u_{1,2}\rvert \le 5$" | 3; $N = 25$ and $T = 0.05$ s stay |
| `Example_dubins.tex` | certificate and controller data | the base design, the scale rule and the two scales (1.2 to 1.5) |
| `Example_dubins.tex` | URA "severely violates actuator limits and delays convergence" | exceeds the limits on 158 of 229 runs (up to 5.4 rad/s and 6.7 m/s$^2$); not slower: 0.80 s against 0.70 s (2.2) |
| `Example_dubins.tex` | sweep sentence on Lemma 2 | see Section 5, last point |
| `method.tex` | captions of the figure* | 213/300, 229/300, 263/300 instead of 245, 197, 255 |
| `method.tex` | Algorithm 1, offline part | the scale rule of Section 7 of `TAC_EXPERIMENT_UPDATE.md` if adopted |
| `experiment.tex` | Table I | Dubins columns, after measuring on the paper's machine (2.3) |
| `experiment.tex` | caption of `fig_dubins_compare` | zoomed $\rho$ axis, URA outside the bounds |
| `figs/` | `fig_dubins_vanilla`, `fig_dubins_uncon`, `fig_dubins_rampc`, `fig_dubins_compare` | replace by the files of `final/figs/` |

## 8. Open points

1. **$T = 0.05$ s against the bounds.** The second-order bound gives
   $T^\star = 10.4$ ms at $s = 0.764$ and the closed form of the current proof
   0.35 ms; the design uses the sampled threshold of the one-step increase
   (75 ms) with a factor 0.8 on $s$, as decided. The paper has to say that
   the period is validated by the sampled check and not by the bound, or
   Lemma 2 has to be tightened before the experiments are described through
   it.
2. **Form of Lemma 2 and of (12).** The time-scaling proposition (Section 7
   of the plan) with $\lambda$ fixed and the $K$ bound; the closed-form bound
   cannot certify any period used.
3. **Sampled quantities.** $s_U$ (sampled 3.02, analytic bound 1.74),
   $s_\varepsilon$, $\varepsilon$, $\sup\ell$ and hence $\omega$, $\eta$, $M_2$
   are extrema over samples (9 295 to 76 896 states), except the analytic
   $s_U$ bound and the certified rate $\lambda\delta$.
4. **Initial guesses.** Both MPCs try a steer-to-target rollout when the cold
   start fails; RA-MPC also the vanilla solution and a continuation over the
   terminal set's scale ($s/4$, $s/2$, $s$). Without it RA-MPC lost 35 far-side
   states to solver failures (224 instead of 263). The paper should mention
   the initialisation in one sentence.
5. **URA scale.** $s = 5$ is a design choice (cruise speed about 3 m/s); it
   exceeds $s_U \approx 3$, which is what makes the comparison. State it.
6. **Corollary 1 is loose** for this terminal weight (4.1).
7. **Online time above $T$** on this machine (2.3).
8. **Manipulator example**: not audited; the scale rule needs the gravity
   term (Section 7.3 of the plan).
9. **Notebooks** `example_dubins_car_*.ipynb` and `make_paper_figs_dubins.ipynb`
   still run the committed controller; `controllers/*.py` is untouched.

## 9. Reproduction

```
cd experiments/strict_xra/final
python solve_scale.py            # Section 1.3: s_U, s_eps, s_T  -> data/scale.json
python verify_inputs_grid.py     # Section 1.3: analytic input bound over X_RA^delta
python run_main.py               # Sections 2 and 4   -> data/main.pkl      (3 min)
python run_lemma2.py             # Section 3          -> data/lemma2.pkl    (5 min)
python run_sweep.py              # Section 5          -> data/sweep.pkl     (15 min)
python run_timing.py             # Section 2.3        -> data/timing.pkl    (1 min, run alone)
python make_figs.py              # all figures        -> figs/
python summarize.py              # prints the numbers and tables of this file
```

`config.py` holds the setting: the terminal scale is read from
`data/scale.json`, the set and URA scales default to 5, the environment
variables `SCALE`, `SET_SCALE`, `URA_SCALE`, `TAG` run a variant into
`data<TAG>/` and `figs<TAG>/`. The MPC class is `design/ra_mpc_cert.py`
(`nsub`, `safe_sub`, `safe_margin` and the fallback guesses are off by default
in the earlier scripts of `design/`).

## 10. Earlier settings (for the record)

| setting | terminal scale | URA | initial velocities | URA / vanilla / RA-MPC | where |
| --- | --- | --- | --- | --- | --- |
| A, 2026-10-05 | 0.159 | the terminal controller | around the field at 0.159, $r_0 = 0.03$ | 229 / 165 / 236 | `final/data_s0159`, `final/figs_s0159` |
| B, 2026-10-06 | 0.159 | $s = 5$ | around the field at 5, $r_0 = 0.94$ | 229 / 213 / 258 | `final/data_s5`, `final/figs_s5` |
| C, this document | 0.764 | $s = 5$ | as B | 229 / 213 / 263 | `final/data`, `final/figs` |

B differs from C only in the terminal scale and the initial-guess strategy.
The time-scaling derivation behind C is Section 7 of `TAC_EXPERIMENT_UPDATE.md`.
