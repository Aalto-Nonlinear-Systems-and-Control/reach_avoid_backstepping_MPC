"""Uniform sampling period vs terminal level delta: second-order bound, simulation, invariance."""
import importlib.util, numpy as np, sympy as sp
s = importlib.util.spec_from_file_location("m", "out/rigorous_alpha2_mu0.1_vanilla_k1.py"); m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
x1, x2, th, v, w, a = sp.symbols("x1 x2 th v w a"); X = [x1, x2, th, v]
Psi = sp.sympify(m.certificate_opt); U = [sp.sympify(e) for e in m.u_opt]
F = sp.Matrix([v * sp.cos(th), v * sp.sin(th), w, a])
Pd = (sp.Matrix([sp.diff(Psi, z) for z in X]).T * F)[0]; Pdd = (sp.Matrix([sp.diff(Pd, z) for z in X]).T * F)[0]
Pf = sp.lambdify(X, Psi, "numpy"); Uf = sp.lambdify(X, U, "numpy"); Pddf = sp.lambdify(X + [w, a], Pdd, "numpy")
Pdc = sp.lambdify(X, Pd.subs({w: U[0], a: U[1]}), "numpy")
r = np.random.default_rng(5); lo = np.array([-2.2, -2.2, 2*np.pi/3, 0.1]); hi = np.array([2.2, 2.2, 4*np.pi/3, 1.0])
S = lo[:, None] + (hi - lo)[:, None] * r.random((4, 1500000)); P = np.asarray(Pf(*S), float)
tg = S[1]**2 + ((S[0] + 1.7) / 0.5)**2 - 0.4; pmax = P.max()
base = (P >= 0) & (tg > 0); X0 = S[:, base][:, :40000]; P0 = P[base][:40000]
U0 = np.array(Uf(*X0), float); D0 = np.asarray(Pdc(*X0), float)

def arc(T, n=50, need_M2=True):
    h = T / n; x = X0.copy(); M = np.abs(Pddf(*x, *U0)) if need_M2 else None
    rhs = lambda z: np.array([z[3]*np.cos(z[2]), z[3]*np.sin(z[2]), U0[0], U0[1]])
    for _ in range(n):
        q1 = rhs(x); q2 = rhs(x + h/2*q1); q3 = rhs(x + h/2*q2); q4 = rhs(x + h*q3); x = x + h/6*(q1 + 2*q2 + 2*q3 + q4)
        if need_M2: M = np.maximum(M, np.abs(Pddf(*x, *U0)))
    return np.asarray(Pf(*x), float), M

Tgrid = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05, 0.1]
P1 = {T: arc(T, need_M2=False)[0] for T in Tgrid}
_, M2loc = arc(0.05)
Tc = 2 * D0 / M2loc                      # pointwise second-order period (local M2)
# where are the slowest states?
slow = np.argsort(Tc)[:200]
print(f"max Psi {pmax:.3f}; {X0.shape[1]} states in X_RA\\X_T")
print(f"slowest 200 states (smallest 2*Psidot/M2): median position ({np.median(X0[0,slow]):.2f}, {np.median(X0[1,slow]):.2f}), "
      f"median Psi {np.median(P0[slow]):.3f} ({np.median(P0[slow])/pmax*100:.0f}% of max), median Psidot {np.median(D0[slow]):.4f}")
print(f"pointwise period 2*Psidot/M2_local over all states: 1% {np.quantile(Tc,0.01):.2e}, 10% {np.quantile(Tc,0.1):.2e}, median {np.median(Tc):.2e} s")
print("\n delta/maxPsi | share of states | bound 2*min(Psidot)/max(M2) | min_x 2*Psidot/M2_local | largest T with Psi increasing at ALL states | largest T keeping {Psi>=delta} invariant")
for frac in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7):
    d = frac * pmax; sel = P0 >= d
    if sel.sum() < 50: continue
    tb = 2 * D0[sel].min() / M2loc[sel].max(); tp = Tc[sel].min()
    inc = max([T for T in Tgrid if np.all(P1[T][sel] > P0[sel])], default=0)
    inv = max([T for T in Tgrid if np.all(P1[T][sel] >= d)], default=0)
    print(f"   {frac:4.2f}       | {sel.mean()*100:5.1f}%          | {tb:.2e} s                  | {tp:.2e} s              | {inc:g} s                                    | {inv:g} s")
