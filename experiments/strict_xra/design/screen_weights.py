"""Screen weights alpha for psi = alpha*h/max: positivity, critical points outside X_T, gradient profile."""
import numpy as np, sympy as sp
from scipy.optimize import root
from scipy.ndimage import minimum_filter
y1, y2 = sp.symbols("y1 y2")
phi = y2**2 + ((y1 + sp.Rational(17, 10)) / sp.Rational(1, 2))**2 - sp.Rational(2, 5)
h = -(y1**4 + y2**4 - 16) * (y1**4 + y2**4 - 4)
R = sp.Rational
W = {
    "W0  (2.5-y1)^2": (R(5, 2) - y1)**2,
    "W1a (2.5-y1)^2 + 0.5 y1 y2^2": (R(5, 2) - y1)**2 + R(1, 2) * y1 * y2**2,
    "W1b (2.5-y1)^2 + 1.0 y1 y2^2": (R(5, 2) - y1)**2 + y1 * y2**2,
    "W1c (2.5-y1)^2 + 1.5 y1 y2^2": (R(5, 2) - y1)**2 + R(3, 2) * y1 * y2**2,
    "W2b (3-y1)^2 + 1.0 y1 y2^2": (3 - y1)**2 + y1 * y2**2,
    "W3  (2.5-y1)^2 + 0.3 y1 y2^2 (2.5-y1)": (R(5, 2) - y1)**2 + R(3, 10) * y1 * y2**2 * (R(5, 2) - y1),
}
g = np.linspace(-2.05, 2.05, 1641); Y1, Y2 = np.meshgrid(g, g); dA = (g[1] - g[0])**2
H = sp.lambdify((y1, y2), h)(Y1, Y2); F = sp.lambdify((y1, y2), phi)(Y1, Y2)
safe = H > 0; ring_out = safe & (F > 0)
for name, a in W.items():
    A = sp.lambdify((y1, y2), a)(Y1, Y2)
    psi = a * h; pmax = float((A * H)[safe].max()); psi = psi / pmax
    P = sp.lambdify((y1, y2), psi)(Y1, Y2)
    gr = [sp.diff(psi, v) for v in (y1, y2)]; G = sp.lambdify((y1, y2), gr); Hs = sp.lambdify((y1, y2), sp.hessian(psi, (y1, y2)))
    g1, g2 = G(Y1, Y2); N = np.hypot(g1, g2)
    Nm = np.where(safe, N, np.inf)
    crit = []
    hf = sp.lambdify((y1, y2), h); pf = sp.lambdify((y1, y2), psi); ff = sp.lambdify((y1, y2), phi)
    gs = np.linspace(-2.0, 2.0, 81)
    for a0 in gs:
        for b0 in gs:
            if hf(a0, b0) <= 0: continue
            sol = root(lambda z: np.array(G(z[0], z[1]), float), [a0, b0], tol=1e-14)
            z = sol.x
            if np.linalg.norm(np.array(G(*z), float)) < 1e-9 and hf(*z) > 1e-6:
                if not any(np.linalg.norm(z - c[0]) < 1e-4 for c in crit):
                    ev = np.linalg.eigvalsh(np.array(Hs(*z), float))
                    kind = "max" if ev[1] < -1e-7 else ("min" if ev[0] > 1e-7 else ("saddle" if ev[0] < -1e-7 and ev[1] > 1e-7 else "degenerate"))
                    crit.append((z, float(pf(*z)), kind, float(ff(*z)) > 0))
    cmax = max([c[1] for c in crit if c[3]], default=0.0)
    print(f"\n{name}:  min alpha on safe set = {A[safe].min():.3f}, alpha(target)/alpha(far) = {float(a.subs({y1:-1.78,y2:0})/a.subs({y1:1.7,y2:0})):.1f}, deg psi = {sp.Poly(sp.expand(psi), y1, y2).total_degree()}")
    for z, pv, k, o in sorted(crit, key=lambda c: -c[1]):
        print(f"    critical ({z[0]:+.3f}, {z[1]:+.3f})  psi/max = {pv:.4f}  {k:10s} {'OUTSIDE' if o else 'inside '} X_T")
    print(f"    bound delta_s/lambda >= {cmax:.4f};  certifiable share of ring outside X_T: {((P > cmax) & (F > 0)).sum() / ring_out.sum() * 100:.0f}%")
    row = []
    for c in (0.05, 0.1, 0.2, 0.3, 0.5, 0.7):
        m = (P >= c) & (F > 0)
        row.append(f"c={c}: share {m.sum()/ring_out.sum()*100:3.0f}%, min|grad| {N[m].min():.3f}")
    print("    " + " | ".join(row))
