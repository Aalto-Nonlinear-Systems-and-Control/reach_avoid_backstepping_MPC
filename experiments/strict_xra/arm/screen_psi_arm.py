"""Screen safe-set functions psi for the two-link arm: psi = w(y) * b(y) - c |y|^4 with b the band of the paper.
Checks: the safe set stays inside the workspace (|y| < 8) and away from the folded singularity, has one maximum inside
the target and no other critical point outside it, and how much of the band is certifiable (psi above the largest
critical value outside the target)."""
import numpy as np, sympy as sp
from scipy.optimize import root
y1, y2 = sp.symbols("y1 y2")
b = -((4 * (y1 - 2) - 2 * y2**3)**2) + sp.Rational(4, 5) * y2**3 + 10
phi = (y1 - sp.Rational(11, 2))**2 / sp.Rational(36, 25) + (y2 - sp.Rational(9, 5))**2 / sp.Rational(4, 25) - 2
r2 = y1**2 + y2**2
L = 8.0
W = {}
for c in (0.005, 0.01, 0.02):
    W[f"b - {c}|y|^4"] = (b - c * r2**2)
    for k in (0.3, 0.6, 1.0):
        W[f"(1+{k}(y2+2.5)) b - {c}|y|^4"] = ((1 + k * (y2 + sp.Rational(5, 2))) * b - c * r2**2)
g1 = np.linspace(-6, 9, 1501); g2 = np.linspace(-4, 4, 801); Y1, Y2 = np.meshgrid(g1, g2); R = np.hypot(Y1, Y2)
F = sp.lambdify((y1, y2), phi)(Y1, Y2); tgt = F < 0
bf = sp.lambdify((y1, y2), b); band = bf(Y1, Y2) > 0
for name, psi in W.items():
    pf = sp.lambdify((y1, y2), psi); P = pf(Y1, Y2); safe = P > 0
    if not safe.any(): print(name, "empty"); continue
    G = sp.lambdify((y1, y2), [sp.diff(psi, v) for v in (y1, y2)]); Hs = sp.lambdify((y1, y2), sp.hessian(psi, (y1, y2)))
    crit = []
    for a0 in np.linspace(-6, 9, 61):
        for b0 in np.linspace(-4, 4, 33):
            if pf(a0, b0) <= 0: continue
            sol = root(lambda z: np.array(G(z[0], z[1]), float), [a0, b0], tol=1e-13); z = sol.x
            if np.linalg.norm(np.array(G(*z), float)) < 1e-8 and pf(*z) > 1e-6 and not any(np.linalg.norm(z - c[0]) < 1e-4 for c in crit):
                ev = np.linalg.eigvalsh(np.array(Hs(*z), float))
                kind = "max" if ev[1] < -1e-7 else ("min" if ev[0] > 1e-7 else ("saddle" if ev[0] < -1e-7 and ev[1] > 1e-7 else "degenerate"))
                crit.append((z, float(pf(*z)), kind, float(sp.lambdify((y1, y2), phi)(*z)) > 0))
    pmax = P[safe].max(); cmax = max([c[1] for c in crit if c[3]], default=0.0)
    out = safe & ~tgt; cert = (P > cmax) & out
    c2 = (R**2 - 32) / 32; q2 = np.degrees(np.arccos(np.clip(c2, -1, 1)))
    print(f"\n{name}: safe set |y| in [{R[safe].min():.2f}, {R[safe].max():.2f}] (q2 in [{q2[safe].min():.0f}, {q2[safe].max():.0f}] deg, min sin q2 {np.abs(np.sin(np.radians(q2[safe]))).min():.2f}), "
          f"y1 in [{Y1[safe].min():.1f}, {Y1[safe].max():.1f}], y2 in [{Y2[safe].min():.1f}, {Y2[safe].max():.1f}], target part inside safe set {100*(safe & tgt).sum()/tgt.sum():.0f} %")
    for z, pv, k, o in sorted(crit, key=lambda c: -c[1]):
        print(f"    critical ({z[0]:+.3f}, {z[1]:+.3f})  psi/max = {pv/pmax:.3f}  {k:10s} {'OUTSIDE' if o else 'inside '} X_T")
    print(f"    certifiable share of the safe set outside X_T: {100*cert.sum()/out.sum():.0f} %  (bound delta_s/lambda >= {cmax/pmax:.3f} of max)")
