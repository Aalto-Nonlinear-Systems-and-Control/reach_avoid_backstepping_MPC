"""Wavy arc around the base with an angular modulation of the radius (zig-zag of uniform radial width):
psi = w [h^2 - (|y|^2 - rho0^2 - A C_n(y) / rho0^(n-2))^2 / (2 rho0)^2] - eps (y1 - rho0)^4,
C_n = Re((y1 + i y2)^n) or Im(...) (phase), w = 1 + kappa (y2 + 6)^2, circular target of radius 1 at the maximum."""
import sys, itertools, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "."); from screen_zigzag import analyse, draw, y1, y2, r2
_a, _b = sp.symbols("a b", real=True)
def Cn(n, phase):
    e = sp.expand((_a + sp.I * _b)**n); c = sp.re(e) if phase == "cos" else sp.im(e)
    return sp.expand(c).subs({_a: y1, _b: y2}, simultaneous=True)
out = []
for rho0, n, ph, A, h, k, eps in itertools.product((4.0, 4.25), (3, 4), ("cos", "sin"), (0.25, 0.3, 0.35), (1.0,), (0.1, 0.3), (0.01, 0.03)):
    mod = A * Cn(n, ph) / rho0**(n - 2)
    psi = (1 + sp.nsimplify(k) * (y2 + 6)**2) * (sp.nsimplify(h)**2 - (r2 - rho0**2 - mod)**2 / (2 * rho0)**2) - sp.nsimplify(eps) * (y1 - rho0)**4
    r = analyse(psi)
    if r is None: continue
    ok = r["n_other"] == 0 and r["rmin"] >= 2.4 and r["rmax"] <= 7.1 and r["cert"] >= 0.95
    name = f"rho0={rho0}, n={n} {ph}, A={A}, h={h}, kappa={k}, eps={eps}"; out.append((name, r, ok))
    print(f"{name}: |y| [{r['rmin']:.2f}, {r['rmax']:.2f}], max ({r['tc'][0]:.2f}, {r['tc'][1]:.2f}), others {r['n_other']}, cert {100*r['cert']:.0f} %, half {100*r['half']:.0f} %, degree {r['deg']}, area {r['area']:.1f} {'<-- ok' if ok else ''}", flush=True)
good = [x for x in out if x[2]]; good.sort(key=lambda x: -x[1]["area"]); sel = good[:: max(1, len(good) // 6)][:6]
fig, axs = plt.subplots(2, 3, figsize=(18, 10.5), constrained_layout=True)
for ax, (name, r, ok) in zip(axs.ravel(), sel): draw(ax, name, r)
fig.suptitle("Top-view arm, wavy arcs around the base (radius modulated by cos/sin n theta), circular target (radius 1) at the maximum; blue safe set, darker psi >= 0.5 max; |y| = 8 solid, 2.5 and 7 dotted", fontsize=10)
fig.savefig("figs/arm_wavy.png", dpi=130); print(f"{len(good)} admissible; saved figs/arm_wavy.png")
