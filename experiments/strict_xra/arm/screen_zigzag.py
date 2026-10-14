"""Zig-zag safe sets of (nearly) uniform width for the top-view arm, circular target at the maximum of psi.
Z1 wavy arc:   psi = w [h^2 - (|y|^2 - rho(y2)^2)^2 / (2 rho0)^2] - eps (y1 - rho0)^4,  rho(y2) = rho0 + A s(y2)
Z2 S corridor: psi = w [h^2 (1 + c'(y2)^2) - (y1 - c(y2))^2] - eps (y2 - yc)^4,     c(y2) = a + A s(y2)
with s(y2) = (3 sqrt3 / 2) (y2 / L) (1 - (y2 / L)^2) (amplitude 1 at y2 = +-L/sqrt3) and w = 1 + kappa (y2 + 6)^2."""
import sys, itertools, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import root
y1, y2 = sp.symbols("y1 y2"); r2 = y1**2 + y2**2
g = np.linspace(-7.5, 8.5, 801); gy = np.linspace(-7, 7, 701); Y1, Y2 = np.meshgrid(g, gy); R = np.hypot(Y1, Y2)
RT = 1.0   # target radius


def shape(L, sign=1):
    return sign * sp.Rational(3, 2) * sp.sqrt(3) * (y2 / L) * (1 - (y2 / L)**2)


def analyse(psi):
    pf = sp.lambdify((y1, y2), psi); P = pf(Y1, Y2); safe = P > 0
    if not safe.any(): return None
    G = sp.lambdify((y1, y2), [sp.diff(psi, v) for v in (y1, y2)]); Hs = sp.lambdify((y1, y2), sp.hessian(psi, (y1, y2))); crit = []
    for a0 in np.linspace(-7, 8, 31):
        for b0 in np.linspace(-7, 7, 29):
            if pf(a0, b0) <= 0: continue
            sol = root(lambda z: np.array(G(z[0], z[1]), float), [a0, b0], tol=1e-13); z = sol.x
            if np.linalg.norm(np.array(G(*z), float)) < 1e-8 and pf(*z) > 1e-6 and not any(np.linalg.norm(z - cc[0]) < 1e-4 for cc in crit):
                ev = np.linalg.eigvalsh(np.array(Hs(*z), float)); crit.append((z, float(pf(*z)), "max" if ev[1] < -1e-7 else ("min" if ev[0] > 1e-7 else "saddle")))
    if not crit: return None
    pmax = P[safe].max(); tc = max(crit, key=lambda cc: cc[1])[0]
    F = (Y1 - tc[0])**2 + (Y2 - tc[1])**2 - RT**2; tgt = F < 0; out = safe & ~tgt
    others = [cc for cc in crit if (cc[0][0] - tc[0])**2 + (cc[0][1] - tc[1])**2 > RT**2]; cmax = max([cc[1] for cc in others], default=0.0)
    return dict(P=P, F=F, pmax=pmax, crit=crit, tc=tc, rmin=R[safe].min(), rmax=R[safe].max(), cert=((P > cmax) & out).sum() / out.sum(), half=((P >= 0.5 * pmax) & out).sum() / out.sum(),
                n_other=len(others), area=safe.sum() * (g[1] - g[0]) * (gy[1] - gy[0]), deg=sp.Poly(sp.expand(psi), y1, y2).total_degree())


def draw(ax, name, r):
    ax.contourf(Y1, Y2, r["P"] / r["pmax"], levels=[0, 0.5, 10], colors=["#dfe9f5", "#b9d3ee"]); ax.contour(Y1, Y2, r["P"], levels=[0], colors="#2166ac", linewidths=1.5)
    ax.contourf(Y1, Y2, -r["F"], levels=[0, 1e9], colors=["#fee0d2"], alpha=0.9); ax.contour(Y1, Y2, r["F"], levels=[0], colors="#d6604d", linewidths=1.5)
    th = np.linspace(0, 2 * np.pi, 400)
    for rr, ls in ((8, "k-"), (7, "k:"), (2.5, "k:")): ax.plot(rr * np.cos(th), rr * np.sin(th), ls, lw=0.9)
    ax.plot(0, 0, "ko"); [ax.plot(cc[0][0], cc[0][1], "r*" if cc[2] == "max" else "rx", ms=10) for cc in r["crit"]]
    ax.set_aspect("equal"); ax.set_xlim(-7.5, 8.5); ax.set_ylim(-7, 7); ax.set_xlabel("$y_1$"); ax.set_ylabel("$y_2$")
    ax.set_title(f"{name}\n|y| in [{r['rmin']:.1f}, {r['rmax']:.1f}], degree {r['deg']}, other critical points {r['n_other']}, certifiable {100*r['cert']:.0f} %, psi>=0.5max {100*r['half']:.0f} %, target centre ({r['tc'][0]:.2f}, {r['tc'][1]:.2f})", fontsize=8.5)


if __name__ == "__main__":
    out = []
    for rho0, A, L, h, k, eps, sg in itertools.product((4.75,), (0.6, 1.0), (4.0, 5.0), (1.0, 1.2), (0.1,), (0.01, 0.03, 0.1), (1, -1)):
        rho = rho0 + A * shape(L, sg); psi = (1 + sp.nsimplify(k) * (y2 + 6)**2) * (sp.nsimplify(h)**2 - (r2 - rho**2)**2 / (2 * rho0)**2) - sp.nsimplify(eps) * (y1 - rho0)**4
        r = analyse(psi)
        if r is None: continue
        ok = r["n_other"] == 0 and r["rmin"] >= 2.4 and r["rmax"] <= 7.1 and r["cert"] >= 0.95
        name = f"Z1 wavy arc: rho0={rho0}, A={A}{'' if sg == 1 else ' (flipped)'}, L={L}, h={h}, eps={eps}"; out.append((name, r, ok))
        print(f"{name}: |y| [{r['rmin']:.2f}, {r['rmax']:.2f}], max ({r['tc'][0]:.2f}, {r['tc'][1]:.2f}), others {r['n_other']}, cert {100*r['cert']:.0f} %, half {100*r['half']:.0f} %, degree {r['deg']} {'<-- ok' if ok else ''}", flush=True)
    for a, A, L, h, k, eps, yc in itertools.product((4.5,), (1.0, 1.5), (4.0, 5.0), (1.0, 1.2), (0.1,), (0.02, 0.05, 0.1), (0.0,)):
        c = a + A * shape(L); cp = sp.diff(c, y2); psi = (1 + sp.nsimplify(k) * (y2 + 6)**2) * (sp.nsimplify(h)**2 * (1 + cp**2) - (y1 - c)**2) - sp.nsimplify(eps) * (y2 - yc)**4
        r = analyse(psi)
        if r is None: continue
        ok = r["n_other"] == 0 and r["rmin"] >= 2.4 and r["rmax"] <= 7.1 and r["cert"] >= 0.95
        name = f"Z2 S corridor: a={a}, A={A}, L={L}, h={h}, eps={eps}"; out.append((name, r, ok))
        print(f"{name}: |y| [{r['rmin']:.2f}, {r['rmax']:.2f}], max ({r['tc'][0]:.2f}, {r['tc'][1]:.2f}), others {r['n_other']}, cert {100*r['cert']:.0f} %, half {100*r['half']:.0f} %, degree {r['deg']} {'<-- ok' if ok else ''}", flush=True)
    good = [x for x in out if x[2]]; z1 = [x for x in good if x[0].startswith("Z1")]; z2 = [x for x in good if x[0].startswith("Z2")]
    sel = (z1[:: max(1, len(z1) // 3)][:3] if z1 else []) + (z2[:: max(1, len(z2) // 3)][:3] if z2 else [])
    fig, axs = plt.subplots(2, 3, figsize=(18, 10.5), constrained_layout=True)
    for ax, (name, r, ok) in zip(axs.ravel(), sel): draw(ax, name, r)
    fig.suptitle("Top-view arm, zig-zag safe sets with a circular target (radius 1) at the maximum of psi; blue safe set, darker psi >= 0.5 max; |y| = 8 solid, 2.5 and 7 dotted", fontsize=10)
    fig.savefig("figs/arm_zigzag.png", dpi=130); print(f"{len(z1)} Z1 and {len(z2)} Z2 admissible; saved figs/arm_zigzag.png")
