"""Candidate corridors for the top-view arm, closed along the corridor instead of by |y|^4:
psi = w(y) [h^2 - (y1 - a - b (y2 - y0)^2)^2] - eps (y2 - y0)^4,  w = 1 + kappa (y2 + 3).
The target ellipse (semi-axes 1.70, 0.57 as in the paper) is centred at the maximum of psi."""
import sys, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import root
y1, y2 = sp.symbols("y1 y2")
g1 = np.linspace(-6, 9, 751); g2 = np.linspace(-4.5, 4.5, 451); Y1, Y2 = np.meshgrid(g1, g2); R = np.hypot(Y1, Y2)


def analyse(a, b, y0, h, k, eps):
    psi = (1 + sp.nsimplify(k) * (y2 + 3)) * (sp.nsimplify(h)**2 - (y1 - sp.nsimplify(a) - sp.nsimplify(b) * (y2 - sp.nsimplify(y0))**2)**2) - sp.nsimplify(eps) * (y2 - sp.nsimplify(y0))**4
    pf = sp.lambdify((y1, y2), psi); P = pf(Y1, Y2); safe = P > 0
    if not safe.any(): return None
    G = sp.lambdify((y1, y2), [sp.diff(psi, v) for v in (y1, y2)]); Hs = sp.lambdify((y1, y2), sp.hessian(psi, (y1, y2))); crit = []
    for a0 in np.linspace(-6, 9, 31):
        for b0 in np.linspace(-4.5, 4.5, 19):
            if pf(a0, b0) <= 0: continue
            sol = root(lambda z: np.array(G(z[0], z[1]), float), [a0, b0], tol=1e-13); z = sol.x
            if np.linalg.norm(np.array(G(*z), float)) < 1e-8 and pf(*z) > 1e-6 and not any(np.linalg.norm(z - cc[0]) < 1e-4 for cc in crit):
                ev = np.linalg.eigvalsh(np.array(Hs(*z), float)); crit.append((z, float(pf(*z)), "max" if ev[1] < -1e-7 else ("min" if ev[0] > 1e-7 else "saddle")))
    pmax = P[safe].max(); mx = max(crit, key=lambda cc: cc[1]); tgt_c = mx[0]
    phi = (y1 - tgt_c[0])**2 / sp.Rational(36, 25) + (y2 - tgt_c[1])**2 / sp.Rational(4, 25) - 2; ff = sp.lambdify((y1, y2), phi); F = ff(Y1, Y2); tgt = F < 0
    others = [cc for cc in crit if ff(*cc[0]) > 0]; cmax = max([cc[1] for cc in others], default=0.0); out = safe & ~tgt
    return dict(psi=psi, phi=phi, P=P, F=F, pmax=pmax, crit=crit, tgt_c=tgt_c, rmin=R[safe].min(), rmax=R[safe].max(), y2min=Y2[safe].min(), y2max=Y2[safe].max(),
                cert=((P > cmax) & out).sum() / out.sum(), half=((P >= 0.5 * pmax) & out).sum() / out.sum(), tin=(safe & tgt).sum() / tgt.sum(), n_other=len(others))


if __name__ == "__main__":
    grid = []
    for a in (3.5, 4.0):
        for b in (0.4, 0.6):
            for y0 in (0.5, 1.0):
                for h in (1.0, 1.2):
                    for k in (1.0, 2.0):
                        for eps in (0.1, 0.2, 0.4):
                            r = analyse(a, b, y0, h, k, eps)
                            if r is None: continue
                            ok = r["n_other"] == 0 and r["rmin"] >= 2.4 and r["rmax"] <= 7.1 and r["cert"] >= 0.95
                            grid.append((ok, a, b, y0, h, k, eps, r))
                            print(f"a={a} b={b} y0={y0} h={h} k={k} eps={eps}: |y| in [{r['rmin']:.2f}, {r['rmax']:.2f}], y2 in [{r['y2min']:.2f}, {r['y2max']:.2f}], max at ({r['tgt_c'][0]:.2f}, {r['tgt_c'][1]:.2f}), other critical {r['n_other']}, certifiable {100*r['cert']:.0f} %, psi>=0.5max {100*r['half']:.0f} %, target in set {100*r['tin']:.0f} % {'<-- ok' if ok else ''}", flush=True)
    good = [g for g in grid if g[0]]; good.sort(key=lambda g: (-(g[7]["y2max"] - g[7]["y2min"]), -g[7]["half"]))
    sel = good[:6] if len(good) >= 6 else good
    fig, axs = plt.subplots(2, 3, figsize=(18, 9.5), constrained_layout=True)
    for ax, (ok, a, b, y0, h, k, eps, r) in zip(axs.ravel(), sel):
        ax.contourf(Y1, Y2, r["P"] / r["pmax"], levels=[0, 0.5, 10], colors=["#dfe9f5", "#b9d3ee"]); ax.contour(Y1, Y2, r["P"], levels=[0], colors="#2166ac", linewidths=1.5)
        ax.contourf(Y1, Y2, -r["F"], levels=[0, 1e9], colors=["#fee0d2"], alpha=0.9); ax.contour(Y1, Y2, r["F"], levels=[0], colors="#d6604d", linewidths=1.5)
        th = np.linspace(0, 2 * np.pi, 400)
        for rr, ls in ((8, "k-"), (7, "k:"), (2.5, "k:")): ax.plot(rr * np.cos(th), rr * np.sin(th), ls, lw=0.9)
        ax.plot(0, 0, "ko"); [ax.plot(cc[0][0], cc[0][1], "r*" if cc[2] == "max" else "rx", ms=10) for cc in r["crit"]]
        ax.set_aspect("equal"); ax.set_xlim(-6, 9); ax.set_ylim(-4.5, 4.5); ax.set_xlabel("$y_1$"); ax.set_ylabel("$y_2$")
        ax.set_title(f"a={a}, b={b}, y0={y0}, h={h}, kappa={k}, eps={eps}\n|y| in [{r['rmin']:.1f}, {r['rmax']:.1f}], length y2 in [{r['y2min']:.1f}, {r['y2max']:.1f}], certifiable {100*r['cert']:.0f} %, psi>=0.5max {100*r['half']:.0f} %, target centre ({r['tgt_c'][0]:.2f}, {r['tgt_c'][1]:.2f})", fontsize=8.5)
    fig.suptitle("Top-view arm: corridor candidates satisfying the requirements (unique maximum inside the target, no other critical point, |y| in [2.4, 7.1]); blue: safe set, darker: psi >= 0.5 max; red: target; |y| = 8 solid, 2.5 and 7 dotted", fontsize=10)
    fig.savefig("figs/arm_corridors.png", dpi=130); print(f"{len(good)} candidates ok; saved figs/arm_corridors.png")
