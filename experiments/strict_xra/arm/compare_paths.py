"""Four corridor options for the top-view arm, drawn with the nominal path, the tolerance corridor, the placement disc and the arm at the start."""
import sys, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from screen_paths import build, base, analyse, path_metrics, y1, y2, Y1, Y2, RT
from cert_arm import ik

OPTIONS = [
    ("A  two bends, 190 deg (current)", dict(harms=[(3, 0.25, 0)], norm="none")),
    ("B  wide bend + tight bend, 190 deg", dict(harms=[(3, 0.25, 0), (2, 0.2, 135)])),
    ("C  wrap-around, 210 deg", dict(harms=[(3, 0.25, 30)], weight="long3", c=4.0, e=(0, 1), d=(0, 1), kappa=0.2, eps=0.0002, m=3)),
    ("D  wrap-around, 250 deg (long tail)", dict(harms=[(3, 0.25, 30)], weight="long3", c=4.0, e=(0, 1), d=(0, 1), kappa=0.2, eps=0.002, m=2)),
]

if __name__ == "__main__":
    fig, axs = plt.subplots(2, 2, figsize=(11, 10.5)); axs = axs.ravel()
    for ax, (name, kw) in zip(axs, OPTIONS):
        cfg = dict(base); cfg.update(kw); p, psi = build(**cfg); r = analyse(psi); pm = path_metrics(p, psi, r, cfg["rho0"]); pp = sp.lambdify((y1, y2), p)
        ax.contourf(Y1, Y2, r["P"] / r["pmax"], levels=[0, 0.5, 10], colors=["#dfe9f5", "#b9d3ee"]); ax.contour(Y1, Y2, r["P"], levels=[0], colors="#2166ac", linewidths=1.8)
        ax.contour(Y1, Y2, np.where(r["P"] > 0, pp(Y1, Y2), np.nan), levels=[0], colors="0.3", linestyles="--", linewidths=1.1)
        ax.add_patch(plt.Circle(r["tc"], RT, color="#fee0d2", ec="#d6604d", lw=1.8, zorder=5))
        th = np.linspace(0, 2 * np.pi, 400); ax.plot(8 * np.cos(th), 8 * np.sin(th), "k-", lw=1)
        for rr in (2.5, 7): ax.plot(rr * np.cos(th), rr * np.sin(th), "k:", lw=0.7)
        X = pm["X"]; ends = (X[0], X[-1]); y0 = max(ends, key=lambda x: np.hypot(*(x - r["tc"]))); pf = sp.lambdify((y1, y2), psi)
        # move the start a little inward along the path so that it is well inside the corridor
        k = 0 if np.allclose(y0, X[0]) else len(X) - 1; step = 1 if k == 0 else -1; j = k
        while pf(*X[j]) < 0.15 * r["pmax"] and 0 <= j + step < len(X): j += step
        y0 = X[j]; q = ik(y0.reshape(2, 1))[:, 0]; jx, jy = 4 * np.cos(q[0]), 4 * np.sin(q[0])
        ax.plot([0, jx], [0, jy], "o-", lw=6, color="0.5", alpha=0.6, zorder=7); ax.plot([jx, y0[0]], [jy, y0[1]], "o-", lw=6, color="0.5", alpha=0.6, zorder=7); ax.plot(y0[0], y0[1], "o", color="#ff8c00", ms=8, zorder=8)
        ax.plot(0, 0, "ks", ms=6); ax.set_aspect("equal"); ax.set_xlim(-8.3, 8.3); ax.set_ylim(-8.3, 8.3); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{name}\n|y| in [{r['rmin']:.1f}, {r['rmax']:.1f}], degree {r['deg']}, path length {pm['length']:.0f}, width {pm['wmed']:.1f} (median), target {min(pm['to_ends']):.1f} before the end, terminal share {100*r['half']:.0f} %", fontsize=8.5)
        print(name, "| start", np.round(y0, 2), "| target", np.round(r["tc"], 2), "| deg", r["deg"], "| span", round(pm["span"]), "| cert", round(100 * r["cert"]), flush=True)
    fig.tight_layout(); fig.savefig("figs/arm_path_options.png", dpi=130); print("saved figs/arm_path_options.png")
