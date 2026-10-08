"""Workcell layout: base at the bench edge, corridor on the bench in front. Options drawn with bench, nominal path, corridor, fixtures at the bends, placement disc and the arm at the start."""
import sys, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from screen_front import build_front, build_lateral, metrics_lateral, Y_EDGE
from screen_paths import path_metrics
from screen_zigzag import analyse, y1, y2, Y1, Y2, RT
from cert_arm import ik

OPTIONS = [
    ("E1  front arc, two gentle bends", "front", dict(rho0=4.5, harms=[(3, 0.25, 60)], h=1.0, kappa=0.05, eps=0.005, m=3)),
    ("E2  front arc, one straight stretch", "front", dict(rho0=4.5, harms=[(3, 0.25, 90)], h=1.0, kappa=0.05, eps=0.005, m=3)),
    ("E3  front arc, detour around a central fixture", "front", dict(rho0=4.5, harms=[(4, 0.15, 0)], h=1.0, kappa=0.05, eps=0.005, m=3)),
    ("E4  front arc, uneven bends", "front", dict(rho0=4.5, harms=[(3, 0.25, 60), (2, 0.12, 0)], h=1.0, kappa=0.05, eps=0.04, m=2)),
]


def fixtures(pf, X, k_bends=2):
    """Grey fixtures beside the corridor at the interior points of largest path curvature, on the concave side."""
    xs, ys = X[:, 0], X[:, 1]; dx, dy = np.gradient(xs), np.gradient(ys); ddx, ddy = np.gradient(dx), np.gradient(dy); curv = (dx * ddy - dy * ddx) / (dx**2 + dy**2)**1.5
    idx = []; order = np.argsort(-abs(curv)); lo, hi = int(0.18 * len(X)), int(0.82 * len(X))
    for i in order:
        if lo < i < hi and all(np.hypot(*(X[i] - X[j])) > 2.5 for j in idx): idx.append(i)
        if len(idx) == k_bends: break
    out = []
    for i in idx:
        n = np.array([-dy[i], dx[i]]); n /= np.linalg.norm(n); n *= np.sign(curv[i])   # towards the centre of curvature
        t = 0.0
        while pf(*(X[i] + t * n)) > 0 and t < 4: t += 0.02
        rf = 0.7; out.append((X[i] + (t + rf + 0.2) * n, rf))
    return out


if __name__ == "__main__":
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.6)); axs = axs.ravel()
    for ax, (name, kind, kw) in zip(axs, OPTIONS):
        if kind == "front": p, psi = build_front(**kw); r = analyse(psi); pm = path_metrics(p, psi, r, kw["rho0"])
        else: p, psi = build_lateral(**kw); r = analyse(psi); pm = metrics_lateral(p, psi, r)
        pf = sp.lambdify((y1, y2), psi); pp = sp.lambdify((y1, y2), p); safe = r["P"] > 0
        ax.add_patch(plt.Rectangle((-6.5, Y_EDGE), 13, 7.0, color="#f5f0e6", zorder=0))
        ax.contourf(Y1, Y2, r["P"] / r["pmax"], levels=[0, 0.5, 10], colors=["#dfe9f5", "#b9d3ee"]); ax.contour(Y1, Y2, r["P"], levels=[0], colors="#2166ac", linewidths=1.8)
        ax.contour(Y1, Y2, np.where(safe, pp(Y1, Y2), np.nan), levels=[0], colors="0.3", linestyles="--", linewidths=1.1)
        ax.add_patch(plt.Circle(r["tc"], RT, color="#fee0d2", ec="#d6604d", lw=1.8, zorder=5))
        X = pm["X"]
        for (cf, rf) in fixtures(pf, X): ax.add_patch(plt.Circle(cf, rf, color="0.6", ec="k", zorder=6))
        th = np.linspace(0, np.pi, 200); ax.plot(8 * np.cos(th), 8 * np.sin(th), "k-", lw=1)
        for rr in (2.5, 7): ax.plot(rr * np.cos(th), rr * np.sin(th), "k:", lw=0.7)
        ends = (X[0], X[-1]); y0 = max(ends, key=lambda x: np.hypot(*(x - r["tc"]))); k = 0 if np.allclose(y0, X[0]) else len(X) - 1; step = 1 if k == 0 else -1; j = k
        while pf(*X[j]) < 0.15 * r["pmax"] and 0 <= j + step < len(X): j += step
        y0 = X[j]; q = ik(y0.reshape(2, 1))[:, 0]; jx, jy = 4 * np.cos(q[0]), 4 * np.sin(q[0])
        # draw the elbow-up branch (links stay above the bench edge in the top view): mirror the elbow across the base-to-hand line
        d = y0 / np.linalg.norm(y0); e = np.array([jx, jy]); e_up = 2 * d * (d @ e) - e
        if e_up[1] > e[1]: jx, jy = e_up
        ax.plot([0, jx], [0, jy], "o-", lw=6, color="0.5", alpha=0.6, zorder=7); ax.plot([jx, y0[0]], [jy, y0[1]], "o-", lw=6, color="0.5", alpha=0.6, zorder=7); ax.plot(y0[0], y0[1], "o", color="#ff8c00", ms=8, zorder=8)
        ax.plot(0, 0, "ks", ms=7); ax.set_aspect("equal"); ax.set_xlim(-7.2, 7.2); ax.set_ylim(-1.2, 8.2); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{name}\n|y| in [{r['rmin']:.1f}, {r['rmax']:.1f}], y2 >= {Y2[safe].min():.1f}, degree {r['deg']}, length {pm['length']:.0f}, width {pm['wmed']:.1f}, target {min(pm['to_ends']):.1f} before the end", fontsize=8.5)
        print(name, "| start", np.round(y0, 2), "| target", np.round(r["tc"], 2), "| deg", r["deg"], "| others", r["n_other"], "| cert", round(100 * r["cert"]), "| y2min", round(Y2[safe].min(), 2), flush=True)
    fig.tight_layout(); fig.savefig("figs/arm_front_options.png", dpi=130); print("saved figs/arm_front_options.png")
