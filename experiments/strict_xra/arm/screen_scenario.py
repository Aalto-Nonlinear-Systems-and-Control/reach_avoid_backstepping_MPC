"""Pick-and-place scenario: object picked at A (right) and placed at B (left), A and B mirror images, the straight line A-B
blocked by obstacles so that the reference path bends. The safe set is the tolerance corridor around the reference path
(front-arc family of screen_front.build_front), the target the placement disc at B = argmax psi, A = mirror of B.
Checks beyond the usual ones: A well inside the corridor, the corridor extends beyond A, every obstacle outside the corridor with a margin."""
import sys, itertools, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from screen_front import build_front, Y_EDGE
from screen_paths import path_metrics
from screen_zigzag import analyse, y1, y2, Y1, Y2, RT
from cert_arm import ik
R_OBS = 0.75; MARGIN = 0.25


def obstacles_for(pm, r, pf):
    """Obstacles that block the straight line A-B: one at the midpoint of A-B; further ones on the concave side of
    additional bends if the path has them. Returns [(centre, radius)] and the clearance of each to the corridor."""
    B = np.array(r["tc"]); A = np.array([-B[0], B[1]]); obs = [((A + B) / 2, R_OBS)]
    safe = r["P"] > 0; pts = np.c_[Y1[safe], Y2[safe]]
    clear = [np.min(np.hypot(pts[:, 0] - c[0], pts[:, 1] - c[1])) - rr for c, rr in obs]
    return A, B, obs, clear


def scenario_checks(p, psi, r, pm):
    pf = sp.lambdify((y1, y2), psi); A, B, obs, clear = obstacles_for(pm, r, pf); X = pm["X"]
    thA = np.degrees(np.arctan2(A[1], A[0])); ths = np.degrees(np.arctan2(X[:, 1], X[:, 0])); beyond = thA - ths.min()
    psiA = pf(*A) / r["pmax"]; line_blocked = abs(obs[0][0][1] - A[1]) < R_OBS   # the segment A-B is horizontal at y2 = A[1]
    return dict(A=A, B=B, obs=obs, clear=clear, beyond=beyond, psiA=psiA, blocked=line_blocked)


CANDS = []
def add(name, **kw): CANDS.append((name, kw))
base = dict(rho0=4.5, h=1.0, kappa=0.05, eps=0.005, m=3)
for dcut, kappa in itertools.product((90, 82, 75), (0.03, 0.05)):
    add(f"bridge W3 A3=0.25 kappa={kappa} dcut={dcut}", harms=[(3, 0.25, 270)], kappa=kappa, dcut=dcut)
    add(f"skewed bridge W3 A3=0.25 phi3=240 kappa={kappa} dcut={dcut}", harms=[(3, 0.25, 240)], kappa=kappa, dcut=dcut)
    add(f"W-shape W4 A4=0.15 kappa={kappa} dcut={dcut}", harms=[(4, 0.15, 0)], kappa=kappa, dcut=dcut)
    add(f"bridge W3 A3=0.2 + W4 A4=0.06 kappa={kappa} dcut={dcut}", harms=[(3, 0.2, 270), (4, 0.06, 0)], kappa=kappa, dcut=dcut)
    add(f"bridge W3 A3=0.25 + W2 A2=0.08 phi2=90 kappa={kappa} dcut={dcut}", harms=[(3, 0.25, 270), (2, 0.08, 90)], kappa=kappa, dcut=dcut)


def draw_scenario(ax, name, p, r, pm, sc):
    pf = sp.lambdify((y1, y2), sp.sympify(0)); safe = r["P"] > 0; pp = sp.lambdify((y1, y2), p)
    ax.add_patch(plt.Rectangle((-6.5, Y_EDGE), 13, 7.0, color="#f5f0e6", zorder=0))
    ax.contourf(Y1, Y2, r["P"] / r["pmax"], levels=[0, 0.5, 10], colors=["#dfe9f5", "#b9d3ee"]); ax.contour(Y1, Y2, r["P"], levels=[0], colors="#2166ac", linewidths=1.8)
    ax.contour(Y1, Y2, np.where(safe, pp(Y1, Y2), np.nan), levels=[0], colors="0.3", linestyles="--", linewidths=1.1)
    A, B = sc["A"], sc["B"]; ax.plot([A[0], B[0]], [A[1], B[1]], color="#d6604d", ls=":", lw=1.2)
    ax.add_patch(plt.Circle(B, RT, color="#fee0d2", ec="#d6604d", lw=1.8, zorder=5)); ax.add_patch(plt.Rectangle(B - 0.3, 0.6, 0.6, color="none", ec="#d6604d", lw=1.2, ls="--", zorder=6))
    ax.add_patch(plt.Rectangle(A - 0.3, 0.6, 0.6, color="#ff8c00", ec="k", lw=1.0, zorder=8))
    for c, rr in sc["obs"]: ax.add_patch(plt.Circle(c, rr, color="0.55", ec="k", zorder=6))
    th = np.linspace(0, np.pi, 200); ax.plot(8 * np.cos(th), 8 * np.sin(th), "k-", lw=1)
    for rr in (2.5, 7): ax.plot(rr * np.cos(th), rr * np.sin(th), "k:", lw=0.7)
    q = ik(A.reshape(2, 1))[:, 0]; e = np.array([4 * np.cos(q[0]), 4 * np.sin(q[0])]); d = A / np.linalg.norm(A); e_up = 2 * d * (d @ e) - e
    oc = sc["obs"][0][0]; e = e_up if np.hypot(*(e_up - oc)) > np.hypot(*(e - oc)) else e   # elbow branch farther from the obstacle
    ax.plot([0, e[0]], [0, e[1]], "o-", lw=6, color="0.5", alpha=0.6, zorder=7); ax.plot([e[0], A[0]], [e[1], A[1]], "o-", lw=6, color="0.5", alpha=0.6, zorder=7)
    ax.text(A[0] + 0.5, A[1] - 0.9, "A (pick)", fontsize=9); ax.text(B[0] - 1.9, B[1] - 1.5, "B (place)", fontsize=9)
    ax.plot(0, 0, "ks", ms=7); ax.set_aspect("equal"); ax.set_xlim(-7.2, 7.2); ax.set_ylim(-1.2, 8.2); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{name}\n|y| in [{r['rmin']:.1f}, {r['rmax']:.1f}], deg {r['deg']}, length {pm['length']:.0f}, width {pm['wmed']:.1f}, psi(A) {100*sc['psiA']:.0f} % of max, corridor beyond A {sc['beyond']:.0f} deg, obstacle clearance {min(sc['clear']):.2f}", fontsize=7.5)


if __name__ == "__main__":
    sel = sys.argv[1:]; good = []
    for name, kw in CANDS:
        if sel and not any(s in name for s in sel): continue
        cfg = dict(base); cfg.update(kw); p, psi = build_front(**cfg); r = analyse(psi)
        if r is None: print(name, "no corridor"); continue
        pm = path_metrics(p, psi, r, cfg["rho0"]); sc = scenario_checks(p, psi, r, pm); safe = r["P"] > 0; y2min = Y2[safe].min()
        ok = r["n_other"] == 0 and r["rmin"] >= 2.6 and r["rmax"] <= 7.0 and r["cert"] >= 0.95 and y2min >= Y_EDGE and sc["psiA"] >= 0.25 and sc["beyond"] >= 6 and min(sc["clear"]) >= MARGIN and sc["blocked"]
        print(f"{name}: |y| [{r['rmin']:.2f},{r['rmax']:.2f}] y2>={y2min:.2f} deg {r['deg']} others {r['n_other']} cert {100*r['cert']:.0f}% half {100*r['half']:.0f}% span {pm['span']:.0f} width {pm['wmin']:.2f}/{pm['wmed']:.2f}/{pm['wmax']:.2f} length {pm['length']:.1f} B ({sc['B'][0]:.2f},{sc['B'][1]:.2f}) psi(A) {100*sc['psiA']:.0f}% beyond A {sc['beyond']:.0f} deg clearance {min(sc['clear']):.2f} {'<-- ok' if ok else ''}", flush=True)
        if ok: good.append((name, p, r, pm, sc))
    if good:
        n = len(good); nc = 4; nr = int(np.ceil(n / nc)); fig, axs = plt.subplots(nr, nc, figsize=(4.6 * nc, 3.5 * nr)); axs = np.atleast_1d(axs).ravel()
        for ax, (name, p, r, pm, sc) in zip(axs, good): draw_scenario(ax, name, p, r, pm, sc)
        for ax in axs[n:]: ax.axis("off")
        fig.tight_layout(); fig.savefig("figs/arm_scenario_cands.png", dpi=110); print(f"saved figs/arm_scenario_cands.png with {n} admissible candidates")
