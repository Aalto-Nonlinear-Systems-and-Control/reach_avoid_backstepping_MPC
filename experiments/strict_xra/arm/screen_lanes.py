"""Kitting scenario: four objects on the right (A_i), four slots on the left (B_i, mirror images), one reference path per object
bridging over a central box obstacle; four nested lanes of equal absolute bulge a and half-width h.
Lane i: p_i = |y|^2 - rho_i^2 - (2 a / rho_i) Re(e^{-i 3 pi/2} z^3) / rho_i,  psi_i = w [h^2 - p_i^2 N_i / (2 rho_i)^2] - eps (y2 - rho_i)^(2 m).
Checks per lane: the usual ones (unique maximum, |y| in [2.6, 7], certifiable, on the bench), A_i well inside, corridor beyond A_i,
box obstacle outside the corridor, objects and slots of the other lanes outside the corridor (with the object half-size)."""
import sys, itertools, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import screen_zigzag
from screen_zigzag import analyse, y1, y2, Y1, Y2
from screen_front import build_front, Y_EDGE
from screen_paths import path_metrics
from cert_arm import ik
RHOS = (3.55, 4.28, 5.02, 5.75); A_BULGE = 0.5; H = 0.45; RT = 0.45; OBJ = 0.2            # object half-size
BOX = ((-0.6, 1.5), 1.2, 1.7)                                                         # (lower-left corner), width, height
COLORS = ("#1b9e77", "#d95f02", "#7570b3", "#e7298a")


def lane(i, kappa, eps, m, h=H):
    rho = RHOS[i]; return build_front(rho0=rho, harms=[(3, 2 * A_BULGE / rho, 270)], h=h, kappa=kappa, eps=eps, m=m)


def box_clearance(r):
    (x0, y0), bw, bh = BOX; safe = r["P"] > 0; px, py = Y1[safe], Y2[safe]
    dx = np.maximum(np.maximum(x0 - px, px - (x0 + bw)), 0); dy = np.maximum(np.maximum(y0 - py, py - (y0 + bh)), 0)
    return float(np.min(np.hypot(dx, dy)))


def solve_lanes(kappa, eps, m, verbose=True):
    screen_zigzag.RT = RT; out = []
    for i in range(4):
        p, psi = lane(i, kappa, eps, m); r = analyse(psi)
        if r is None: out.append(None); continue
        pm = path_metrics(p, psi, r, RHOS[i]); pf = sp.lambdify((y1, y2), psi); B = np.array(r["tc"]); A = np.array([-B[0], B[1]])
        X = pm["X"]; ths = np.degrees(np.arctan2(X[:, 1], X[:, 0])); beyond = np.degrees(np.arctan2(A[1], A[0])) - ths.min()
        out.append(dict(p=p, psi=psi, r=r, pm=pm, pf=pf, A=A, B=B, beyond=beyond, psiA=pf(*A) / r["pmax"], box=box_clearance(r)))
    # cross checks: the objects and slots of the other lanes must lie outside this lane's corridor by the object half-size
    for i, L in enumerate(out):
        if L is None: continue
        safe = L["r"]["P"] > 0; px, py = Y1[safe], Y2[safe]; dmin = np.inf
        for j, M in enumerate(out):
            if j == i or M is None: continue
            for pt in (M["A"], M["B"]): dmin = min(dmin, np.min(np.hypot(px - pt[0], py - pt[1])))
        L["cross"] = dmin; r = L["r"]; safe = r["P"] > 0; y2min = Y2[safe].min()
        L["ok"] = r["n_other"] == 0 and r["rmin"] >= 2.6 and r["rmax"] <= 7.05 and r["cert"] >= 0.95 and y2min >= Y_EDGE and L["psiA"] >= 0.25 and L["beyond"] >= 5 and L["box"] >= 0.2 and dmin >= OBJ + 0.05
        if verbose: print(f"  lane {i+1} rho={RHOS[i]}: |y| [{r['rmin']:.2f},{r['rmax']:.2f}] y2>={y2min:.2f} deg {r['deg']} others {r['n_other']} cert {100*r['cert']:.0f}% half {100*r['half']:.0f}% width {L['pm']['wmin']:.2f}/{L['pm']['wmed']:.2f}/{L['pm']['wmax']:.2f} length {L['pm']['length']:.1f} B ({L['B'][0]:.2f},{L['B'][1]:.2f}) psi(A) {100*L['psiA']:.0f}% beyond {L['beyond']:.0f} deg box {L['box']:.2f} cross {dmin:.2f} {'ok' if L['ok'] else 'FAIL'}", flush=True)
    return out


def draw_lanes(ax, lanes, title):
    ax.add_patch(plt.Rectangle((-6.5, Y_EDGE), 13, 7.0, color="#f5f0e6", zorder=0))
    (x0, y0), bw, bh = BOX; ax.add_patch(plt.Rectangle((x0, y0), bw, bh, color="0.55", ec="k", zorder=6))
    for i, L in enumerate(lanes):
        if L is None: continue
        r = L["r"]; col = COLORS[i]; pp = sp.lambdify((y1, y2), L["p"])
        ax.contourf(Y1, Y2, r["P"], levels=[0, 1e9], colors=[col], alpha=0.22); ax.contour(Y1, Y2, r["P"], levels=[0], colors=[col], linewidths=1.4)
        ax.contour(Y1, Y2, np.where(r["P"] > 0, pp(Y1, Y2), np.nan), levels=[0], colors=[col], linestyles="--", linewidths=0.9)
        A, B = L["A"], L["B"]; ax.add_patch(plt.Circle(B, RT, color="none", ec=col, lw=1.4, zorder=5)); ax.add_patch(plt.Rectangle(B - OBJ, 2 * OBJ, 2 * OBJ, color="none", ec=col, lw=1.2, ls="--", zorder=6))
        ax.add_patch(plt.Rectangle(A - OBJ, 2 * OBJ, 2 * OBJ, color=col, ec="k", lw=0.8, zorder=8))
    th = np.linspace(0, np.pi, 200); ax.plot(8 * np.cos(th), 8 * np.sin(th), "k-", lw=1)
    for rr in (2.5, 7): ax.plot(rr * np.cos(th), rr * np.sin(th), "k:", lw=0.7)
    L0 = lanes[1] or lanes[0]
    if L0 is not None:
        A = L0["A"]; q = ik(A.reshape(2, 1))[:, 0]; e = np.array([4 * np.cos(q[0]), 4 * np.sin(q[0])])
        ax.plot([0, e[0]], [0, e[1]], "o-", lw=5, color="0.5", alpha=0.5, zorder=7); ax.plot([e[0], A[0]], [e[1], A[1]], "o-", lw=5, color="0.5", alpha=0.5, zorder=7)
    ax.text(5.2, 1.2, "objects", fontsize=9); ax.text(-6.4, 1.2, "slots", fontsize=9)
    ax.plot(0, 0, "ks", ms=7); ax.set_aspect("equal"); ax.set_xlim(-7.2, 7.2); ax.set_ylim(-1.2, 8.2); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)


if __name__ == "__main__":
    results = []
    for kappa, (m, eps) in itertools.product((0.03,), ((3, 0.003),)):
        print(f"kappa={kappa} m={m} eps={eps}", flush=True); lanes = solve_lanes(kappa, eps, m); nok = sum(1 for L in lanes if L and L["ok"]); results.append((kappa, m, eps, lanes, nok))
    good = [x for x in results if x[4] == 4] or sorted(results, key=lambda x: -x[4])[:4]
    n = len(good); nc = min(2, n); nr = int(np.ceil(n / nc)); fig, axs = plt.subplots(nr, nc, figsize=(7.5 * nc, 5.2 * nr)); axs = np.atleast_1d(axs).ravel()
    for ax, (kappa, m, eps, lanes, nok) in zip(axs, good): draw_lanes(ax, lanes, f"Kitting scenario: four objects (right) to four slots (left), one tolerance corridor per task\nhalf-width h={H}, bulge a={A_BULGE}, rho={RHOS}, kappa={kappa}, m={m}, eps={eps}: {nok}/4 lanes pass all checks")
    for ax in axs[n:]: ax.axis("off")
    fig.tight_layout(); fig.savefig("figs/arm_lanes.png", dpi=120); print("saved figs/arm_lanes.png")
