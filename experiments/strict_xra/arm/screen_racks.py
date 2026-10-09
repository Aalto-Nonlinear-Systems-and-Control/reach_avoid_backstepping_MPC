"""Cyclic kitting with two vertical racks: objects in the left rack (column y1 = -X0, heights YS), slots in the right rack
(y1 = +X0, same heights). Object k goes to slot SIGMA[k] (crossing paths) and later back.
Path k: y2 = c_k(y1), x = y1 / X1 with X1 = X0 + DX (the path is flattest DX inside the compartment, so it stays level there),
    c_k = mid_k + d_k (3x - x^3)/2 + b (1 - x^2)^2,  mid_k, d_k such that c_k(-X0) = yL and c_k(X0) = yR exactly (bridge of height b)
Certificate (forward; reverse uses w(-y1)):
    psi_k = w(y1) [h^2 (1 + gam c_k'^2) - (y2 - c_k)^2] - eps (y1 / X0)^(2m),   w = 1 + kappa (y1 + a)^2,
    eps = kappa (X0 + a) h^2 X0 / m  (puts dpsi/dy1 = 0 at the slot).
Critical points satisfy y2 = c_k(y1) and a 1D condition in y1, so the maximum is exactly at the destination slot."""
import os, sys, itertools, pickle, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import screen_zigzag
from screen_zigzag import analyse, y1, y2, Y1, Y2
from screen_front import Y_EDGE
from cert_arm import ik

X0 = 3.2; DX = 0.4; YS = (1.8, 2.8, 3.8, 4.8); SIGMA = (2, 0, 3, 1); WALL = 0.5; RACK_OPEN = X0 - 0.25
NAMES = ("red", "blue", "green", "orange"); COLORS = ("#d62728", "#1f77b4", "#2ca02c", "#ff7f0e")
H = 0.38; OBJ = 0.2; RT = 0.35; M = 5; KAPPA = 0.02; A_W = 4.5
FIX_W = 1.2; FIX_Y0 = 1.1; RACK_D = 1.5; VIEW = ((-5.3, 5.3), (-2.9, 6.0)); BRIDGE = 1.0
LAYOUT = os.environ.get("LAYOUT", "wide")
if LAYOUT == "square":      # compact cell for square figures: racks closer and higher; weight and closing scaled with X0 (a ~ X0, kappa ~ 1/X0^2)
    X0 = 2.4; DX = 0.3; YS = (2.3, 3.3, 4.3, 5.3); RACK_OPEN = X0 - 0.25; A_W = 4.5 * X0 / 3.2; KAPPA = 0.02 * (3.2 / X0)**2
    FIX_W = 1.0; FIX_Y0 = 1.2; RACK_D = 1.2; VIEW = ((-3.7, 3.7), (-0.7, 6.7)); BRIDGE = float(os.environ.get("BRIDGE", 0.8))
screen_zigzag.RT = RT


def path(k, b):
    X1 = X0 + DX; x = y1 / X1; x0 = X0 / X1; yL, yR = YS[k], YS[SIGMA[k]]
    S = lambda t: (3 * t - t**3) / 2; B = lambda t: (1 - t**2)**2
    d = (yR - yL) / (2 * S(x0)); mid = (yL + yR) / 2 - b * B(x0)
    return sp.nsimplify(mid, rational=True, tolerance=1e-12) + sp.nsimplify(d, rational=True, tolerance=1e-12) * S(x) + sp.nsimplify(b, rational=True) * B(x)


def wall_margin(t):
    """Smallest distance from the corridor to the compartment walls of its two slots, inside the racks (negative = crosses a wall)."""
    safe = t["r"]["P"] > 0; worst = 9.0
    for side, slot in ((1, t["dst"][1] if t["sgn"] > 0 else t["src"][1]), (-1, t["src"][1] if t["sgn"] > 0 else t["dst"][1])):
        sel = safe & (side * Y1 >= RACK_OPEN)
        if sel.any(): worst = min(worst, float(np.min(WALL - np.abs(Y2[sel] - slot))))
    return worst


def certificate(c, sgn, gam):
    w = 1 + KAPPA * (sgn * y1 + A_W)**2; eps = KAPPA * (X0 + A_W) * H**2 * X0 / M; cp = sp.diff(c, y1)
    return sp.expand(w * (H**2 * (1 + gam * cp**2) - (y2 - c)**2) - eps * (y1 / X0)**(2 * M))


def along_unique(c, sgn, gam):
    """Number of sign changes of d/dy1 of psi along the path inside the corridor (must be 1)."""
    f = sp.lambdify(y1, certificate(c, sgn, gam).subs(y2, c)); xx = np.linspace(-2 * X0, 2 * X0, 8001); fx = f(xx); pos = xx[fx > 0]
    xi = np.linspace(pos[0], pos[-1], 8001); d = np.diff(f(xi)); s = np.sign(d); return int((s[1:] * s[:-1] < 0).sum()), pos[0], pos[-1]


def fixture_clear(r, top):
    safe = r["P"] > 0; px, py = Y1[safe], Y2[safe]; x0, x1 = -FIX_W / 2, FIX_W / 2
    dx = np.maximum(np.maximum(x0 - px, px - x1), 0); dy = np.maximum(np.maximum(FIX_Y0 - py, py - top), 0); return float(np.min(np.hypot(dx, dy)))


def far_positive(psi):
    pf = sp.lambdify((y1, y2), psi); g = np.linspace(-25, 25, 501); A, B = np.meshgrid(g, g); P = pf(A, B); return bool(((P > 0) & (np.hypot(A, B) > 7.6)).any())


def run(b, gam, verbose=True):
    tasks = []
    for k in range(4):
        c = path(k, b)
        for sgn in (1, -1):
            psi = certificate(c, sgn, gam); r = analyse(psi); nchg, xs, xe = along_unique(c, sgn, gam)
            src = np.array([-sgn * X0, YS[k] if sgn > 0 else YS[SIGMA[k]]]); dst = np.array([sgn * X0, YS[SIGMA[k]] if sgn > 0 else YS[k]])
            tasks.append(dict(k=k, sgn=sgn, c=c, psi=psi, r=r, src=src, dst=dst, nchg=nchg, ends=(xs, xe)))
    floors = []
    for t in tasks:
        if t["r"] is None: continue
        sel = (t["r"]["P"] > 0) & (abs(Y1) <= FIX_W / 2 + 0.3); floors.append(Y2[sel].min() if sel.any() else 9.0)
    top = min(floors) - 0.3
    stations = [np.array([-X0, yy]) for yy in YS] + [np.array([X0, yy]) for yy in YS]
    for t in tasks:
        r = t["r"]
        if r is None: t["ok"] = False; continue
        pf = sp.lambdify((y1, y2), t["psi"]); safe = r["P"] > 0; px, py = Y1[safe], Y2[safe]
        others = [s for s in stations if np.linalg.norm(s - t["src"]) > 1e-6 and np.linalg.norm(s - t["dst"]) > 1e-6]
        t.update(cross=min(np.min(np.hypot(px - s[0], py - s[1])) for s in others), fix=fixture_clear(r, top), far=far_positive(t["psi"]), y2min=float(py.min()), wall=wall_margin(t),
                 off=float(np.linalg.norm(np.array(r["tc"]) - t["dst"])), psi_src=pf(*t["src"]) / r["pmax"])
        t["ok"] = (r["n_other"] == 0 and r["cert"] >= 0.95 and t["nchg"] == 1 and r["rmin"] >= 2.6 and r["rmax"] <= 7.0 and not t["far"] and t["y2min"] >= Y_EDGE
                   and t["off"] <= 0.05 and t["psi_src"] >= 0.3 and t["cross"] >= OBJ + 0.05 and t["fix"] >= 0.25 and t["wall"] >= 0.0)
        if verbose:
            print(f"  {NAMES[t['k']]:6s} {'fwd' if t['sgn'] > 0 else 'rev'}: deg {r['deg']} crit {len(r['crit'])} others {r['n_other']} 1D-changes {t['nchg']} cert {100*r['cert']:.0f}% |y| [{r['rmin']:.2f},{r['rmax']:.2f}] "
                  f"y2>={t['y2min']:.2f} far {t['far']} max-slot {t['off']:.3f} psi(src) {100*t['psi_src']:.0f}% cross {t['cross']:.2f} wall {t['wall']:+.3f} fixture {t['fix']:.2f} ends y1 [{t['ends'][0]:.2f},{t['ends'][1]:.2f}] {'ok' if t['ok'] else 'FAIL'}", flush=True)
    return tasks, top


def draw_cell(ax, tasks, top, phase, title, show_all=True, focus=None):
    ax.add_patch(plt.Rectangle((-7.2, Y_EDGE), 14.4, 6.9, color="#f3eee3", zorder=0)); ax.plot([-7.2, 7.2], [Y_EDGE, Y_EDGE], color="#b9ab8f", lw=2, zorder=1)
    for side in (-1, 1):                                                     # racks: back wall, side walls, dividers
        xo = side * RACK_OPEN; xb = side * (X0 + RACK_D); ylo, yhi = YS[0] - WALL, YS[-1] + WALL
        ax.add_patch(plt.Rectangle((min(xo, xb), ylo), abs(xb - xo), yhi - ylo, color="#e4dccb", zorder=1))
        ax.plot([xb, xb], [ylo, yhi], color="#7d705a", lw=2.5, zorder=2)
        for yy in [ylo] + [(YS[i] + YS[i + 1]) / 2 for i in range(3)] + [yhi]: ax.plot([xo, xb], [yy, yy], color="#7d705a", lw=1.6, zorder=2)
    ax.add_patch(plt.Rectangle((-FIX_W / 2, FIX_Y0), FIX_W, top - FIX_Y0, color="0.55", ec="k", zorder=6)); ax.text(0, (FIX_Y0 + top) / 2, "fixture", ha="center", va="center", fontsize=8, color="w", zorder=7)
    sel = [t for t in tasks if (t["sgn"] > 0) == (phase == "fwd")]
    if focus is not None: sel = [t for t in sel if t["k"] == focus]
    for t in sel:
        r = t["r"]; col = COLORS[t["k"]]; cf = sp.lambdify(y1, t["c"]); xs, xe = t["ends"]; xx = np.linspace(xs, xe, 400)
        ax.contourf(Y1, Y2, r["P"], levels=[0, 1e9], colors=[col], alpha=0.14, zorder=2); ax.contour(Y1, Y2, r["P"], levels=[0], colors=[col], linewidths=1.0, zorder=3)
        ax.plot(xx, cf(xx), color=col, lw=1.5, ls="--", zorder=4)
        xm = 0.35 * t["sgn"] * (1 if t["k"] % 2 else -1); ax.annotate("", xy=(xm + 0.25 * t["sgn"], cf(xm + 0.25 * t["sgn"])), xytext=(xm, cf(xm)), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6, mutation_scale=14), zorder=5)
    for k in range(4):
        if phase == "fwd": src, dst = np.array([-X0, YS[k]]), np.array([X0, YS[SIGMA[k]]])
        else: src, dst = np.array([X0, YS[SIGMA[k]]]), np.array([-X0, YS[k]])
        ax.add_patch(plt.Rectangle(src - OBJ, 2 * OBJ, 2 * OBJ, color=COLORS[k], ec="k", lw=0.8, zorder=8))
        ax.add_patch(plt.Rectangle(dst - OBJ, 2 * OBJ, 2 * OBJ, color="none", ec=COLORS[k], lw=1.0, ls=":", zorder=8))
    t0 = sel[0]; A = t0["src"]; q = ik(A.reshape(2, 1))[:, 0]; e = np.array([4 * np.cos(q[0]), 4 * np.sin(q[0])])      # branch q2 > 0, as in the dynamics
    ax.plot([0, e[0]], [0, e[1]], "-", lw=8, color="0.35", alpha=0.6, solid_capstyle="round", zorder=9); ax.plot([e[0], A[0]], [e[1], A[1]], "-", lw=6.5, color="0.35", alpha=0.6, solid_capstyle="round", zorder=9)
    ax.add_patch(plt.Circle(e, 0.22, color="0.25", zorder=10)); ax.add_patch(plt.Circle((0, 0), 0.5, color="0.2", zorder=10))
    ax.set_aspect("equal"); ax.set_xlim(*VIEW[0]); ax.set_ylim(*VIEW[1]); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=10)


def draw_one(ax, t):
    r = t["r"]; col = COLORS[t["k"]]; P = r["P"] / r["pmax"]
    ax.add_patch(plt.Rectangle((-7.2, Y_EDGE), 14.4, 6.9, color="#f3eee3", zorder=0))
    ax.contourf(Y1, Y2, P, levels=[0, 0.5, 2], colors=["#dfe9f5", "#b9d3ee"], zorder=1); ax.contour(Y1, Y2, P, levels=[0], colors="#2166ac", linewidths=1.6, zorder=2)
    ax.contour(Y1, Y2, P, levels=[0.2, 0.4, 0.6, 0.8], colors="0.5", linewidths=0.6, linestyles="-", zorder=2)
    ax.add_patch(plt.Circle(r["tc"], RT, color="#fee0d2", ec="#d6604d", lw=1.6, zorder=5)); ax.plot(*t["src"], "s", color=col, ms=9, mec="k", zorder=6)
    ax.set_aspect("equal"); ax.set_xlim(*VIEW[0]); ax.set_ylim(*VIEW[1]); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{NAMES[t['k']]} object, forward task: X_S (blue), psi >= 0.5 max (dark), level sets, target disc at the slot", fontsize=9)


if __name__ == "__main__":
    configs = [(b, g) for b in (0.9, 1.1) for g in (0.0, 0.5, 1.0)]
    if len(sys.argv) > 2: configs = [(float(sys.argv[1]), float(sys.argv[2]))]
    res = []
    for b, g in configs:
        print(f"bulge={b} gamma={g}", flush=True); tasks, top = run(b, g); nok = sum(t["ok"] for t in tasks); print(f"  -> {nok}/8 pass, fixture top {top:.2f}", flush=True); res.append((b, g, tasks, top, nok))
    b, g, tasks, top, nok = max(res, key=lambda x: (x[4], x[1], x[0]))
    SUF = "_square" if LAYOUT == "square" else ""
    fig, axs = plt.subplots(1, 3, figsize=(20, 5.2) if not SUF else (19, 6.6), gridspec_kw=dict(wspace=0.04))
    draw_cell(axs[0], tasks, top, "fwd", f"forward: objects from the left rack to the right rack (bridge {b}, gamma {g}, {nok}/8 certificates pass)")
    draw_cell(axs[1], tasks, top, "rev", "reverse: back to the left rack along the same tubes")
    draw_one(axs[2], [t for t in tasks if t["k"] == 0 and t["sgn"] > 0][0])
    fig.savefig(f"figs/arm_racks{SUF}.png", dpi=130, bbox_inches="tight"); print(f"saved figs/arm_racks{SUF}.png", b, g)
    with open(f"data/racks_best{SUF}.pkl", "wb") as f: pickle.dump(dict(b=b, gam=g, top=top, X0=X0, YS=YS, SIGMA=SIGMA, H=H, M=M, KAPPA=KAPPA, A_W=A_W, psi={(t["k"], t["sgn"]): str(t["psi"]) for t in tasks}), f)
