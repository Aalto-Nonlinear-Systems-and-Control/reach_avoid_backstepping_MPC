"""Corridors on a bench in front of the base (base at the bench edge, bench = {y2 >= Y_EDGE}), for a workcell layout.
E  front arc:  p = |y|^2 - rho0^2 - sum_n A_n Re(e^{-i phi_n} z^n) / rho0^(n-2); closing at the back -eps (y2 - rho0)^(2m); weight to the left.
F  lateral:    p = y2 - c(y1), c = y0 + A s(y1 / L), s(x) = (3 sqrt3 / 2) x (1 - x^2) (S) or x (1 - x^2)(1 - x^2 / b^2) scaled (double S);
               tube h^2 - p^2, closing -eps (y1 - yc)^(2m), weight w = 1 + kappa (6 - y1)^2 (increasing to the left, the placement side)."""
import sys, itertools, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from screen_zigzag import analyse, y1, y2, r2, Y1, Y2, RT
from screen_paths import harmonic, path_metrics
Y_EDGE = 0.8


def build_front(rho0, harms, h, kappa, eps, m, norm=True):
    p = r2 - rho0**2 - sum(A * harmonic(n, np.radians(phi)) / rho0**(n - 2) for n, A, phi in harms)
    u = r2 / rho0**2; N = (3 - 3 * u + u**2) if norm else 1; w = 1 + kappa * (6 - y1)**2
    return p, sp.expand(w * (h**2 - p**2 * N / (2 * rho0)**2) - eps * (y2 - rho0)**(2 * m))


def build_lateral(y0, A, L, h, kappa, eps, m, yc, shape="S", b=0.6):
    x = y1 / L
    if shape == "S": s = sp.Rational(3, 2) * sp.sqrt(3) * x * (1 - x**2)
    else:
        xx = sp.symbols("xx"); raw = xx * (1 - xx**2) * (1 - xx**2 / b**2); f = sp.lambdify(xx, raw); s = (raw / float(max(abs(f(np.linspace(-1, 1, 2001)))))).subs(xx, x)
    c = y0 + A * s; p = y2 - c; w = 1 + kappa * (6 - y1)**2
    return p, sp.expand(w * (h**2 - p**2) - eps * (y1 - yc)**(2 * m))


def metrics_lateral(p, psi, r):
    """Path metrics for a graph path y2 = c(y1): sample y1 and measure the width along the normal."""
    pf = sp.lambdify((y1, y2), psi); c = sp.solve(p, y2)[0]; cf = sp.lambdify(y1, c); cp = sp.lambdify(y1, sp.diff(c, y1))
    xs = np.linspace(-7, 7, 561); X = np.array([[x, cf(x)] for x in xs]); inside = np.array([pf(*x) > 0 for x in X]); X = X[inside]
    if len(X) < 10: return None
    vals = np.array([pf(*x) for x in X]); core = vals >= 0.3 * vals.max(); widths = []
    from scipy.optimize import brentq
    for x in X[core]:
        n = np.array([-cp(x[0]), 1.0]); n /= np.linalg.norm(n); ws = []
        for sgn in (1, -1):
            f = lambda t: pf(*(x + sgn * t * n)); ts = np.linspace(0, 3, 61); v = np.array([f(t) for t in ts]); j = np.where(v <= 0)[0]
            ws.append(brentq(f, ts[j[0] - 1], ts[j[0]]) if len(j) and j[0] > 0 else np.nan)
        widths.append(sum(ws))
    widths = np.array(widths); seg = np.r_[0, np.cumsum(np.hypot(np.diff(X[:, 0]), np.diff(X[:, 1])))]; k = np.argmin(np.hypot(X[:, 0] - r["tc"][0], X[:, 1] - r["tc"][1]))
    return dict(span=float(X[:, 0].max() - X[:, 0].min()), wmin=np.nanmin(widths), wmed=np.nanmedian(widths), wmax=np.nanmax(widths), bends=0, length=seg[-1], to_ends=(seg[k], seg[-1] - seg[k]), X=X)


def report(name, p, psi, r, pm, extra=""):
    safe = r["P"] > 0; y2min = Y2[safe].min(); ok = r["n_other"] == 0 and r["rmin"] >= 2.6 and r["rmax"] <= 7.0 and r["cert"] >= 0.95 and y2min >= Y_EDGE
    print(f"{name}: |y| [{r['rmin']:.2f},{r['rmax']:.2f}] y2>={y2min:.2f} deg {r['deg']} others {r['n_other']} cert {100*r['cert']:.0f}% half {100*r['half']:.0f}% span {pm['span']:.0f} width {pm['wmin']:.2f}/{pm['wmed']:.2f}/{pm['wmax']:.2f} length {pm['length']:.1f} target->ends {pm['to_ends'][0]:.1f}/{pm['to_ends'][1]:.1f} max ({r['tc'][0]:.2f},{r['tc'][1]:.2f}){extra} {'<-- ok' if ok else ''}", flush=True)
    return ok, y2min


def draw_bench(ax, name, p, r, pm, y2min):
    ax.add_patch(plt.Rectangle((-6.5, Y_EDGE), 13, 7.2, color="#f5f0e6", zorder=0)); ax.contourf(Y1, Y2, r["P"] / r["pmax"], levels=[0, 0.5, 10], colors=["#dfe9f5", "#b9d3ee"]); ax.contour(Y1, Y2, r["P"], levels=[0], colors="#2166ac", linewidths=1.5)
    pp = sp.lambdify((y1, y2), p); ax.contour(Y1, Y2, np.where(r["P"] > 0, pp(Y1, Y2), np.nan), levels=[0], colors="0.3", linestyles="--", linewidths=1.0)
    ax.add_patch(plt.Circle(r["tc"], RT, color="#fee0d2", ec="#d6604d", lw=1.5, zorder=5)); th = np.linspace(0, np.pi, 200)
    for rr, ls in ((8, "k-"), (7, "k:"), (2.5, "k:")): ax.plot(rr * np.cos(th), rr * np.sin(th), ls, lw=0.8)
    ax.plot(0, 0, "ks", ms=6); ax.set_aspect("equal"); ax.set_xlim(-7.5, 7.5); ax.set_ylim(-1.5, 8.2); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{name}\n|y| in [{r['rmin']:.1f}, {r['rmax']:.1f}], y2 >= {y2min:.1f}, deg {r['deg']}, span {pm['span']:.0f}, width {pm['wmin']:.1f}..{pm['wmax']:.1f} (med {pm['wmed']:.1f}), length {pm['length']:.0f}, target {min(pm['to_ends']):.1f} before the end", fontsize=7.5)


if __name__ == "__main__":
    sel = sys.argv[1:]; out = []
    E, F = [], []
    for rho0, harms, kappa, (m, eps) in itertools.product((4.5,), ([(3, 0.25, 30)], [(3, 0.25, 60)], [(3, 0.25, 90)], [(3, 0.25, 60), (2, 0.12, 0)], [(3, 0.25, 60), (2, 0.12, 90)], [(4, 0.15, 0)], [(4, 0.15, 45)]), (0.05, 0.1),
                                              ((2, 0.02), (2, 0.04), (3, 0.002), (3, 0.005))):
        E.append((f"E rho0={rho0} harms={harms} kappa={kappa} m={m} eps={eps}", dict(rho0=rho0, harms=harms, h=1.0, kappa=kappa, eps=eps, m=m)))
    for y0, A, L, kappa, (m, eps), yc, shape in itertools.product((4.5, 5.0), (0.8, 1.2), (5.0,), (0.05,), ((2, 0.002), (2, 0.005), (3, 0.0002)), (0.0, 1.0), ("S", "DS")):
        F.append((f"F {shape} y0={y0} A={A} L={L} kappa={kappa} m={m} eps={eps} yc={yc}", dict(y0=y0, A=A, L=L, h=1.0, kappa=kappa, eps=eps, m=m, yc=yc, shape=shape)))
    for name, kw in E + F:
        if sel and not any(s in name for s in sel): continue
        try:
            if name.startswith("E"): p, psi = build_front(**kw); r = analyse(psi); pm = path_metrics(p, psi, r, kw["rho0"]) if r else None
            else: p, psi = build_lateral(**kw); r = analyse(psi); pm = metrics_lateral(p, psi, r) if r else None
        except Exception as ex: print(name, "error", ex); continue
        if r is None or pm is None: print(name, "no corridor"); continue
        ok, y2min = report(name, p, psi, r, pm); out.append((name, p, r, pm, ok, y2min))
    good = [o for o in out if o[4]]
    if good:
        n = len(good); nc = 4; nr = int(np.ceil(n / nc)); fig, axs = plt.subplots(nr, nc, figsize=(4.4 * nc, 3.4 * nr)); axs = np.atleast_1d(axs).ravel()
        for ax, (name, p, r, pm, _, y2min) in zip(axs, good): draw_bench(ax, name, p, r, pm, y2min)
        for ax in axs[n:]: ax.axis("off")
        fig.tight_layout(); fig.savefig("figs/arm_front.png", dpi=110); print(f"saved figs/arm_front.png with {n} admissible candidates")
