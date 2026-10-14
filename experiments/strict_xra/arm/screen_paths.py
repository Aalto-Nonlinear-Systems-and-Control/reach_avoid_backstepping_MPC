"""Task corridors for the top-view arm: tube of uniform half-width h around a polynomial nominal path p(y) = 0.
p(y) = |y|^2 - rho0^2 - sum_n A_n Re(e^{-i phi_n} z^n) / rho0^(n-2),  z = y1 + i y2   (radius rho0 modulated by harmonics n)
psi  = w(y) [h^2 - p^2 N(|y|^2 / rho0^2) / (2 rho0)^2] - eps (y . d - rho0)^(2 m),  w = 1 + kappa (y . e + 6)^2
N(u) = 1 (plain tube, width ~ 2 h rho0 / |y|), or N(u) = 3 - 3 u + u^2 (second-order Taylor of 1 / u at u = 1, positive for all u):
the width then stays 2 h within about 1 % for |y| / rho0 in [0.85, 1.15] while psi on the path stays exactly w h^2 (unimodal).
Metrics per candidate: admissibility (unique maximum, annulus 2.6..7, certifiable), angular span, width min/median/max
on the core of the path, number of inflections (bends), distance from the target centre to the far end of the path."""
import sys, itertools, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
from screen_zigzag import analyse, y1, y2, r2, Y1, Y2, RT
_a, _b = sp.symbols("a b", real=True)


def harmonic(n, phi):
    e = sp.expand((_a + sp.I * _b)**n * (sp.cos(phi) - sp.I * sp.sin(phi)))
    return sp.expand(sp.re(e)).subs({_a: y1, _b: y2}, simultaneous=True)


def build(rho0, harms, h, kappa, e, eps, d, m, norm="taylor2", weight="linear", c=5.0):
    """weight 'linear': w = 1 + kappa (y.e + 6)^2, monotone along at most a half circle.
    weight 'long': w = 1 + kappa (W + c)^2 with W = (y.e_perp)(rho0 - y.e) / rho0^2 (= sin(phi)(1 - cos(phi)) on the circle,
    phi the angle from e), increasing along the arc phi in [-120, 120] deg: a 240-degree corridor with a single maximum."""
    p = r2 - rho0**2 - sum(A * harmonic(n, np.radians(phi)) / rho0**(n - 2) for n, A, phi in harms)
    u = r2 / rho0**2; N = {"none": 1, "taylor2": 3 - 3 * u + u**2, "grad": None}[norm]
    pen = (sp.diff(p, y1)**2 + sp.diff(p, y2)**2) if norm == "grad" else None
    ep = (-e[1], e[0]); sn = (ep[0] * y1 + ep[1] * y2) / rho0; cs = (e[0] * y1 + e[1] * y2) / rho0   # sin, cos of the angle from e on |y| = rho0
    if weight == "linear": w = 1 + kappa * (e[0] * y1 + e[1] * y2 + 6)**2
    elif weight == "long": w = 1 + kappa * (sn * (1 - cs) + c)**2          # derivative (2 cos + 1)(1 - cos): double zero at phi = 0
    else: w = 1 + kappa * (2 * sn - sp.Rational(2, 3) * sn**3 - sn * cs / 2 + c)**2   # "long3": derivative (cos + 1/2)(2 cos^2 - 2 cos + 1) > 0 on (-120, 120) deg
    core = (h**2 * pen - p**2) if norm == "grad" else (h**2 * (2 * rho0)**2 - p**2 * N)
    psi = sp.expand(w * core / (2 * rho0)**2 - eps * (d[0] * y1 + d[1] * y2 - rho0)**(2 * m))
    return p, psi


def path_metrics(p, psi, r, rho0):
    pf = sp.lambdify((y1, y2), psi); pp = sp.lambdify((y1, y2), p); G = sp.lambdify((y1, y2), [sp.diff(p, y1), sp.diff(p, y2)])
    th = np.linspace(-np.pi, np.pi, 721); pts = []
    for t in th:
        f = lambda rr: pp(rr * np.cos(t), rr * np.sin(t)); rs = np.linspace(2.0, 7.6, 57); v = f(rs); idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
        if len(idx) == 0: continue
        i = idx[np.argmin(abs(rs[idx] - rho0))]; rr = brentq(f, rs[i], rs[i + 1]); x = np.array([rr * np.cos(t), rr * np.sin(t)])
        if pf(*x) > 0: pts.append((t, x))
    if len(pts) < 10: return None
    thetas = np.array([t for t, _ in pts]); X = np.array([x for _, x in pts])
    # angular span; order the points along the arc (the inactive angles form one gap, possibly across +-180 deg)
    occ = np.zeros(len(th), bool); occ[np.searchsorted(th, thetas)] = True; span = occ.sum() * (th[1] - th[0])
    dth = np.diff(np.r_[thetas, thetas[0] + 2 * np.pi]); k = int(np.argmax(dth)); order = np.r_[np.arange(k + 1, len(pts)), np.arange(0, k + 1)]; thetas, X = thetas[order], X[order]
    # width along the normal at core points (psi on the path >= 0.3 of its maximum on the path)
    vals = np.array([pf(*x) for x in X]); core = vals >= 0.3 * vals.max(); widths = []
    for x in X[core]:
        n = np.array(G(*x), float); n /= np.linalg.norm(n); ws = []
        for sgn in (1, -1):
            f = lambda t: pf(*(x + sgn * t * n)); ts = np.linspace(0, 3, 61); v = np.array([f(t) for t in ts]); j = np.where(v <= 0)[0]
            ws.append(brentq(f, ts[j[0] - 1], ts[j[0]]) if len(j) and j[0] > 0 else np.nan)
        widths.append(sum(ws))
    widths = np.array(widths)
    # inflections of the path within the active arc: sign changes of the signed curvature
    xs, ys = X[:, 0], X[:, 1]; dx, dy = np.gradient(xs), np.gradient(ys); ddx, ddy = np.gradient(dx), np.gradient(dy); curv = (dx * ddy - dy * ddx) / (dx**2 + dy**2)**1.5
    sc = np.sign(curv[core]); bends = int((sc[1:] * sc[:-1] < 0).sum())
    # arc length from the target centre to the two ends of the active path
    seg = np.r_[0, np.cumsum(np.hypot(np.diff(xs), np.diff(ys)))]; k = np.argmin(np.hypot(xs - r["tc"][0], ys - r["tc"][1]))
    return dict(span=np.degrees(span), wmin=np.nanmin(widths), wmed=np.nanmedian(widths), wmax=np.nanmax(widths), bends=bends, length=seg[-1], to_ends=(seg[k], seg[-1] - seg[k]), X=X)


def draw(ax, name, p, r, pm):
    ax.contourf(Y1, Y2, r["P"] / r["pmax"], levels=[0, 0.5, 10], colors=["#dfe9f5", "#b9d3ee"]); ax.contour(Y1, Y2, r["P"], levels=[0], colors="#2166ac", linewidths=1.5)
    pp = sp.lambdify((y1, y2), p); ax.contour(Y1, Y2, np.where(r["P"] > 0, pp(Y1, Y2), np.nan), levels=[0], colors="0.3", linestyles="--", linewidths=1.0)
    ax.add_patch(plt.Circle(r["tc"], RT, color="#fee0d2", ec="#d6604d", lw=1.5, zorder=5))
    th = np.linspace(0, 2 * np.pi, 400)
    for rr, ls in ((8, "k-"), (7, "k:"), (2.5, "k:")): ax.plot(rr * np.cos(th), rr * np.sin(th), ls, lw=0.8)
    ax.plot(0, 0, "ks", ms=5); ax.set_aspect("equal"); ax.set_xlim(-7.5, 8.5); ax.set_ylim(-7, 7); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{name}\n|y| in [{r['rmin']:.1f}, {r['rmax']:.1f}], deg {r['deg']}, span {pm['span']:.0f} deg, width {pm['wmin']:.2f}..{pm['wmax']:.2f}, bends {pm['bends']}, end dist {min(pm['to_ends']):.1f}", fontsize=7.5)


CANDIDATES = []
def add(name, **kw): CANDIDATES.append((name, kw))
base = dict(rho0=4.25, h=1.0, kappa=0.3, e=(0, 1), eps=0.02, d=(1, 0), m=2, norm="taylor2")
add("W3 plain tube (current)", harms=[(3, 0.25, 0)], norm="none")
add("W3 normalised tube", harms=[(3, 0.25, 0)])
for phi in (20, 40): add(f"W3 phase {phi}", harms=[(3, 0.25, phi)])
for A2, phi in itertools.product((0.12, 0.2), (0, 45, 90, 135)): add(f"W3 + W2 A2={A2} phi2={phi}", harms=[(3, 0.25, 0), (2, A2, phi)])
for A4, phi in itertools.product((0.06, 0.1), (0, 22, 45, 67)): add(f"W3 + W4 A4={A4} phi4={phi}", harms=[(3, 0.25, 0), (4, A4, phi)])
for A5, phi in itertools.product((0.04, 0.06), (0, 36)): add(f"W3 + W5 A5={A5} phi5={phi}", harms=[(3, 0.2, 0), (5, A5, phi)])
for A1, phi in itertools.product((0.3, 0.5), (0, 90, 180, 270)): add(f"W3 shifted A1={A1} phi1={phi}", harms=[(3, 0.25, 0), (1, A1, phi)])
for A4, phi in itertools.product((0.12, 0.18), (0, 22, 45)): add(f"W4 A4={A4} phi4={phi}", harms=[(4, A4, phi)])
for A2, A4 in itertools.product((0.1, 0.15), (0.06, 0.1)): add(f"W2 + W3 + W4 A2={A2} A4={A4}", harms=[(2, A2, 45), (3, 0.22, 0), (4, A4, 22)])
for kappa, eps, phi3, norm in itertools.product((0.05, 0.1, 0.2), (0.002, 0.005, 0.01), (0, 30, 60), ("none", "taylor2")):
    add(f"L W3 kappa={kappa} eps={eps} phi3={phi3} {norm}", harms=[(3, 0.25, phi3)], weight="long", e=(0, 1), d=(0, 1), kappa=kappa, eps=eps, norm=norm)
for kappa, (m, eps), phi3 in itertools.product((0.05, 0.1, 0.2), ((2, 0.002), (2, 0.005), (3, 0.0002), (3, 0.0005)), (0, 30, 60)):
    add(f"L3 W3 kappa={kappa} m={m} eps={eps} phi3={phi3}", harms=[(3, 0.25, phi3)], weight="long3", c=4.0, e=(0, 1), d=(0, 1), kappa=kappa, eps=eps, m=m, norm="taylor2")
for kappa, (m, eps), (A2, phi2) in itertools.product((0.1, 0.2), ((2, 0.002), (3, 0.0002)), ((0.12, 45), (0.12, 135), (0.12, 225), (0.12, 315))):
    add(f"L3 W3+W2 kappa={kappa} m={m} eps={eps} A2={A2} phi2={phi2}", harms=[(3, 0.25, 60), (2, A2, phi2)], weight="long3", c=4.0, e=(0, 1), d=(0, 1), kappa=kappa, eps=eps, m=m, norm="taylor2")

if __name__ == "__main__":
    sel = sys.argv[1:]; out = []
    for name, kw in CANDIDATES:
        if sel and not any(s in name for s in sel): continue
        cfg = dict(base); cfg.update(kw); p, psi = build(**cfg); r = analyse(psi)
        if r is None: print(f"{name}: no safe set", flush=True); continue
        pm = path_metrics(p, psi, r, cfg["rho0"])
        if pm is None: print(f"{name}: path not inside the safe set", flush=True); continue
        ok = r["n_other"] == 0 and r["rmin"] >= 2.6 and r["rmax"] <= 7.0 and r["cert"] >= 0.95
        print(f"{name}: |y| [{r['rmin']:.2f},{r['rmax']:.2f}] deg {r['deg']} others {r['n_other']} cert {100*r['cert']:.0f}% half {100*r['half']:.0f}% span {pm['span']:.0f} width {pm['wmin']:.2f}/{pm['wmed']:.2f}/{pm['wmax']:.2f} bends {pm['bends']} length {pm['length']:.1f} target->ends {pm['to_ends'][0]:.1f}/{pm['to_ends'][1]:.1f} max ({r['tc'][0]:.2f},{r['tc'][1]:.2f}) {'<-- ok' if ok else ''}", flush=True)
        out.append((name, p, r, pm, ok))
    good = [o for o in out if o[4]]
    if good:
        n = len(good); nc = 4; nr = int(np.ceil(n / nc)); fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.9 * nr)); axs = np.atleast_1d(axs).ravel()
        for ax, (name, p, r, pm, _) in zip(axs, good): draw(ax, name, p, r, pm)
        for ax in axs[n:]: ax.axis("off")
        fig.tight_layout(); fig.savefig("figs/arm_paths.png", dpi=120); print(f"saved figs/arm_paths.png with {n} admissible candidates")
