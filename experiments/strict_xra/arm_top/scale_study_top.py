"""Base design checks and the choice of mu0 for the top-view arm (s = 1): ECGBF margin on the grid, |k1|, share of rest
states of the corridor inside X_RA, largest admissible period T_eps(1) (sampled one-step increase on the terminal set),
torque of the terminal controller on the terminal set and on X_RA. With g = 0: T_eps(s) = T_eps(1)/s, tau_s = s^2 tau_1."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_top import CertArm, fk, ik, phi_np, BOX, TASK
K1 = os.path.join(HERE, "k1", sys.argv[1] if len(sys.argv) > 1 else f"{TASK}_d8_l{os.environ.get('LAM', '0.05')}.txt"); DFRAC = 0.5
c0 = CertArm(K1, mu=1.0); g = c0.margin_on_grid(n=601)
print(f"{os.path.basename(K1)}: lambda {c0.lam:.4f}, delta_s {c0.ds:.2e}, shift {c0.shift:.2e} | grid: min ECGBF margin {g['min_margin']:+.3e}, max |k1_i| {g['kmax']:.3f}, |y| in [{g['rmin']:.2f}, {g['rmax']:.2f}]")
r = np.random.default_rng(5); Y = np.array([r.uniform(*BOX[0], 80000), r.uniform(*BOX[1], 80000)]); keep = (np.broadcast_to(c0.psi(*Y), (Y.shape[1],)) > 0) & (phi_np(*Y) > 0); Y = Y[:, keep][:, :3000]
X_rest = np.vstack([ik(Y), np.zeros((2, Y.shape[1]))])


def T_eps(c, x, lo=1e-3, hi=3.0):
    tau = c.uc(x); P0 = c.Psi(x)
    def ok(T):
        xn = c.zoh(x, tau, T, n=max(10, int(T / 0.002))); d = c.Psi(xn) - P0; k = phi_np(*fk(xn[:2])) > 0; return (d[k].min() > 0) if k.any() else True
    if not ok(lo): return 0.0
    for _ in range(14):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


rows = []
for mu in (0.1, 0.3, 1.0, 3.0):
    c = CertArm(K1, mu=mu); delta = DFRAC * (1 - c.shift); share = float((c.Psi(X_rest) >= 0).mean())
    xd, _ = c.sample(delta, n_y=60000, seed=11); xd = xd[:, phi_np(*fk(xd[:2])) > 0]; x0, _ = c.sample(0.0, n_y=60000, seed=12); x0 = x0[:, phi_np(*fk(x0[:2])) > 0]
    te = T_eps(c, xd[:, :4000]); td = np.abs(c.uc(xd)).max(1); t0 = np.abs(c.uc(x0)).max(1); tr = np.abs(c.uc(X_rest)).max(1)
    rows.append(dict(mu=mu, share_rest=share, T_eps1=te, tau_terminal=td.tolist(), tau_xra=t0.tolist(), tau_rest=tr.tolist(), n_terminal=xd.shape[1]))
    print(f"mu0 {mu:<4g}: rest states in X_RA {100*share:5.1f}% | T_eps(1) {1e3*te:8.1f} ms | max|tau_1| terminal set {np.round(td, 2)}, X_RA {np.round(t0, 2)}, rest states {np.round(tr, 2)}", flush=True)
json.dump(dict(k1=os.path.basename(K1), lam=c0.lam, ds=c0.ds, grid=g, rows=rows), open(os.path.join(HERE, "data", f"scale_study_{TASK}.json"), "w"), indent=1)
