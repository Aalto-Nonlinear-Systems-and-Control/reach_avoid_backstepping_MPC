"""Scale and torque bounds of the top-view arm (g = 0, exact time scaling): s = 0.8 T_eps(1) / T with T_eps(1) the largest
sampled period with a one-step increase of Psi on the terminal set at s = 1; bounds U_MAX = 1.25 s^2 max|tau_1| on the
terminal set (rounded up), so s_U = 1.12 s; URA at URA_FAC * s. Checks T_eps(s) = T_eps(1)/s and reports the
second-order bound. Writes data/scale_<TASK>.json."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_top import CertArm, fk, phi_np, TASK
K1 = os.path.join(HERE, "k1", f"{TASK}_d8_l{os.environ.get('LAM', '0.05')}.txt"); MU0, DFRAC, T, URA_FAC = 0.1, 0.5, 0.01, 2.0
c1 = CertArm(K1, mu=MU0); delta = DFRAC * (1 - c1.shift); xd, _ = c1.sample(delta, n_y=200000, seed=11); xd = xd[:, phi_np(*fk(xd[:2])) > 0]
x0, _ = c1.sample(0.0, n_y=200000, seed=12); x0 = x0[:, phi_np(*fk(x0[:2])) > 0]


def T_eps(c, x, lo=1e-4, hi=3.0):
    tau = c.uc(x); P0 = c.Psi(x)
    def ok(T):
        xn = c.zoh(x, tau, T, n=max(10, int(T / 0.0005))); d = c.Psi(xn) - P0; k = phi_np(*fk(xn[:2])) > 0; return (d[k].min() > 0) if k.any() else True
    if not ok(lo): return 0.0
    for _ in range(16):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


Te1 = T_eps(c1, xd[:, :8000]); s = 0.8 * Te1 / T; t1 = np.abs(c1.uc(xd)).max(1); U_MAX = np.ceil(1.25 * s * s * t1 * 2) / 2
cs = c1.scaled(s); xs = xd.copy(); xs[2:] *= s; Tes = T_eps(cs, xs[:, :8000], lo=1e-5, hi=1.0)
tau_s = cs.uc(xs); s_U = s * float(np.sqrt((U_MAX[:, None] / np.maximum(np.abs(tau_s), 1e-12)).min()))
pd = cs.Psidot(xs, tau_s); eta = float(pd.min()); _, M2, _ = cs.zoh(xs[:, :8000], tau_s[:, :8000], T, n=20, track=True)
x0s = x0.copy(); x0s[2:] *= s; tau0 = np.abs(cs.uc(x0s)).max(1)
print(f"terminal set {xd.shape[1]} states | T_eps(1) {1e3*Te1:.2f} ms -> s = 0.8 T_eps(1)/T = {s:.3f} | T_eps(s) {1e3*Tes:.2f} ms (T_eps(1)/s = {1e3*Te1/s:.2f}) | "
      f"U_MAX {U_MAX} | s_U {s_U:.3f} | at s: lambda {cs.lam:.4f}, mu {cs.mu:.4f}, max|tau| terminal {np.round(np.abs(tau_s).max(1), 2)}, X_RA {np.round(tau0, 2)} | "
      f"eta {eta:.3e} (lambda delta {cs.lam*delta:.3e}), M2 {M2.max():.3f}, second-order T* {2e3*eta/M2.max():.2f} ms | URA scale {URA_FAC*s:.3f}")
json.dump(dict(T=T, MU0=MU0, DFRAC=DFRAC, T_eps1_ms=1e3 * Te1, s=s, T_eps_s_ms=1e3 * Tes, U_MAX=U_MAX.tolist(), s_U=s_U, URA_SCALE=URA_FAC * s, lam=cs.lam, mu=cs.mu, delta=delta,
               eta=eta, M2=float(M2.max()), Tstar2_ms=2e3 * eta / float(M2.max()), tau_terminal=np.abs(tau_s).max(1).tolist(), tau_xra=tau0.tolist()),
          open(os.path.join(HERE, "data", f"scale_{TASK}.json"), "w"), indent=1)
