"""Scale and torque bounds of the Panda terminal controller (exact time scaling): T_eps(1) from the sampled one-step
increase of Psi on the terminal set at s = 1, s = 0.8 T_eps(1)/T, task bounds U = 1.25 s^2 max|tau_1| on the terminal set
(rounded up, checked against the Panda's hardware limits 87/87/87/12 N m), URA at 2 s. Writes data/scale_<TASK>.json."""
import os, json, numpy as np
from cert_franka import CertFranka, sample, zoh, phi_xy, TASK, HERE
from franka_model import Y_f
MU0, DFRAC, T, URA_FAC = 0.001, 0.5, 0.01, 2.0; HW = np.array([87.0, 87.0, 87.0, 12.0])
c1 = CertFranka(mu=MU0); delta = DFRAC * (1 - c1.shift); xd, Yd = sample(c1, delta, n=40000, seed=11); print(f"terminal-set samples {xd.shape[1]}", flush=True)


def T_eps(c, x, lo=1e-3, hi=2.0):
    tau = c.uc(x); P0 = c.Psi(x)
    def ok(Tp):
        xn = zoh(x, tau, Tp); d = c.Psi(xn) - P0; k = phi_xy(c.y(xn)) > 0; return (d[k].min() > 0) if k.any() else True
    if not ok(lo): return 0.0
    for _ in range(14):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


x = xd; Te1 = T_eps(c1, x); s = 0.8 * Te1 / T; t1 = np.abs(c1.uc(xd)).max(1)
U = np.ceil(1.25 * s * s * t1 * 10) / 10; cs = c1.scaled(s); xs = xd.copy(); xs[4:] *= s; Tes = T_eps(cs, xs, lo=1e-4, hi=0.5)
taus = cs.uc(xs); s_U = s * float(np.sqrt((U[:, None] / np.maximum(np.abs(taus), 1e-12)).min())); x0, _ = sample(c1, 0.0, n=20000, seed=12); x0[4:] *= s
print(f"T_eps(1) {1e3*Te1:.2f} ms -> s = {s:.3f} | T_eps(s) {1e3*Tes:.2f} ms (T_eps(1)/s = {1e3*Te1/s:.2f}) | task bounds U {U} (hardware {HW}, ratio {np.round(U/HW, 3)}) | s_U {s_U:.3f} | "
      f"terminal controller on the terminal set {np.round(np.abs(taus).max(1), 2)}, on X_RA {np.round(np.abs(cs.uc(x0)).max(1), 2)} | lambda_s {cs.lam:.3f}, mu_s {cs.mu:.4f} | URA scale {URA_FAC*s:.3f}")
os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
json.dump(dict(T=T, MU0=MU0, DFRAC=DFRAC, T_eps1_ms=1e3 * Te1, s=s, T_eps_s_ms=1e3 * Tes, U_MAX=U.tolist(), HW=HW.tolist(), s_U=s_U, URA_SCALE=URA_FAC * s, lam=cs.lam, mu=cs.mu, delta=delta,
               tau_terminal=np.abs(taus).max(1).tolist()), open(os.path.join(HERE, "data", f"scale_{TASK}.json"), "w"), indent=1)
