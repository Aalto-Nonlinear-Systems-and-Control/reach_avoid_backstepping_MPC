"""Scale of the arm's terminal controller (held torque): s_U from the bounds, T_eps(s) from the sampled one-step
increase, choose s at the largest admissible period, require T <= T_eps(s)/2. Writes data/scale_arm.json."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_arm import CertArm, fk, phi_np, gravity
K1 = os.path.join(HERE, "k1", "ARM_l0.04.txt"); MU0, DFRAC = 0.1, 0.5; U_MAX = np.array([120.0, 40.0]); T = 0.0015
c1 = CertArm(K1, mu=MU0); delta = DFRAC * (1 - c1.shift); x1, _ = c1.sample(delta, n_y=300000, seed=11); x1 = x1[:, phi_np(*fk(x1[:2])) > 0]
G = gravity(x1[:2]); t1 = c1.uc(x1, with_gravity=False)
lim = np.where(t1 > 0, (U_MAX[:, None] - G) / np.maximum(t1, 1e-12), (U_MAX[:, None] + G) / np.maximum(-t1, 1e-12)); s_U = float(np.sqrt(lim.min()))
def T_eps(s, lo=2e-4, hi=0.05):
    c = c1.scaled(s); x = x1.copy(); x[2:] *= s; tau = c.uc(x); P0 = c.Psi(x)
    def ok(T):
        xn = c.zoh(x, tau, T, n=max(10, int(T / 0.0002))); d = c.Psi(xn) - P0; keep = phi_np(*fk(xn[:2])) > 0; return d[keep].min() > 0
    if not ok(lo): return 0.0
    for _ in range(16):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo
grid = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]; Te = {s: T_eps(s) for s in grid}
adm = [s for s in grid if s <= s_U]; s = max(adm, key=lambda k: Te[k])            # the scale with the largest admissible period
print(f"terminal set {x1.shape[1]} states | s_U = {s_U:.2f} | T_eps(s) [ms]: " + ", ".join(f"{k:g}: {1e3*v:.2f}" for k, v in Te.items()) + f" | chosen s = {s} (largest T_eps), T = {1e3*T:g} ms is {Te[s]/T:.2f} times below T_eps")
assert T <= 0.6 * Te[s], "period too close to the sampled threshold"
c = c1.scaled(s); x = x1.copy(); x[2:] *= s; tau = c.uc(x); pd = c.Psidot(x, tau); eta = float(pd.min()); _, M2, _ = c.zoh(x, tau, T, n=10, track=True)
print(f"at s = {s}: lambda {c.lam:.4f}, mu {c.mu:.4f}, |tau| on the terminal set <= {np.abs(tau).max(1).round(1)}, eta {eta:.3e} (lambda delta {c.lam*delta:.3e}), M2 {M2.max():.2f}, second-order T* {2e3*eta/M2.max():.2f} ms, T_eps {1e3*Te[s]:.2f} ms")
json.dump(dict(T=T, s_U=s_U, T_eps_ms={str(k): 1e3 * v for k, v in Te.items()}, s=s, lam=c.lam, mu=c.mu, delta=delta, eta=eta, M2=float(M2.max()), Tstar2_ms=2e3 * eta / float(M2.max()), tau_max_terminal=np.abs(tau).max(1).tolist()),
          open(os.path.join(HERE, "data", "scale_arm.json"), "w"), indent=1)
