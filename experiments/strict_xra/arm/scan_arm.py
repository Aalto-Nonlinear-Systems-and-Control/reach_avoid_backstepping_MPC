"""Scan mu_0 for the arm's base design: rest states inside X_RA, torques on the terminal set (gravity + s^2 part),
admissible scale s_U from the bounds, sampled one-step increase thresholds, second-order T*."""
import sys, json, numpy as np
sys.path.insert(0, ".")
from cert_arm import CertArm, fk, ik, phi_np, gravity
K1 = sys.argv[1] if len(sys.argv) > 1 else "k1/ARM_l0.02.txt"; DFRAC = 0.2; TAU_MAX = np.array([680.0, 500.0])
out = {}
for mu in [float(a) for a in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["0.1", "0.3", "1", "3"])]:
    c = CertArm(K1, mu=mu); M = 1 - c.shift; delta = DFRAC * M
    x, info = c.sample(delta, n_y=200000, seed=3); n = x.shape[1]
    y = info["y"]; rest_in = (info["kn"]**2 <= 2 * mu * (np.broadcast_to(c.psi(*y), (n,)) - c.shift)).mean()        # rest states (dq = 0) inside X_RA^0 at these positions
    rest_inT = (info["kn"] <= info["rad"]).mean()                                                                     # rest states inside X_RA^delta
    tau = c.uc(x); t1 = c.uc(x, with_gravity=False); G = gravity(x[:2])
    # s_U: |G_i + s^2 t1_i| <= tau_max_i  ->  s^2 <= (tau_max - G) / t1 if t1 > 0 else (tau_max + G) / (-t1)
    lim = np.where(t1 > 0, (TAU_MAX[:, None] - G) / np.maximum(t1, 1e-12), (TAU_MAX[:, None] + G) / np.maximum(-t1, 1e-12))
    s_U = float(np.sqrt(lim.min())); i = np.unravel_index(lim.argmin(), lim.shape)
    speeds = np.hypot(*(np.einsum("ijn,jn->in", __import__("cert_arm").jac(x[:2]), x[2:])))
    pd = c.Psidot(x, tau); eta = float(pd.min()); _, M2, _ = c.zoh(x, tau, 0.01, track=True); M2 = float(M2.max()); Ts1 = 2 * eta / M2
    P0 = c.Psi(x); eps = {}
    for T in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
        xn = c.zoh(x, tau, T, n=max(20, int(T / 0.001))); d = c.Psi(xn) - P0; keep = phi_np(*fk(xn[:2])) > 0
        eps[T] = float(d[keep].min())
    tmax = max([T for T, e in eps.items() if e > 0], default=0.0)
    print(f"mu_0 = {mu:<4g}: terminal set {n} states, speeds [{speeds.min():.3f}, {speeds.max():.3f}] m/s | rest states inside X_RA at these positions: {100*rest_in:.0f} % (inside X_RA^delta {100*rest_inT:.0f} %) | "
          f"tau at s=1: |tau1| <= {np.abs(tau[0]).max():.1f}, |tau2| <= {np.abs(tau[1]).max():.1f} (gravity part up to {np.abs(G).max(1).round(1)}) | s_U = {s_U:.2f} (channel {i[0]+1}) | "
          f"eta_1 {eta:.3e} (lambda delta {c.lam*delta:.3e}), M2 {M2:.3f}, T*_1 {1e3*Ts1:.1f} ms | one-step increase positive up to T = {tmax} s at s=1 | eps: " + ", ".join(f"{T}: {e:+.1e}" for T, e in eps.items()), flush=True)
    out[mu] = dict(n=n, rest_in=rest_in, s_U=s_U, eta=eta, M2=M2, Tstar1=Ts1, eps=eps, tau_max=[float(np.abs(tau[0]).max()), float(np.abs(tau[1]).max())])
json.dump({str(k): v for k, v in out.items()}, open("scan_arm.json", "w"), indent=1)
