"""Split M2(s) = c0 + c2 s^2 of the arm (c0: frozen gravity compensation, c2: the scalable part) and the resulting
largest one-step-increase period T_max = eta_1 / sqrt(c0 c2) at the best scale s* = sqrt(c0/c2)."""
import sys, numpy as np
sys.path.insert(0, ".")
from cert_arm import CertArm, fk, phi_np, jac, gravity
for K1 in ("k1/ARM_l0.02.txt", "k1/ARM_l0.04.txt"):
    for mu in (0.1, 0.3, 1.0):
        for df in (0.2, 0.5):
            c1 = CertArm(K1, mu=mu); M = 1 - c1.shift; delta = df * M
            x1, info = c1.sample(delta, n_y=120000, seed=3); n = x1.shape[1]
            rows = []
            for s in (0.1, 0.3, 1.0, 3.0):
                c = c1.scaled(s); x = x1.copy(); x[2:] *= s; tau = c.uc(x)
                pd = c.Psidot(x, tau); eta = float(pd.min()); _, M2, _ = c.zoh(x, tau, 0.002, n=20, track=True); rows.append((s, eta, float(M2.max()), float(np.quantile(M2, 0.99))))
            S = np.array([r[0] for r in rows]); M2s = np.array([r[2] for r in rows]); A = np.vstack([np.ones_like(S), S**2]).T
            c0, c2 = np.linalg.lstsq(A, M2s, rcond=None)[0]; eta1 = rows[2][1]
            Tmax = eta1 / np.sqrt(max(c0, 1e-12) * max(c2, 1e-12)); sstar = np.sqrt(max(c0, 1e-12) / max(c2, 1e-12))
            print(f"{K1[3:-4]} mu_0={mu:<4g} delta={df}M: eta_1 {eta1:.3e} | M2(s): " + ", ".join(f"s={r[0]}: {r[2]:.2f} (99% {r[3]:.2f}), eta {r[1]:.2e}" for r in rows) +
                  f" | fit c0 {c0:.2f}, c2 {c2:.2f} -> s* {sstar:.2f}, T_max {1e3*Tmax:.2f} ms", flush=True)
