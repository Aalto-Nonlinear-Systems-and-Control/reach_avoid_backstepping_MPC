"""Base design checks for the Panda certificate at s = 1 and the choice of mu0."""
import time, numpy as np, mujoco
from cert_franka import CertFranka, sample, phi_xy, ik_batch, BOX, Z0, CZ, CO, TARGET, ALPHA
from franka_model import MJ, NQ
c = CertFranka(mu=0.001); print(f"lambda {c.lam}, delta_s {c.ds:.2e}, shift {c.shift:.2e}")
# first-layer margin on a grid (planar part, z = Z0, y4 = 0) and the separable bound
g1 = np.linspace(*BOX[0], 401); g2 = np.linspace(*BOX[1], 201); A, B = np.meshgrid(g1, g2); Y = np.vstack([A.ravel(), B.ravel(), Z0 * np.ones(A.size), np.zeros(A.size)])
P = c.psi_of_y(Y); K = c.k1_of_y(Y); m = (P >= 0) & (phi_xy(Y) > 0)
import casadi as ca; gfun = ca.Function("g", [ca.SX.sym("y", 4)], [ca.gradient(c.psi_y(ca.SX.sym("y", 4)), ca.SX.sym("y", 4))]) if False else None
ys = ca.SX.sym("y", 4); G = np.array(ca.Function("g", [ys], [ca.gradient(c.psi_y(ys), ys)]).map(Y.shape[1])(Y))
marg = (G * K).sum(0) - c.lam * P + c.ds; print(f"ECGBF margin on the planar grid (in X_S minus X_T): min {marg[m].min():+.3e}; max |k1_xy| {np.abs(K[:2, m]).max():.3f} m/s")
# rest states of the corridor in X_RA for several mu0, and joint-limit margins of the sampled states
r = np.random.default_rng(5); Yr = np.array([r.uniform(*BOX[0], 40000), r.uniform(*BOX[1], 40000), Z0 * np.ones(40000), np.zeros(40000)])
Yr = Yr[:, (c.psi_of_y(Yr) > 0) & (phi_xy(Yr) > 0)][:, :3000]; Qr, er = ik_batch(Yr); Xr = np.vstack([Qr, np.zeros_like(Qr)])
lo, hi = MJ.jnt_range[:, 0], MJ.jnt_range[:, 1]
for mu in (0.0003, 0.001, 0.003, 0.01):
    cm = CertFranka(mu=mu); t0 = time.time(); x, _ = sample(cm, 0.5 * (1 - cm.shift), n_y=60000, seed=11); u = cm.uc(x)
    jm = np.min(np.r_[x[:4] - lo[:, None], hi[:, None] - x[:4]]); print(f"mu0 {mu:<7g}: rest states in X_RA {100*(cm.Psi(Xr) >= 0).mean():5.1f}% | terminal-set samples {x.shape[1]}, joint-limit margin {jm:+.3f} rad | max|tau_1| {np.round(np.abs(u).max(1), 3)} ({time.time()-t0:.1f}s)", flush=True)
