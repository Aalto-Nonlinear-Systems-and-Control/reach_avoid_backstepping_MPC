"""Arm, gravity compensated continuously: terminal-set size, rest states, s_eps(T) by bisection, URA from rest at
several scales with torque peaks (to choose the bounds)."""
import sys, numpy as np
sys.path.insert(0, ".")
from cert_arm import CertArm, fk, ik, phi_np, gravity, jac, f_np, _mul2
from scipy.integrate import solve_ivp
from joblib import Parallel, delayed
K1 = sys.argv[1] if len(sys.argv) > 1 else "k1/ARM_l0.04.txt"; T = 0.01
r = np.random.default_rng(5); Yb = np.array([r.uniform(-4, 8.5, 60000), r.uniform(-3, 3, 60000)])


def s_eps(c1, x1, T, lo=0.05, hi=8.0):
    P0 = c1.Psi(x1)
    def ok(s):
        c = c1.scaled(s); x = x1.copy(); x[2:] *= s; tau = c.uc(x)
        xn = c.zoh(x, tau, T, n=max(20, int(T / 0.0005)), comp_gravity=True); d = c.Psi(xn) - c.Psi(x); keep = phi_np(*fk(xn[:2])) > 0
        return d[keep].min() > 0
    if not ok(lo): return 0.0
    for _ in range(18):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def ura_rest(c, y0, T_max=60.0):
    q = ik(y0.reshape(2, 1))[:, 0]; x0 = np.array([q[0], q[1], 0.0, 0.0])
    def rhs(t, x): return f_np(x.reshape(4, 1), c.uc(x.reshape(4, 1)))[:, 0]
    e1 = lambda t, x: float(phi_np(*fk(x[:2]))); e1.terminal = True; e1.direction = -1
    e2 = lambda t, x: float(c.psi(*fk(x[:2]))); e2.terminal = True; e2.direction = -1
    e3 = lambda t, x: abs(np.sin(x[1])) - 0.1; e3.terminal = True; e3.direction = -1
    s = solve_ivp(rhs, (0, T_max), x0, events=[e1, e2, e3], rtol=1e-7, atol=1e-9, max_step=0.02)
    tau = c.uc(s.y); st = "reached" if s.t_events[0].size else "unsafe" if s.t_events[1].size else "singular" if s.t_events[2].size else "timeout"
    return st, float(s.t[-1]), np.abs(tau).max(1), float(c.Psi(x0.reshape(4, 1))[0])


for mu in (0.1, 0.3):
    for df in (0.3, 0.5):
        c1 = CertArm(K1, mu=mu); M = 1 - c1.shift; delta = df * M
        ps = np.broadcast_to(c1.psi(*Yb), (Yb.shape[1],)); out = (ps > c1.shift) & (phi_np(*Yb) > 0)
        share = ((ps > c1.shift + delta) & out).sum() / out.sum()
        x1, info = c1.sample(delta, n_y=120000, seed=3)
        se = s_eps(c1, x1, T)
        print(f"{K1[3:-4]} mu_0={mu} delta={df}M: terminal-set position share {100*share:.0f} % of the band outside X_T | s_eps(T={T}) = {se:.3f} (sampled; continuous gravity compensation)", flush=True)
    # URA from rest at band positions, scale s_URA; also the rest-state membership
    c1 = CertArm(K1, mu=mu); Y = Yb[:, (np.broadcast_to(c1.psi(*Yb), (Yb.shape[1],)) > c1.shift) & (phi_np(*Yb) > 0)][:, :150]
    for sU in (1.0, 2.0, 3.0):
        c = c1.scaled(sU); res = Parallel(n_jobs=-1)(delayed(ura_rest)(c, Y[:, j]) for j in range(Y.shape[1]))
        st = np.array([a[0] for a in res]); tt = np.array([a[1] for a in res]); pk = np.array([a[2] for a in res]); inX = np.array([a[3] >= 0 for a in res])
        ok = st == "reached"
        print(f"   URA at s={sU:g}, mu_0={mu}: from 150 rest positions: reached {ok.sum()}, unsafe {(st=='unsafe').sum()}, singular {(st=='singular').sum()}, timeout {(st=='timeout').sum()} | rest states in X_RA {inX.sum()} (reached among them {(ok&inX).sum()}) | "
              f"time med {np.median(tt[ok]) if ok.any() else np.nan:.1f} s | peak |tau1| med {np.median(pk[ok,0]) if ok.any() else np.nan:.0f} max {pk[:,0].max():.0f}, |tau2| med {np.median(pk[ok,1]) if ok.any() else np.nan:.0f} max {pk[:,1].max():.0f}", flush=True)
