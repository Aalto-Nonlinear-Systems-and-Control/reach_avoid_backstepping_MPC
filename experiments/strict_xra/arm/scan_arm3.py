"""Option A (held torque): largest period with positive one-step increase of Psi on X_RA^delta as a function of the
scale s, and URA torque peaks from rest at larger scales."""
import sys, json, numpy as np
sys.path.insert(0, ".")
from cert_arm import CertArm, fk, ik, phi_np, f_np
from scipy.integrate import solve_ivp
from joblib import Parallel, delayed
K1 = "k1/ARM_l0.04.txt"


def T_eps(c1, x1, s, lo=2e-4, hi=0.05):
    c = c1.scaled(s); x = x1.copy(); x[2:] *= s; tau = c.uc(x); P0 = c.Psi(x)
    def ok(T):
        xn = c.zoh(x, tau, T, n=max(10, int(T / 0.0002))); d = c.Psi(xn) - P0; keep = phi_np(*fk(xn[:2])) > 0
        return d[keep].min() > 0
    if not ok(lo): return 0.0
    for _ in range(16):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


out = {}
for mu in (0.1,):
    for df in (0.3, 0.5):
        c1 = CertArm(K1, mu=mu); delta = df * (1 - c1.shift); x1, _ = c1.sample(delta, n_y=200000, seed=3)
        row = {s: 1e3 * T_eps(c1, x1, s) for s in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)}
        out[f"mu{mu}_d{df}"] = row
        print(f"mu_0={mu}, delta={df}M ({x1.shape[1]} states): largest period with Psi increasing on every state [ms]: " + ", ".join(f"s={s:g}: {t:.2f}" for s, t in row.items()), flush=True)
# URA from rest at band positions: torque peaks for larger scales
r = np.random.default_rng(5); Yb = np.array([r.uniform(-4, 8.5, 60000), r.uniform(-3, 3, 60000)])
c1 = CertArm(K1, mu=0.1); Y = Yb[:, (np.broadcast_to(c1.psi(*Yb), (Yb.shape[1],)) > c1.shift) & (phi_np(*Yb) > 0)][:, :150]
def ura_rest(c, y0, T_max=60.0):
    q = ik(y0.reshape(2, 1))[:, 0]; x0 = np.array([q[0], q[1], 0.0, 0.0])
    rhs = lambda t, x: f_np(x.reshape(4, 1), c.uc(x.reshape(4, 1)))[:, 0]
    e1 = lambda t, x: float(phi_np(*fk(x[:2]))); e1.terminal = True; e1.direction = -1
    e2 = lambda t, x: float(c.psi(*fk(x[:2]))); e2.terminal = True; e2.direction = -1
    s = solve_ivp(rhs, (0, T_max), x0, events=[e1, e2], rtol=1e-7, atol=1e-9, max_step=0.01)
    return ("reached" if s.t_events[0].size else "unsafe"), float(s.t[-1]), np.abs(c.uc(s.y)).max(1)
for sU in (3.0, 4.0, 5.0, 6.0):
    res = Parallel(n_jobs=-1)(delayed(ura_rest)(c1.scaled(sU), Y[:, j]) for j in range(150)); ok = np.array([a[0] == "reached" for a in res]); pk = np.array([a[2] for a in res]); tt = np.array([a[1] for a in res])
    line = f"URA s={sU:g}: reached {ok.sum()}/150, time med {np.median(tt[ok]):.2f} s, peaks |tau1| med {np.median(pk[ok,0]):.0f} max {pk[ok,0].max():.0f}, |tau2| med {np.median(pk[ok,1]):.0f} max {pk[ok,1].max():.0f}"
    for B in ((120, 40), (150, 50), (200, 60)):
        v = (pk > np.array(B)).any(1); line += f" | exceeding {B}: {v.sum()} ({(v&ok).sum()} of successes)"
    print(line, flush=True)
json.dump(out, open("scan_arm3.json", "w"), indent=1)
