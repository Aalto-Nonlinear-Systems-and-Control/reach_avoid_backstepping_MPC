"""Baseline: continuous reach-avoid controller u^c alone (no MPC) from the 300 initial states and from X_RA."""
import numpy as np
from scipy.integrate import solve_ivp
from cert import Cert
import ra_mpc_cert as R, ra_mpc_dubins as rm
from eval_mpc import S, cand

def run(c, x0, T_max):
    def rhs(t, x):
        u = c.uc(x.reshape(4, 1))[:, 0]
        return [x[3] * np.cos(x[2]), x[3] * np.sin(x[2]), u[0], u[1]]
    ev_t = lambda t, x: float(rm.tgt_np(*x)); ev_t.terminal = True; ev_t.direction = -1
    ev_s = lambda t, x: float(rm.safe_np(*x)); ev_s.terminal = True; ev_s.direction = -1
    ev_v = lambda t, x: abs(x[3]) - 1e-4; ev_v.terminal = True; ev_v.direction = -1
    try:
        s = solve_ivp(rhs, (0, T_max), x0, events=[ev_t, ev_s, ev_v], rtol=1e-7, atol=1e-9, max_step=0.05)
    except Exception:
        return "error", np.nan, np.zeros(2)
    u = np.abs(c.uc(s.y)).max(1) if s.y.shape[1] else np.zeros(2)
    if s.t_events[0].size: return "reached", s.t[-1], u
    if s.t_events[1].size: return "unsafe", s.t[-1], u
    if s.t_events[2].size: return "v=0", s.t[-1], u
    return "timeout", s.t[-1], u

c1 = Cert("k1/W3_l0.06.txt", mu=0.02)
for name, c, Tm in (("R2 (s=0.159)", c1.scaled(0.159), 400.0), ("unscaled (s=1)", c1, 100.0)):
    P = c.Psi(S[:, cand]); inside = P >= 0
    out = [run(c, S[:, i], Tm) for i in cand]
    st = [o[0] for o in out]
    from collections import Counter
    print(f"{name}: initial states inside X_RA^0: {inside.sum()}/300;  u^c alone from the 300 states: {dict(Counter(st))}")
    # from states inside the terminal set (intrinsic samples)
    x, _ = c.sample(0.2 * (1 - c.shift), n_y=4000, seed=5); x = x[:, :60]
    o2 = [run(c, x[:, j], Tm) for j in range(x.shape[1])]
    tt = np.array([o[1] for o in o2 if o[0] == "reached"]); uu = np.array([o[2] for o in o2 if o[0] == "reached"])
    print(f"   from 60 states of {{Psi >= 0.2 M}}: {dict(Counter(o[0] for o in o2))}; reach time median {np.median(tt):.1f} s, max {tt.max():.1f} s; max|omega| {uu[:,0].max():.2f}, max|a| {uu[:,1].max():.3f}")
