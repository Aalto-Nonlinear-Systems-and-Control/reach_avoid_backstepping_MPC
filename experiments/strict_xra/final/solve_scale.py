"""Choose the time scale s of the terminal controller from the base design (Section 7 of TAC_EXPERIMENT_UPDATE.md):
s_U from the input bounds on X_RA^delta (7.6), s_eps(T) from the sampled one-step increase of Psi at the period T,
and, for reference only, the conservative s_T of the second-order bound. Writes data/scale.json."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); DESIGN = os.path.abspath(os.path.join(HERE, "..", "design"))
sys.path.insert(0, DESIGN); sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "python")))
from cert import Cert
import ra_mpc_dubins as rm

K1_FILE = os.path.join(DESIGN, "k1", "W3_l0.06.txt"); MU0, DFRAC = 0.02, 0.2
U_MAX = np.array([3.0, 3.0]); T = float(os.environ.get("T_MPC", "0.05")); SAFETY = float(os.environ.get("SAFETY", "0.8"))

C1 = Cert(K1_FILE, mu=MU0)                       # base design, s = 1
M_PSI = 1.0 - C1.shift; DELTA = DFRAC * M_PSI
x, info = C1.sample(DELTA, n_y=400000, seed=11); n = x.shape[1]
tgt = rm.tgt_np(*x) > 0; x = x[:, tgt]; n = x.shape[1]
u1 = C1.uc(x)                                   # (omega_1, a_1) at s = 1 on X_RA^delta (outside the target)

# ── s_U: smallest admissible root over the set, Dubins closed form (7.7) ──
rho = np.minimum(U_MAX[0] / np.abs(u1[0]), np.sqrt(U_MAX[1] / np.abs(u1[1])))
s_U = float(rho.min()); i_U = int(rho.argmin())
print(f"X_RA^delta sample: {n} states, speeds at s=1 in [{x[3].min():.3f}, {x[3].max():.3f}] | sup|omega_1| {np.abs(u1[0]).max():.4f}, sup|a_1| {np.abs(u1[1]).max():.4f}")
print(f"s_U = {s_U:.4f} (binding channel: {'omega' if U_MAX[0]/abs(u1[0, i_U]) < np.sqrt(U_MAX[1]/abs(u1[1, i_U])) else 'a'})", flush=True)

# ── s_eps(T): one-step increase of Psi under ZOH at the base scale over an arc of length tau = s T ──
P0 = C1.Psi(x)
def eps1(tau, n_sub=None):
    xn = C1.zoh(x, u1, tau, n=max(40, int(tau / 0.0025)))
    d = C1.Psi(xn) - P0; out = rm.tgt_np(*xn) > 0          # arcs that enter the target are not required to increase
    return float(d[out].min()), float((d[out] > 0).mean())
taus = np.geomspace(0.005, 0.5, 60); e = np.array([eps1(t)[0] for t in taus])
ok = e > 0; tau_max = float(taus[ok].max()) if ok.any() else 0.0
# refine the crossing by bisection
lo, hi = tau_max, float(taus[~ok].min()) if (~ok).any() else tau_max * 2
for _ in range(20):
    mid = 0.5 * (lo + hi); (lo, hi) = (mid, hi) if eps1(mid)[0] > 0 else (lo, mid)
tau_eps = lo; s_eps = tau_eps / T
print(f"one-step increase at the base scale positive up to arc length tau_eps = {tau_eps*1e3:.1f} ms  ->  s_eps(T={T}) = {s_eps:.4f}", flush=True)

# ── reference: conservative second-order bound ──
pd = C1.Psidot(x, u1); eta1 = float(pd.min()); _, M2, _ = C1.zoh(x, u1, 0.05, track=True); M21 = float(M2.max())
Tstar1 = 2 * eta1 / M21; s_T = Tstar1 / T
print(f"reference (conservative): eta_1 {eta1:.4e} (lambda*delta {C1.lam*DELTA:.4e}), M_2,1 {M21:.3f}, T*_1 {Tstar1*1e3:.2f} ms, s_T = {s_T:.4f}", flush=True)

# ── choice and verification at the final scale ──
s = min(s_U, SAFETY * s_eps); C = C1.scaled(s)
xs = x.copy(); xs[3] *= s; us = C.uc(xs)                   # the same states at scale s
d_fin, frac = (lambda r: (r[0], r[1]))(eps1(s * T))
print(f"chosen s = {s:.4f} = min(s_U, {SAFETY} * s_eps) | at s: lambda {C.lam:.4f}, mu {C.mu:.4e}, max|k1_i| on the ring {C.margin_on_grid(801)['kmax']:.3f} m/s")
print(f"verification at s, T={T}: sup|omega| {np.abs(us[0]).max():.3f} <= {U_MAX[0]}, sup|a| {np.abs(us[1]).max():.4f} <= {U_MAX[1]}; smallest one-step increase {d_fin:+.3e} (increasing {100*frac:.2f} %)", flush=True)
json.dump(dict(T=T, safety=SAFETY, s_U=s_U, s_eps=s_eps, tau_eps=tau_eps, s_T_conservative=s_T, eta1=eta1, M21=M21, s=s, lam=C.lam, mu=C.mu, delta=DELTA,
               sup_omega1=float(np.abs(u1[0]).max()), sup_a1=float(np.abs(u1[1]).max()), eps_at_s=d_fin, n_samples=n), open(os.path.join(HERE, "data", "scale.json"), "w"), indent=1)
print("saved data/scale.json")
