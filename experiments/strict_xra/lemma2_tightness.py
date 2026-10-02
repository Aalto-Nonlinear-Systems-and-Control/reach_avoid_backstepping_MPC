"""Where does Lemma 2 lose? Constants on different domains + factor-by-factor comparison along ZOH arcs."""
import sys, importlib.util, numpy as np, sympy as sp
path = sys.argv[1] if len(sys.argv) > 1 else "out/rigorous_alpha2_mu0.1_vanilla_k1.py"
lam = float(sys.argv[2]) if len(sys.argv) > 2 else 0.005
T = 0.05
s = importlib.util.spec_from_file_location("m", path); m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
x1, x2, th, v = sp.symbols("x1 x2 th v"); X = [x1, x2, th, v]
Psi = sp.sympify(m.certificate_opt); U = [sp.sympify(e) for e in m.u_opt]
gp = [sp.diff(Psi, a) for a in X]
Pf = sp.lambdify(X, Psi, "numpy"); Uf = sp.lambdify(X, U, "numpy")
Lgf = sp.lambdify(X, [gp[2], gp[3]], "numpy")
GU = [sp.lambdify(X, [sp.diff(e, a) for a in X], "numpy") for e in U]
tgt = lambda x: x[1]**2 + ((x[0] + 1.7) / 0.5)**2 - 0.4
B = lambda f, x: np.array([np.broadcast_to(np.asarray(a, float), x.shape[1:]) for a in f(*x)])

def sample(lo, hi, n, seed=0):
    r = np.random.default_rng(seed); S = np.asarray(lo)[:, None] + (np.asarray(hi) - np.asarray(lo))[:, None] * r.random((4, n))
    p = np.asarray(Pf(*S), float); return S, p

def consts(S, sel):
    x = S[:, sel]; u = B(Uf, x); Lg = B(Lgf, x)
    Mg = np.abs(Lg).max(1); lu = np.array([np.sqrt((B(g, x)**2).sum(0)).max() for g in GU])
    Mv = np.sqrt(x[3]**2 + (u**2).sum(0)).max()
    return Mg, lu, Mv

lsys = 1.0  # |df/dx| <= max(1, |v|) and |v| <= 1 on these domains
def Tstar(Mg, lu, Mv, delta):
    a = (Mg * lu).sum() * Mv; return np.log(1 + lsys * lam * delta / a) / lsys, a

doms = {
    "full state, theta in [-pi,pi], v in [0.1,1.5]  (cap of X_RA)": ([-2.2, -2.2, -np.pi, 0.1], [2.2, 2.2, np.pi, 1.5]),
    "synthesis box, theta in [2pi/3,4pi/3], v in [0.1,1]  (used so far)": ([-2.2, -2.2, 2 * np.pi / 3, 0.1], [2.2, 2.2, 4 * np.pi / 3, 1.0]),
}
print(f"certificate: {path}, lambda = {lam}")
for name, (lo, hi) in doms.items():
    S, p = sample(lo, hi, 600000)
    for lab, sel in (("{Psi>=0}", p >= 0), ("{Psi>=0.1 maxPsi}", p >= 0.1 * p.max())):
        Mg, lu, Mv = consts(S, sel)
        Ts, a = Tstar(Mg, lu, Mv, 0.1 * p.max())
        print(f"  {name}\n     on {lab:18s}: n={sel.sum():6d}  M_g={np.round(Mg,3)}  l_uc={np.round(lu,2)}  M_v={Mv:.2f}  "
              f"sum M_g l_uc * M_v = {a:.1f}  ->  T*(0.1 maxPsi) = {Ts:.2e} s")

# factor-by-factor comparison along ZOH arcs, starting in X_RA \ X_T (synthesis box)
S, p = sample(*doms["synthesis box, theta in [2pi/3,4pi/3], v in [0.1,1]  (used so far)"], 400000, seed=1)
sel = (p >= 0) & (tgt(S) > 0); xk = S[:, sel][:, :20000]; uk = B(Uf, xk)
Mg, lu, Mv = consts(S, p >= 0)
n = 50; h = T / n; x = xk.copy()
drift = np.zeros(xk.shape[1]); emax = np.zeros(xk.shape[1]); lipr = np.zeros(xk.shape[1])
rhs = lambda z: np.array([z[3] * np.cos(z[2]), z[3] * np.sin(z[2]), uk[0], uk[1]])
for _ in range(n):
    a1 = rhs(x); a2 = rhs(x + h / 2 * a1); a3 = rhs(x + h / 2 * a2); a4 = rhs(x + h * a3); x = x + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
    d = np.linalg.norm(x - xk, axis=0); drift = np.maximum(drift, d)
    ut = B(Uf, x); Lg = B(Lgf, x)
    e = np.abs((Lg * (ut - uk)).sum(0)); emax = np.maximum(emax, e)            # true ZOH mismatch in Psidot
    lipr = np.maximum(lipr, np.linalg.norm(ut - uk, axis=0) / np.maximum(d, 1e-12))
xdot = np.sqrt(xk[3]**2 + (uk**2).sum(0))
drift_bound = (np.exp(lsys * T) - 1) / lsys * xdot
beta_k = (Mg * lu).sum() * drift_bound            # beta_k(T) = gamma(T) |xdot_k|
beta = (Mg * lu).sum() * (np.exp(lsys * T) - 1) / lsys * Mv
pk = np.asarray(Pf(*xk), float)
q = lambda a: f"median {np.median(a):.3g}, 99% {np.quantile(a, 0.99):.3g}, max {a.max():.3g}"
print(f"\nalong ZOH arcs of length T = {T} s from {xk.shape[1]} states in X_RA \\ X_T:")
print(f"  drift: bound / actual            {q(drift_bound / np.maximum(drift, 1e-12))}")
print(f"  local Lipschitz |du|/|dx| on arc {q(lipr)}   vs  l_uc = {np.round(lu,1)}")
print(f"  Psidot mismatch: beta_k / actual {q(beta_k / np.maximum(emax, 1e-12))}")
print(f"  Psidot mismatch: beta   / actual {q(beta / np.maximum(emax, 1e-12))}")
print(f"  actual mismatch max_t e(t)       {q(emax)}   vs beta(T) = {beta:.3g}")
print(f"  margin lambda*Psi_k               {q(lam * pk)}")
print(f"  states with actual mismatch < lambda*Psi_k: {(emax < lam * pk).mean()*100:.1f}%   "
      f"(Lemma 2 needs beta(T) < lambda*Psi_k: {(beta < lam * pk).mean()*100:.1f}%)")
