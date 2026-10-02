"""Rigour checks for the strict X_RA candidates produced by run_strict.m."""
import sys, os, importlib.util, numpy as np, sympy as sp
from scipy.integrate import solve_ivp
S = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(S, "out")
cases = [int(a) for a in sys.argv[1:]] or [1, 2, 3]

y1, y2 = sp.symbols("y1 y2")
phi_y = y2**2 + ((y1 + 1.7) / 0.5)**2 - 0.4
psi_y = sp.Rational(1, 1000) * (-phi_y + 300) * (-(y1**4 + y2**4 - 16) * (y1**4 + y2**4 - 4))
x1, x2, th, v = sp.symbols("x1 x2 th v"); X = [x1, x2, th, v]
f = sp.Matrix([v * sp.cos(th), v * sp.sin(th), 0, 0]); G = sp.Matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
tgt = lambda x: x[1]**2 + ((x[0] + 1.7) / 0.5)**2 - 0.4
safe = sp.lambdify((y1, y2), psi_y, "numpy")

def load(path):
    spec = importlib.util.spec_from_file_location(os.path.basename(path)[:-3], path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

g = np.linspace(-2.2, 2.2, 1501); Y1, Y2 = np.meshgrid(g, g)
PSI = safe(Y1, Y2); PHI = Y2**2 + ((Y1 + 1.7) / 0.5)**2 - 0.4
dA = (g[1] - g[0])**2

for i in cases:
    txt = open(os.path.join(OUT, f"strict_{i}_vanilla_k1.txt")).read().split("\n")
    k1 = [sp.sympify(t.replace("^", "**")) for t in txt[:2]]; lam, dl = float(txt[2]), float(txt[3])
    print(f"\n===== case {i}: lambda = {lam:g}, vanilla delta/lambda = {dl/lam:.5f}  (far-side local max psi = 9.0742)")
    # 1) k1-level ECGBF condition  grad(psi).k1 >= lam*psi - delta  on {psi>=0, phi>=0}
    c = sp.lambdify((y1, y2), sp.diff(psi_y, y1) * k1[0] + sp.diff(psi_y, y2) * k1[1] - lam * psi_y + dl, "numpy")
    Cv = c(Y1, Y2); reg = (PSI >= 0) & (PHI >= 0)
    print(f"  k1 condition on {{psi>=0, phi>=0}}: min = {Cv[reg].min():.3e} (relative to lam*max psi: {Cv[reg].min()/(lam*PSI.max()):.2e}), "
          f"frac < 0: {(Cv[reg] < 0).mean():.2e}")
    print(f"  vanilla certified set {{psi >= delta/lambda}} outside target: area {((PSI >= dl/lam) & (PHI > 0)).sum()*dA:.3f} of safe-set area {(PSI >= 0).sum()*dA:.3f}")
    # 2) full-state certificate from the SOP stage
    m = load(os.path.join(OUT, f"tmp_strict_{i}_result.py"))
    Psi = sp.sympify(m.certificate_opt); U = [sp.sympify(e) for e in m.u_opt]
    gp = sp.Matrix([sp.diff(Psi, a) for a in X])
    Pf = sp.lambdify(X, Psi, "numpy"); Uf = sp.lambdify(X, U, "numpy")
    Pdf = sp.lambdify(X, (gp.T * (f + G * sp.Matrix(U)))[0], "numpy")
    rng = np.random.default_rng(0)
    lo = np.array([-2.2, -2.2, 2 * np.pi / 3, 0.1]); hi = np.array([2.2, 2.2, 4 * np.pi / 3, 1.0])
    Sx = lo[:, None] + (hi - lo)[:, None] * rng.random((4, 400000))
    P = np.broadcast_to(np.asarray(Pf(*Sx), float), Sx.shape[1:]); T_ = tgt(Sx)
    ra = P >= 0; out = ra & (T_ > 0)
    Pd = np.broadcast_to(np.asarray(Pdf(*Sx[:, out]), float), (out.sum(),))
    u = np.array([np.broadcast_to(np.asarray(a, float), Sx.shape[1:]) for a in Uf(*Sx[:, ra])])
    print(f"  full certificate: {{Psi>=0}} holds {ra.mean()*100:.2f}% of synthesis box (v in [0.1,1]); max Psi = {P.max():.3f}")
    viol = Pd - lam * P[out]
    print(f"  ECGBF Psidot >= lam*Psi on {{Psi>=0}}\\X_T: min (Psidot - lam Psi) = {viol.min():.3e}, frac < 0: {(viol < 0).mean():.2e}, "
          f"median Psidot/Psi = {np.median(Pd/P[out]):.3e}")
    print(f"  u^c on {{Psi>=0}}: max|omega| = {np.abs(u[0]).max():.2f}, max|a| = {np.abs(u[1]).max():.2f}, within |u|<=5: {(np.abs(u) <= 5).all(0).mean()*100:.1f}%")
    # MPC candidate set (same seed / box as the notebooks)
    np.random.seed(42); lo2 = np.array([-2, -2, 2*np.pi/3, -1.0]); hi2 = np.array([2, 2, 4*np.pi/3, 1.0])
    C = np.random.rand(4, 1000) * (hi2 - lo2)[:, None] + lo2[:, None]
    cand = np.where((safe(C[0], C[1]) >= 0) & (tgt(C) > 0))[0][:300]
    PC = np.broadcast_to(np.asarray(Pf(*C[:, cand]), float), cand.shape)
    print(f"  MPC candidates already in X_RA = {{Psi>=0}}: {(PC >= 0).sum()}/300")
    # 3) continuous closed loop from X_RA
    idx = np.where(out)[0][:40]; ok = 0; mono = 0; times = []
    def rhs(t, x):
        uu = np.array(Uf(*x), float).flatten(); return [x[3]*np.cos(x[2]), x[3]*np.sin(x[2]), uu[0], uu[1]]
    ev = lambda t, x: tgt(x); ev.terminal = True; ev.direction = -1
    for j in idx:
        s = solve_ivp(rhs, (0, 60), Sx[:, j], events=ev, rtol=1e-8, atol=1e-10, max_step=0.01)
        p = np.asarray(Pf(*s.y), float); ok += s.status == 1; mono += np.all(np.diff(p) >= -1e-9)
        if s.status == 1: times.append(s.t[-1])
    print(f"  continuous u^c from {len(idx)} states in X_RA\\X_T: reached {ok}, Psi nondecreasing in {mono}, "
          f"reach time median {np.median(times) if times else float('nan'):.2f}s max {max(times) if times else float('nan'):.2f}s")
    # 4) ZOH one step and Lemma 2 constants
    xs = Sx[:, out][:, :60000]; us = np.array([np.broadcast_to(np.asarray(a, float), xs.shape[1:]) for a in Uf(*xs)])
    def zoh(x, u, T, n=50):
        h = T / n
        for _ in range(n):
            r = lambda z: np.array([z[3]*np.cos(z[2]), z[3]*np.sin(z[2]), u[0], u[1]])
            k1_ = r(x); k2_ = r(x + h/2*k1_); k3_ = r(x + h/2*k2_); k4_ = r(x + h*k3_); x = x + h/6*(k1_ + 2*k2_ + 2*k3_ + k4_)
        return x
    p0 = np.asarray(Pf(*xs), float)
    for T in (0.05, 0.01):
        d = np.asarray(Pf(*zoh(xs, us, T)), float) - p0
        print(f"  ZOH one step T={T}: frac dPsi < 0 = {(d < 0).mean():.3f}, frac dPsi < (e^(lam T)-1)Psi = {(d < (np.exp(lam*T)-1)*p0).mean():.3f}")
    gpn = [sp.lambdify(X, a, "numpy") for a in gp]
    Lg = [np.broadcast_to(np.asarray(gpn[2](*xs), float), p0.shape), np.broadcast_to(np.asarray(gpn[3](*xs), float), p0.shape)]
    Mg = [np.abs(a).max() * 1.1 for a in Lg]
    gu = [[sp.lambdify(X, sp.diff(e, a), "numpy") for a in X] for e in U]
    lu = [np.sqrt(sum(np.broadcast_to(np.asarray(fn(*xs), float), p0.shape)**2 for fn in row)).max() * 1.1 for row in gu]
    Mv = np.sqrt((xs[3])**2 + us[0]**2 + us[1]**2).max() * 1.1; lsys = 1.1
    a = (Mg[0]*lu[0] + Mg[1]*lu[1]) * Mv
    for dfrac in (0.1, 0.5):
        delta = dfrac * P.max()
        Tst = np.log(1 + lsys * lam * delta / a) / lsys
        print(f"  Lemma 2: M_g={np.round(Mg,3)}, l_uc={np.round(lu,2)}, M_v={Mv:.2f};  delta={dfrac}*maxPsi -> T*(delta) = {Tst:.3e} s")
