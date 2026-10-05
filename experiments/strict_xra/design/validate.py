"""Check cert.py formulas against the MATLAB export (weight (2.5-y1)^2/600, lambda=0.005, mu=0.1)."""
import importlib.util, numpy as np, sympy as sp
from cert import Cert, H, y1, y2
c = Cert("../out/alpha_d8_l0.005_k1.txt", mu=0.1, psi=(sp.Rational(5, 2) - y1)**2 * H / 600)
s = importlib.util.spec_from_file_location("m", "../out/rigorous_alpha2_mu0.1_vanilla_k1.py"); m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
x1, x2, th, v = sp.symbols("x1 x2 th v"); X = [x1, x2, th, v]
Pf = sp.lambdify(X, sp.sympify(m.certificate_opt), "numpy"); Uf = sp.lambdify(X, [sp.sympify(e) for e in m.u_opt], "numpy")
r = np.random.default_rng(0); x = np.array([r.uniform(-2, 2, 5000), r.uniform(-2, 2, 5000), r.uniform(-np.pi, np.pi, 5000), r.uniform(0.1, 1.5, 5000)])
dP = np.abs(c.Psi(x) - Pf(*x)); u_ref = np.array(Uf(*x), float); dU = np.abs(c.uc(x) - u_ref)
print(f"max |Psi - Psi_matlab| = {dP.max():.2e} (scale {np.abs(Pf(*x)).max():.1f});  max |u - u_matlab| = {dU.max():.2e} (scale {np.abs(u_ref).max():.1f})")
