"""Study: k1 from a sampled linear program instead of the SOS program (12), red_fwd task, lambda fixed.
    max gamma  s.t.  grad psi(y_i) . k1(y_i) - lambda psi(y_i) >= gamma   at samples of X_S minus X_T   (delta_s = 0)
                     |k1_j(y_l)| <= K                                   at samples of {psi >= 0}
k1_j in monomials of y/4 up to degree D, coefficients bounded. Exchange loop: solve on the samples, evaluate on a grid
10 times finer, add every fine point with margin below gamma/2, repeat. Final check on a grid 20 times finer than the start."""
import os, sys, time, numpy as np, sympy as sp
from scipy.optimize import linprog
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_top import CertArm, phi_np, BOX, TASK, SPEC, y1, y2
LAM, K = 0.05, 1.0; D = int(sys.argv[1]) if len(sys.argv) > 1 else 8; H0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
psi = sp.sympify(SPEC["psi"]) / SPEC["pmax"]; P = sp.lambdify((y1, y2), psi); G = sp.lambdify((y1, y2), [sp.diff(psi, v) for v in (y1, y2)])
mons = [(i, j) for i in range(D + 1) for j in range(D + 1 - i)]; nb = len(mons)
basis = lambda A, B: np.stack([(A / 4)**i * (B / 4)**j for i, j in mons], axis=-1)


def grid(h):
    A, B = np.meshgrid(np.arange(BOX[0][0], BOX[0][1] + h, h), np.arange(BOX[1][0], BOX[1][1] + h, h)); A, B = A.ravel(), B.ravel(); return A, B, P(A, B), phi_np(A, B)


def margin(c, A, B):
    Mb = basis(A, B); g = np.array(G(A, B)); k = np.stack([Mb @ c[:nb], Mb @ c[nb:]]); return (g * k).sum(0) - LAM * P(A, B), np.abs(k).max(0)


def solve(La, Lb, Ka, Kb):
    Mb = basis(La, Lb); g = np.array(G(La, Lb)); Mk = basis(Ka, Kb); Z = np.zeros_like(Mk); z = np.zeros((len(Ka), 1))
    Rl = np.hstack([-(g[0][:, None] * Mb), -(g[1][:, None] * Mb), np.ones((len(La), 1))])
    Rk = np.vstack([np.hstack([Mk, Z, z]), np.hstack([-Mk, Z, z]), np.hstack([Z, Mk, z]), np.hstack([Z, -Mk, z])])
    res = linprog(np.r_[np.zeros(2 * nb), -1.0], A_ub=np.vstack([Rl, Rk]), b_ub=np.r_[-LAM * P(La, Lb), K * np.ones(4 * len(Ka))],
                  bounds=[(-50, 50)] * (2 * nb) + [(None, 1.0)], method="highs")
    if res.x is None: raise RuntimeError(res.message)
    return res.x[:-1], res.x[-1]


t0 = time.time(); A, B, p, f = grid(H0); L = (p >= 0) & (f >= 0); Kk = p >= -0.02; La, Lb, Ka, Kb = A[L], B[L], A[Kk], B[Kk]
Af, Bf, pf, ff = grid(H0 / 10); mf = (pf >= 0) & (ff >= 0); Af, Bf = Af[mf], Bf[mf]
for it in range(30):
    c, gam = solve(La, Lb, Ka, Kb); mg, kn = margin(c, Af, Bf); bad = mg < gam / 2
    i = int(np.argmin(mg)); print(f"D={D} iter {it}: {len(La)} Lie samples | gamma {gam:+.3e} | fine grid min margin {mg.min():+.3e} at y = ({Af[i]:.2f}, {Bf[i]:.2f}), psi {P(Af[i], Bf[i]):.3f} | points below gamma/2: {bad.sum()} | max|k1| {kn.max():.3f}", flush=True)
    if mg.min() >= 0 and kn.max() <= K * 1.02: break
    La, Lb = np.r_[La, Af[bad]], np.r_[Lb, Bf[bad]]; Ka, Kb = np.r_[Ka, Af[kn > K]], np.r_[Kb, Bf[kn > K]]
A3, B3, p3, f3 = grid(H0 / 20); m3 = (p3 >= 0) & (f3 >= 0); mg3, kn3 = margin(c, A3[m3], B3[m3])
print(f"final check on spacing {H0/20:g} ({m3.sum()} points): min margin {mg3.min():+.3e}, max|k1| {kn3.max():.3f} | total {time.time()-t0:.0f}s")
c_sos = CertArm(os.path.join(HERE, "k1", f"{TASK}_d8_l0.05.txt"), mu=1.0); kk = c_sos.k1(A3[m3], B3[m3]); gg = np.array(G(A3[m3], B3[m3]))
print(f"SOS k1 (degree 8) on the same grid: min margin {((gg[0]*kk[0] + gg[1]*kk[1]) - LAM*p3[m3]).min():+.3e}")
np.save(os.path.join(HERE, "data", f"lp_k1_{TASK}_D{D}.npy"), c)
