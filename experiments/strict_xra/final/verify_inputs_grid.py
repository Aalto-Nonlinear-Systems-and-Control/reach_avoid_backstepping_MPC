"""Certified-by-construction bound on the terminal-controller inputs over X_RA^delta (Dubins), position by position:
for each y the velocity fibre is the disc |ydot - k1(y)| <= r(y), r = sqrt(2 mu (psi - shift - delta)); with
ydd = b(y) + M(y) ydot, b = mu grad psi - (lambda/2) k1, M = J_k1 + (lambda/2) I, and v >= |k1| - r on the fibre,
  |a| <= |b + M k1| + ||M|| r,   |omega| <= (|b + M k1| + ||M|| r) / (|k1| - r).
Evaluated on a fine position grid (no sampling in velocity); the grid gap is the only remaining approximation."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "design"))
from cert import Cert

K1_FILE = os.path.join(HERE, "..", "design", "k1", "W3_l0.06.txt"); MU0, DFRAC = 0.02, 0.2
U_MAX = np.array([3.0, 3.0]); NG = int(os.environ.get("NG", "2001"))
C = Cert(K1_FILE, mu=MU0); shift = C.shift; delta = DFRAC * (1 - shift)
g = np.linspace(-2.05, 2.05, NG); Y1, Y2 = np.meshgrid(g, g); n = Y1.size; y = np.array([Y1.ravel(), Y2.ravel()])
B = lambda a: np.broadcast_to(np.asarray(a, float), (n,))
psi = B(C.psi(*y)); h = B(C.h(*y)); phi = B(C.phi(*y))
on = (psi - shift - delta > 0) & (h > 0) & (phi > 0)              # positions of X_RA^delta outside the target
y = y[:, on]; n = y.shape[1]; psi = psi[on]
k = np.array([np.broadcast_to(np.asarray(a, float), (n,)) for a in C.k1(*y)])
gp = np.array([np.broadcast_to(np.asarray(a, float), (n,)) for a in C.gpsi(*y)])
Jm = C.J(*y); J = np.array([[np.broadcast_to(np.asarray(Jm[i][j], float), (n,)) for j in range(2)] for i in range(2)])
r = np.sqrt(2 * C.mu * (psi - shift - delta)); kn = np.hypot(k[0], k[1]); vmin = kn - r
assert vmin.min() > 0, "the velocity fibre contains v = 0 somewhere"
b = C.mu * gp - 0.5 * C.lam * k; M = J + 0.5 * C.lam * np.eye(2)[:, :, None]
c = b + np.einsum("ijn,jn->in", M, k); cn = np.hypot(c[0], c[1])
Mn = np.sqrt(((M**2).sum((0, 1)) + np.sqrt(((M**2).sum((0, 1)))**2 - 4 * (M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0])**2)) / 2)   # spectral norm of a 2x2
a_bd = cn + Mn * r; w_bd = a_bd / vmin
s_U = min(U_MAX[0] / w_bd.max(), np.sqrt(U_MAX[1] / a_bd.max()))
print(f"grid {NG}x{NG}: {n} positions of X_RA^delta; v_min over the set {vmin.min():.4f} | bounds at s = 1: |omega| <= {w_bd.max():.4f}, |a| <= {a_bd.max():.4f} "
      f"(sampled suprema 0.994, 0.237) -> certified-by-bound s_U = {s_U:.4f} (sampled 3.019)")
sj = os.path.join(HERE, "data", "scale.json")
if os.path.exists(sj):
    d = json.load(open(sj)); d.update(s_U_bound=float(s_U), omega1_bound=float(w_bd.max()), a1_bound=float(a_bd.max()), grid=NG); json.dump(d, open(sj, "w"), indent=1)
    print(f"chosen s = {d['s']:.4f} is {'within' if d['s'] <= s_U else 'OUTSIDE'} the certified-by-bound range")
