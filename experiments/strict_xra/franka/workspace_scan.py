"""Choose the length scale alpha (metres per layout unit) and the tool height z0 for the square kitting layout:
inverse kinematics with the hand pointing down (y4 = 0) on the corridors of all four tasks; joint-limit margins and the
smallest singular value of the output Jacobian."""
import os, sys, itertools, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); os.environ["LAYOUT"] = "square"; sys.path.insert(0, os.path.join(HERE, "..", "arm"))
import screen_racks as R
from screen_zigzag import y1, y2
from franka_model import ik, J_f, MJ

# corridor points (layout units) of the four tubes, from the forward certificates
pts = []
for k in range(4):
    c = R.path(k, R.BRIDGE); pf = sp.lambdify((y1, y2), R.certificate(c, 1, 0.0)); g1 = np.linspace(-3.6, 3.6, 145); g2 = np.linspace(1.0, 6.0, 101); A, B = np.meshgrid(g1, g2)
    P = pf(A, B); pts.append(np.c_[A[P > 0], B[P > 0]])
pts = np.unique(np.round(np.vstack(pts), 3), axis=0)[::3]; print(f"{len(pts)} corridor points of the four tubes")
lo = MJ.jnt_range[:, 0]; hi = MJ.jnt_range[:, 1]
for alpha, z0 in itertools.product((0.09, 0.10, 0.11), (0.05, 0.10, 0.15, 0.20, 0.25, 0.30)):
    ok, marg, smin, fails = 0, np.inf, np.inf, 0
    for p in pts:
        y = np.array([alpha * p[0], alpha * p[1], z0, 0.0]); X, Y = y[1], -y[0]
        q0 = np.array([np.arctan2(Y, X), 0.3, -2.0, 2.3]); qq, err = ik(y, q0)
        if err > 1e-8 or not np.all(np.isfinite(qq)): fails += 1; continue
        m = np.min(np.r_[qq - lo, hi - qq]); s = np.linalg.svd(np.array(J_f(qq)), compute_uv=False).min()
        marg, smin = min(marg, m), min(smin, s); ok += m > 0
    print(f"alpha {alpha:.2f} m/unit, z0 {z0:.2f} m: IK fails {fails}, within joint limits {ok}/{len(pts)}, smallest joint-limit margin {marg:+.3f} rad, smallest singular value {smin:.3f}", flush=True)
