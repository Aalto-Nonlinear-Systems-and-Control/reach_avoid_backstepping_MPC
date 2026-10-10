"""Where the one-step increase of Psi fails at a given period, and the variation rate M2 of Psidot during the hold."""
import os, sys, json, numpy as np
sys.path.insert(0, "/Users/jianqiang/Downloads/reach_avoid_backstepping_MPC/experiments/strict_xra/arm_top"); os.chdir(sys.path[0])
from cert_top import CertArm, fk, phi_np, TASK, TARGET_POS
SC = json.load(open(f"data/scale_{TASK}.json")); s = float(sys.argv[1]); Tp = float(sys.argv[2]); C0 = CertArm(f"k1/{TASK}_d8_l0.1.txt", mu=SC["MU0"])
designs = {"A (12)+scale": C0.scaled(s)}
for tag in sys.argv[3:]:
    d = json.load(open(f"data/{tag}.json")); designs[tag.split("/")[-1]] = CertArm(f"k1/{d['k1_file']}", mu=d["mu"])
for name, c in designs.items():
    x, _ = c.sample(0.5 * (1 - c.shift), n_y=200000, seed=7); x = x[:, phi_np(*fk(x[:2])) > 0][:, :8000]; y = fk(x[:2])
    tau = c.uc(x); xn, M2, pmin = c.zoh(x, tau, Tp, n=max(20, int(Tp / 0.0005)), track=True); d = c.Psi(xn) - c.Psi(x); ok = phi_np(*fk(xn[:2])) > 0; bad = ok & (d <= 0)
    J = np.array(c.J1(*y), dtype=object); Jn = np.sqrt(sum(np.broadcast_to(np.asarray(J[i][j], float), y[0].shape)**2 for i in range(2) for j in range(2)))
    dist = np.hypot(y[0] - TARGET_POS[0], y[1] - TARGET_POS[1]); spd = np.hypot(*np.array([np.broadcast_to(np.asarray(a, float), y[0].shape) for a in c.k1(*y)]))
    print(f"{name:32s} T {1e3*Tp:.1f} ms: failing {bad.sum()}/{ok.sum()} | M2 median {np.median(M2):.2f}, p99 {np.percentile(M2, 99):.2f}, max {M2.max():.2f} | |J_k1| median {np.median(Jn):.3f}, max {Jn.max():.3f}"
          + (f" | failures at distance to target {np.median(dist[bad]):.2f} (all {np.median(dist):.2f}), |J_k1| there {np.median(Jn[bad]):.3f}, |k1| {np.median(spd[bad]):.2f}" if bad.any() else ""))
