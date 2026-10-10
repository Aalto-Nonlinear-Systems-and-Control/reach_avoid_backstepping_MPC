"""Rate margin on the terminal set: Psidot under u^c minus lambda Psi, and the excess grad psi . k1 - lambda psi, A against LP designs."""
import os, sys, json, numpy as np
sys.path.insert(0, "/Users/jianqiang/Downloads/reach_avoid_backstepping_MPC/experiments/strict_xra/arm_top"); os.chdir(sys.path[0])
from cert_top import CertArm, fk, phi_np, TASK
SC = json.load(open(f"data/scale_{TASK}.json")); s = float(sys.argv[1]); C0 = CertArm(f"k1/{TASK}_d8_l0.1.txt", mu=SC["MU0"])
designs = {"A (12)+scale": C0.scaled(s)}
for tag in sys.argv[2:]:
    d = json.load(open(f"data/{tag}.json")); designs[tag.split("/")[-1]] = CertArm(f"k1/{d['k1_file']}", mu=d["mu"])
for name, c in designs.items():
    x, info = c.sample(0.5 * (1 - c.shift), n_y=200000, seed=3); x = x[:, phi_np(*fk(x[:2])) > 0][:, :20000]; y = fk(x[:2])
    pd = c.Psidot(x, c.uc(x)); P = c.Psi(x); g = np.array(c.gpsi(*y)); k = np.array([np.broadcast_to(np.asarray(a, float), y[0].shape) for a in c.k1(*y)])
    exc = (g * k).sum(0) - c.lam * np.broadcast_to(c.psi(*y), y[0].shape)
    print(f"{name:34s} eta = min Psidot {pd.min():.4f} (median {np.median(pd):.3f}) | min (Psidot - lam Psi) {np.min(pd - c.lam * P):+.4f} | excess grad psi.k1 - lam psi: min {exc.min():+.4f}, median {np.median(exc):.4f}")
