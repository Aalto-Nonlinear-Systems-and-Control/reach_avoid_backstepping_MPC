"""Write the task spec (psi, target disc, sampling box) for one kitting task of the top-view arm, from ../arm/screen_racks.py.
Usage: python make_task.py red fwd  -> tasks/red_fwd.json"""
import os, sys, json, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "arm"))
import screen_racks as R
from screen_zigzag import analyse, y1, y2

name, d = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("red", "fwd")
k = R.NAMES.index(name); sgn = 1 if d == "fwd" else -1; B = R.BRIDGE; SUF = "_sq" if R.LAYOUT == "square" else ""
import pickle; _bp = os.path.join(HERE, "..", "arm", "data", f"racks_best{'_square' if SUF else ''}.pkl"); FIX_TOP = float(pickle.load(open(_bp, "rb"))["top"]) if os.path.exists(_bp) else 2.28
c = R.path(k, B); psi = R.certificate(c, sgn, 0.0); r = analyse(psi)
src = [-sgn * R.X0, R.YS[k] if sgn > 0 else R.YS[R.SIGMA[k]]]; dst = [sgn * R.X0, R.YS[R.SIGMA[k]] if sgn > 0 else R.YS[k]]
assert r["n_other"] == 0 and np.linalg.norm(np.array(r["tc"]) - dst) < 1e-6
safe = r["P"] > 0; box = [[float(R.Y1[safe].min()) - 0.1, float(R.Y1[safe].max()) + 0.1], [float(R.Y2[safe].min()) - 0.1, float(R.Y2[safe].max()) + 0.1]]
spec = dict(task=f"{name}_{d}{SUF}", layout=R.LAYOUT, psi=str(sp.expand(psi)), pmax=float(r["pmax"]), path=str(c), target=dst, source=src, RT=R.RT, h=R.H, box=box, bridge=B,
            X0=R.X0, YS=list(R.YS), SIGMA=list(R.SIGMA), color=R.COLORS[k], fixture=dict(w=R.FIX_W, y0=R.FIX_Y0, top=FIX_TOP), rack_d=R.RACK_D, view=[list(v) for v in R.VIEW],
            weight=dict(kappa=R.KAPPA, a=R.A_W, m=R.M), degree=int(r["deg"]))
json.dump(spec, open(os.path.join(HERE, "tasks", f"{name}_{d}{SUF}.json"), "w"), indent=1)
print(f"tasks/{name}_{d}{SUF}.json: degree {r['deg']}, pmax {r['pmax']:.4f}, source {src}, target {dst}, box {np.round(box, 2).tolist()}")
