"""URA at another time scale on the same initial set (for the figure baseline). Writes data/ura_s<scale>.pkl."""
import sys, pickle, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config as cf
from config import rm

if __name__ == "__main__":
    M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0 = M["X0"]; n = X0.shape[1]
    fv = M["fv"].astype(bool); ra = np.array([L["status"] == "reached" for L in M["logs"]]); d2 = ((X0[:2] - rm.TARGET_POS[:, None])**2).sum(0)
    for s in [float(a) for a in sys.argv[1].split(",")]:
        out = Parallel(n_jobs=-1, prefer="processes")(delayed(cf.ura_run)(X0[:, i], 300.0, False, s) for i in range(n))
        ok = np.array([u["status"] == "reached" for u in out]); up = np.array([u["upeak"] for u in out]); viol = (up > np.asarray(cf.U_MAX)).any(1); tt = np.array([u["t"] for u in out])
        print(f"URA at s={s:g}: {dict(Counter(u['status'] for u in out))} | in its X_RA {sum(u['psi0'] >= 0 for u in out)} | exceeding the bounds: {viol.sum()} runs (omega {(up[:,0]>3).sum()}, a {(up[:,1]>3).sum()}), among successes {(viol & ok).sum()} | "
              f"peak omega median {np.median(up[:,0]):.2f} max {up[:,0].max():.1f}, a max {up[:,1].max():.2f} | time med {np.median(tt[ok]) if ok.any() else np.nan:.1f} s | all three {(ok&fv&ra).sum()}, with violation {(ok&fv&ra&viol).sum()}, URA&~RA {(ok&~ra).sum()}, RA&~URA {(ra&~ok).sum()}", flush=True)
        cand = np.where(ok & fv & ra & viol)[0]
        for i in cand[np.argsort(-d2[cand])][:5]:
            print("   ", int(i), np.round(X0[:, i], 3).tolist(), "dist", round(float(np.sqrt(d2[i])), 2), "peaks", np.round(up[i], 2).tolist(), "t", round(float(tt[i]), 1), "time above bounds", round(float(out[i]["t_viol"]), 2), flush=True)
        pickle.dump(dict(s=s, out=out), open(f"{cf.DATA}/ura_s{s:g}.pkl", "wb"))
