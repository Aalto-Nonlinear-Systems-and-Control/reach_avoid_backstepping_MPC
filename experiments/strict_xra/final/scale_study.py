"""URA at a faster time scale than the bounds allow: initial set around the field of scale s_set, URA at s_ura,
terminal controller and MPC unchanged (scale 0.159, bounds (3, 3), N = 25, T = 0.05 s). Writes data/scale_study.pkl."""
import sys, pickle, time, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config as cf
from config import rm, Cert
from alt_set import build

if __name__ == "__main__":
    pairs = [tuple(float(a) for a in p.split(":")) for p in (sys.argv[1].split(",") if len(sys.argv) > 1 else ["3:3", "4:4", "1:4"])]
    par = Parallel(n_jobs=-1, prefer="processes"); res = {}
    for s_set, s_ura in pairs:
        t0 = time.time(); r0 = cf.R0 * s_set / cf.SCALE; X0 = build(r0, s_set); n = X0.shape[1]
        ura = par(delayed(cf.ura_run)(X0[:, i], 120.0, False, s_ura) for i in range(n))
        ok = np.array([u["status"] == "reached" for u in ura]); up = np.array([u["upeak"] for u in ura]); viol = (up > np.asarray(cf.U_MAX)).any(1); tv = np.array([u["t_viol"] for u in ura]); tt = np.array([u["t"] for u in ura])
        fv = np.array(par(delayed(lambda x: cf.vanilla_solve(x)[0])(X0[:, i]) for i in range(n)))
        logs = par(delayed(cf.ra_run)(X0[:, i]) for i in range(n)); ra = np.array([L["status"] == "reached" for L in logs]); tr = np.array([L["t"][-1] if L["status"] == "reached" else np.nan for L in logs])
        d2 = ((X0[:2] - rm.TARGET_POS[:, None])**2).sum(0); c3 = ok & fv & ra & viol
        print(f"set s={s_set:g} (r0={r0:.3f}, speeds [{X0[3].min():.2f}, {X0[3].max():.2f}]), URA s={s_ura:g}: URA {ok.sum()} {dict(Counter(u['status'] for u in ura))}, in its X_RA {sum(u['psi0']>=0 for u in ura)}, "
              f"exceeding (3,3): {viol.sum()} runs (omega {(up[:,0]>3).sum()}, a {(up[:,1]>3).sum()}), among successes {(viol&ok).sum()}, peak omega med {np.median(up[:,0]):.2f} max {up[:,0].max():.1f}, a med {np.median(up[:,1]):.2f} max {up[:,1].max():.2f}, "
              f"time above bounds med {np.median(tv[viol&ok]) if (viol&ok).any() else 0:.2f} s | URA time med {np.median(tt[ok]):.2f} s | vanilla {fv.sum()} | RA-MPC {ra.sum()} {dict(Counter(L['status'] for L in logs))} later inf {sum(L['status']=='infeasible' and len(L['J'])>0 for L in logs)}, time med {np.nanmedian(tr):.2f} s | "
              f"van&~RA {(fv&~ra).sum()} RA&~van {(ra&~fv).sum()} RA&~URA {(ra&~ok).sum()} URA&~RA {(ok&~ra).sum()} | all three {(ok&fv&ra).sum()}, with URA violation {c3.sum()} ({time.time()-t0:.0f}s)", flush=True)
        for i in np.where(c3)[0][np.argsort(-d2[c3])][:6]:
            print("   ", int(i), np.round(X0[:, i], 3).tolist(), "dist", round(float(np.sqrt(d2[i])), 2), "URA peaks", np.round(up[i], 2).tolist(), "t_URA", round(float(tt[i]), 2), "above bounds", round(float(tv[i]), 2), "s | t_RA", round(float(tr[i]), 2), flush=True)
        res[(s_set, s_ura)] = dict(X0=X0, ura=ura, fv=fv, ra=ra, tr=tr, r0=r0)
        pickle.dump(res, open(f"{cf.DATA}/scale_study.pkl", "wb"))
