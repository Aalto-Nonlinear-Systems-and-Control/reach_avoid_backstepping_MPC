"""Exploratory: initial velocities around the field of the unscaled controller (URA at s = 1), r0 scaled accordingly."""
import pickle, time, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config as cf
from config import rm, Cert

S_URA = 1.0; R0 = cf.R0 / cf.SCALE          # 0.189: same relative size as r0 = 0.03 for the scaled field

def build(r0, s, seed=cf.SEED_VEL):
    P = cf.positions(); n = P.shape[1]; r = np.random.default_rng(seed); c = Cert(cf.K1_FILE, mu=cf.MU_BASE).scaled(s)
    w = np.sqrt(r.random(n)) * np.exp(1j * r.uniform(0, 2 * np.pi, n))
    k = np.array([np.broadcast_to(np.asarray(a, float), (n,)) for a in c.k1(*P)]); yd = k + r0 * np.array([w.real, w.imag])
    return np.vstack([P, np.arctan2(yd[1], yd[0]), np.hypot(yd[0], yd[1])])

if __name__ == "__main__":
    t0 = time.time(); X0 = build(R0, S_URA); n = X0.shape[1]; par = Parallel(n_jobs=-1, prefer="processes")
    print(f"set around the s={S_URA} field, r0={R0:.3f}: speeds [{X0[3].min():.3f}, {X0[3].max():.3f}]", flush=True)
    ura = par(delayed(cf.ura_run)(X0[:, i], 300.0, False, S_URA) for i in range(n))
    ok = np.array([u["status"] == "reached" for u in ura]); up = np.array([u["upeak"] for u in ura]); viol = (up > np.asarray(cf.U_MAX)).any(1); tv = np.array([u["t_viol"] for u in ura])
    print(f"URA(s={S_URA}): {dict(Counter(u['status'] for u in ura))} | in its X_RA {sum(u['psi0'] >= 0 for u in ura)} | exceed the bounds: {viol.sum()} runs (omega {(up[:,0]>3).sum()}, a {(up[:,1]>3).sum()}), among successes {(viol&ok).sum()} | "
          f"peak omega median {np.median(up[:,0]):.2f} max {up[:,0].max():.1f}; a max {up[:,1].max():.2f} | time med {np.median([u['t'] for u in ura if u['status']=='reached']):.1f} s ({time.time()-t0:.0f}s)", flush=True)
    fv = np.array(par(delayed(lambda x: cf.vanilla_solve(x)[0])(X0[:, i]) for i in range(n)))
    logs = par(delayed(cf.ra_run)(X0[:, i]) for i in range(n)); ra = np.array([L["status"] == "reached" for L in logs])
    P0 = cf.C.Psi(X0)
    print(f"vanilla {fv.sum()} | RA-MPC {ra.sum()} {dict(Counter(L['status'] for L in logs))} (later infeasible {sum(L['status']=='infeasible' and len(L['J'])>0 for L in logs)}) | in R2 X_RA {(P0>=0).sum()} | "
          f"van&~RA {(fv&~ra).sum()} RA&~van {(ra&~fv).sum()} RA&~URA {(ra&~ok).sum()} URA&~RA {(ok&~ra).sum()} | all three {(ok&fv&ra).sum()}, with URA violation {(ok&fv&ra&viol).sum()} ({time.time()-t0:.0f}s)", flush=True)
    d2 = ((X0[:2] - rm.TARGET_POS[:, None])**2).sum(0); cand = np.where(ok & fv & ra & viol)[0]
    for i in cand[np.argsort(-d2[cand])][:10]:
        print("  ", int(i), np.round(X0[:, i], 3).tolist(), "dist", round(float(np.sqrt(d2[i])), 2), "URA peaks", np.round(up[i], 2).tolist(), "t_URA", round(float(ura[i]["t"]), 1), "time above bounds", round(float(tv[i]), 2), "t_RA", round(float(logs[i]["t"][-1]), 2), flush=True)
    pickle.dump(dict(X0=X0, ura=ura, fv=fv, logs=logs, R0=R0, S_URA=S_URA), open(f"{cf.DATA}/alt_set.pkl", "wb")); print("saved data/alt_set.pkl")
