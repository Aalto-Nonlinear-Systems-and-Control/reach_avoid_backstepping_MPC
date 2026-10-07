"""Three methods on the arm's initial set (at rest). Writes data/main.pkl."""
import pickle, time, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config_arm as cf
from config_arm import C, ra

if __name__ == "__main__":
    t0 = time.time(); X0 = cf.initial_set(); n = X0.shape[1]; P0 = C.Psi(X0); par = Parallel(n_jobs=-1, prefer="processes")
    ura = par(delayed(cf.ura_run)(X0[:, i]) for i in range(n)); uc = par(delayed(cf.ura_run)(X0[:, i], cf.T_MAX_URA, cf.SCALE) for i in range(n))
    zs, zt, zd = cf.uc_zoh(X0)
    ok_u = np.array([u["status"] == "reached" for u in ura]); up = np.array([u["upeak"] for u in ura]); viol = (up > np.asarray(cf.U_MAX)).any(1)
    print(f"set: {n} rest states, speeds 0 | terminal scale {cf.SCALE:g}, URA scale {cf.URA_SCALE:g} | URA {ok_u.sum()} {dict(Counter(u['status'] for u in ura))}, exceeding the bounds {viol.sum()} (among successes {(viol&ok_u).sum()}) | "
          f"terminal controller alone: continuous {sum(u['status']=='reached' for u in uc)}, sampled {(zs==1).sum()} | in X_RA {(P0>=0).sum()}, in terminal set {(P0>=cf.DELTA).sum()} ({time.time()-t0:.0f}s)", flush=True)
    van = par(delayed(cf.mpc_run)(X0[:, i], "target") for i in range(n)); fv = np.array([len(v["J"]) > 0 for v in van])
    print(f"vanilla: feasible at k=0 {fv.sum()}, closed loop {dict(Counter(v['status'] for v in van))} ({time.time()-t0:.0f}s)", flush=True)
    logs = par(delayed(cf.mpc_run)(X0[:, i], "cert") for i in range(n)); ra_ok = np.array([L["status"] == "reached" for L in logs])
    for L in logs + van:
        L["h_min"], L["entered"] = cf.safety_until_target(L["x"], L["u"], cf.T) if len(L["u"]) else (np.nan, False)
    print(f"RA-MPC: {ra_ok.sum()} {dict(Counter(L['status'] for L in logs))}, later infeasible {sum(L['status']=='infeasible' and len(L['J'])>0 for L in logs)}, "
          f"min psi along closed loops: RA {min(L['h_min'] for L in logs if L['status']=='reached'):.4f}, vanilla {min(v['h_min'] for v in van if v['status']=='reached'):.4f} ({time.time()-t0:.0f}s)", flush=True)
    pickle.dump(dict(X0=X0, Psi0=P0, ura=ura, uc=uc, zoh=dict(status=zs, t=zt, dmin=zd), van=van, logs=logs,
                     cfg=dict(T=cf.T, N=cf.N, U_MAX=cf.U_MAX, DELTA=cf.DELTA, M_PSI=cf.M_PSI, shift=C.shift, lam=C.lam, mu=C.mu, SCALE=cf.SCALE, URA_SCALE=cf.URA_SCALE,
                              omega=logs[0]["omega"], eps=logs[0]["eps"], SAFE_MARGIN=cf.SAFE_MARGIN)), open(f"{cf.DATA}/main.pkl", "wb"))
    print(f"saved data/main.pkl ({time.time()-t0:.0f}s)")
