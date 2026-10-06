"""Three methods on the one initial set (final setting) with full RA-MPC logs. Writes data/main.pkl."""
import pickle, time, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config as cf
from config import rm, C

if __name__ == "__main__":
    t0 = time.time(); X0 = cf.initial_set(); n = X0.shape[1]; P0 = C.Psi(X0)
    par = Parallel(n_jobs=-1, prefer="processes")
    ura = par(delayed(cf.ura_run)(X0[:, i], cf.T_MAX_URA, False, cf.URA_SCALE) for i in range(n))      # URA baseline
    uc = ura if abs(cf.URA_SCALE - cf.SCALE) < 1e-12 else par(delayed(cf.ura_run)(X0[:, i]) for i in range(n))   # terminal controller alone
    zs, zt, zd = cf.ura_zoh(X0)
    fv = np.array(par(delayed(lambda x: cf.vanilla_solve(x)[0])(X0[:, i]) for i in range(n)))
    logs = par(delayed(cf.ra_run)(X0[:, i]) for i in range(n))
    for L in logs:      # smallest h of the continuous closed loop up to the first entry into X_T
        L["h_min"], L["entered"] = cf.safety_until_target(L["x"], L["u"], cf.T) if len(L["u"]) else (np.nan, False)
    ok_u = np.array([u["status"] == "reached" for u in ura]); ok_r = np.array([L["status"] == "reached" for L in logs])
    up = np.array([u["upeak"] for u in ura]); viol = (up > np.asarray(cf.U_MAX)).any(1)
    print(f"set scale {cf.SET_SCALE:g} (r0 {cf.R0:.3f}, speeds [{X0[3].min():.3f}, {X0[3].max():.3f}]), URA scale {cf.URA_SCALE:g}: URA exceeding the bounds {viol.sum()} (among successes {(viol & ok_u).sum()}), "
          f"terminal controller alone reaches {sum(u['status'] == 'reached' for u in uc)}", flush=True)
    print(f"n={n} | in X_RA {(P0>=0).sum()}, in terminal set {(P0>=cf.DELTA).sum()} | URA {ok_u.sum()} {dict(Counter(u['status'] for u in ura))} | "
          f"sampled URA {(zs==1).sum()} | vanilla {fv.sum()} | RA-MPC {ok_r.sum()} {dict(Counter(L['status'] for L in logs))} ({time.time()-t0:.0f}s)", flush=True)
    # vanilla receding-horizon rollouts on the vanilla-feasible states (closed-loop check + solve times)
    van = par(delayed(cf.vanilla_rollout)(X0[:, i]) for i in np.where(fv)[0])
    for v in van:
        v["h_min"], v["entered"] = cf.safety_until_target(v["x"], v["u"], cf.T) if len(v["u"]) else (np.nan, False)
    RLs = [L for L in logs if L["status"] == "reached"]
    print("vanilla closed loop on its feasible states:", dict(Counter(v["status"] for v in van)), "| min h: RA-MPC",
          min(L["h_min"] for L in RLs), "vanilla", min(v["h_min"] for v in van if v["status"] == "reached"), flush=True)
    pickle.dump(dict(X0=X0, Psi0=P0, ura=ura, uc=uc, zoh=dict(status=zs, t=zt, dmin=zd), fv=fv, logs=logs, van=dict(zip(np.where(fv)[0].tolist(), van)),
                     cfg=dict(T=cf.T, N=cf.N, U_MAX=cf.U_MAX, R0=cf.R0, DELTA=cf.DELTA, M_PSI=cf.M_PSI, shift=C.shift, lam=C.lam, mu=C.mu,
                              omega=logs[0]["omega"], eps=logs[0]["eps"], H_SUB=cf.H_SUB, SAFE_MARGIN=cf.SAFE_MARGIN, SET_SCALE=cf.SET_SCALE, URA_SCALE=cf.URA_SCALE)), open(f"{cf.DATA}/main.pkl", "wb"))
    print(f"saved data/main.pkl ({time.time()-t0:.0f}s)")
