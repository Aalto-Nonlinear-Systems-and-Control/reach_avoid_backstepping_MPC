"""Extends the period sweep of run_sweep.py to coarser periods, with N T = 1.25 s kept, so that the loss of the
guarantee beyond the sampling-period bound shows in the success rate. New rows are merged into data/sweep.pkl
(the file before the extension is kept as data/sweep_base.pkl); rows already present are not recomputed."""
import os, shutil, pickle, time, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config as cf
from config import rm

EXT = [(3, 1.25 / 3), (1, 1.25)]

if __name__ == "__main__":
    src, base = f"{cf.DATA}/sweep.pkl", f"{cf.DATA}/sweep_base.pkl"
    if not os.path.exists(base): shutil.copy(src, base)
    out = pickle.load(open(src, "rb")); have = {round(r["T"], 5) for r in out["rows"]}; omega = out["omega"]
    X0 = cf.initial_set(); n = X0.shape[1]; par = Parallel(n_jobs=-1, prefer="processes")
    print(f"extending sweep: omega = {omega:.4g}, have T = {sorted(have)}", flush=True)
    for N, T in EXT:
        if round(T, 5) in have: continue
        t0 = time.time()
        logs = par(delayed(cf.ra_run)(X0[:, i], dt=T, N=N, omega=omega) for i in range(n))
        fv = np.array(par(delayed(lambda x, T=T, N=N: cf.vanilla_solve(x, dt=T, N=N)[0])(X0[:, i]) for i in range(n)))
        zs, zt, zd = cf.ura_zoh(X0, dt=T)
        ra = np.array([L["status"] == "reached" for L in logs]); st = Counter(L["status"] for L in logs)
        inf0 = np.array([L["status"] == "infeasible" and len(L["J"]) == 0 for L in logs]); later = np.array([L["status"] == "infeasible" and len(L["J"]) > 0 for L in logs])
        dec = np.concatenate([np.diff(L["J"]) + rm.ELL_LOWER for L in logs if len(L["J"]) > 1] or [np.zeros(0)])
        row = dict(N=N, T=T, ra=ra, fv=fv, zoh=zs == 1, status=dict(st), inf0=int(inf0.sum()), later=int(later.sum()), unsafe=st.get("unsafe", 0), timeout=st.get("timeout", 0),
                   dec_viol=int((dec > 1e-6).sum()), steps=int(dec.size), bound_ok=int(sum(L["check"]["bound_ok"] for L in logs if L["status"] == "reached")),
                   eps=float(logs[0]["eps"]), t_med=float(np.median([L["check"]["t_entry"] for L in logs if L["status"] == "reached"])) if ra.any() else np.nan,
                   solve_ms=float(np.mean(np.concatenate([L["solve_ms"] for L in logs if len(L["solve_ms"])]))))
        out["rows"].append(row)
        print(f"T={T:<7g} N={N:<3d} RA-MPC {ra.sum():3d}/{n} (feasible at k=0: {n-inf0.sum()}, later infeasible {later.sum()}, unsafe {row['unsafe']}, timeout {row['timeout']}) | "
              f"vanilla {fv.sum():3d} | sampled URA {(zs==1).sum():3d} | decrease viol {row['dec_viol']}/{row['steps']} | eps(T) {row['eps']:+.2e} | "
              f"t_entry med {row['t_med']:.2f} s | {row['solve_ms']:.0f} ms/solve ({time.time()-t0:.0f}s)", flush=True)
        out["rows"].sort(key=lambda r: r["T"]); out["SW"] = [(r["N"], r["T"]) for r in out["rows"]]
        pickle.dump(out, open(src, "wb"))
    print("saved data/sweep.pkl")
