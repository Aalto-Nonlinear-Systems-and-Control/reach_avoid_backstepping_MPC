"""Sampling-period sweep with the terminal design fixed at its base values (delta, omega) and N*T = 1.25 s.
Writes data/sweep.pkl."""
import pickle, time, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config as cf
from config import rm, C, R

SW = [(50, 0.025), (40, 0.03125), (25, 0.05), (20, 0.0625), (10, 0.125), (5, 0.25), (4, 0.3125), (2, 0.625)]

if __name__ == "__main__":
    X0 = cf.initial_set(); n = X0.shape[1]; par = Parallel(n_jobs=-1, prefer="processes")
    base = R.CertMPC(C, dt=cf.T, N=cf.N, dfrac=cf.DFRAC, u_max=np.asarray(cf.U_MAX)); omega = base.omega
    print(f"base design: delta = {cf.DELTA:.4f}, omega = {omega:.4g}, eps(T=0.05) = {base.eps:.3e}", flush=True)
    out = dict(SW=SW, omega=omega, rows=[])
    for N, T in SW:
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
        pickle.dump(out, open(f"{cf.DATA}/sweep.pkl", "wb"))
    print("saved data/sweep.pkl")
