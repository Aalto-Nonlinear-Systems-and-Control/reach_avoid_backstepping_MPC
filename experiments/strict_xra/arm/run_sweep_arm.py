"""Period sweep for the arm with N T = 0.2 s and the terminal design fixed (delta, omega); subset of the initial set."""
import pickle, time, sys, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config_arm as cf
from config_arm import C, ra

SW = [(200, 0.001), (133, 0.0015), (100, 0.002), (67, 0.003), (40, 0.005), (20, 0.01)]
N_SUB = int(sys.argv[1]) if len(sys.argv) > 1 else 40

if __name__ == "__main__":
    X0 = cf.initial_set()[:, :N_SUB]; n = X0.shape[1]; par = Parallel(n_jobs=-1, prefer="processes")
    M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); omega = M["cfg"]["omega"]; out = dict(SW=SW, omega=omega, idx=list(range(n)), rows=[])
    print(f"terminal design fixed: delta {cf.DELTA:.4f}, omega {omega:.4g}", flush=True)
    for N, T in SW:
        t0 = time.time()
        logs = par(delayed(cf.mpc_run)(X0[:, i], "cert", T, N, omega) for i in range(n)); van = par(delayed(cf.mpc_run)(X0[:, i], "target", T, N) for i in range(n))
        zs, _, _ = cf.uc_zoh(X0, dt=T)
        ra_ok = np.array([L["status"] == "reached" for L in logs]); fv = np.array([len(v["J"]) > 0 for v in van]); st = Counter(L["status"] for L in logs)
        later = sum(L["status"] == "infeasible" and len(L["J"]) > 0 for L in logs); dec = np.concatenate([np.diff(L["J"]) + ra.ELL_LOWER for L in logs if len(L["J"]) > 1] or [np.zeros(0)])
        row = dict(N=N, T=T, ra=ra_ok, fv=fv, van_reached=np.array([v["status"] == "reached" for v in van]), zoh=zs == 1, status=dict(st), inf0=int(sum(len(L["J"]) == 0 for L in logs)), later=later,
                   dec_viol=int((dec > 1e-6).sum()), steps=int(dec.size), eps=float(logs[0]["eps"]), t_med=float(np.median([L["check"]["t_entry"] for L in logs if L["status"] == "reached"])) if ra_ok.any() else np.nan,
                   solve_ms=float(np.mean(np.concatenate([L["solve_ms"] for L in logs if len(L["solve_ms"])]))))
        out["rows"].append(row)
        print(f"T={1e3*T:<4g} ms N={N:<3d} RA-MPC {ra_ok.sum():3d}/{n} (feasible at k=0 {n-row['inf0']}, later infeasible {later}) | vanilla {fv.sum():3d} (reached {row['van_reached'].sum()}) | u^c sampled {(zs==1).sum()} | "
              f"decrease viol {row['dec_viol']}/{row['steps']} | eps(T) {row['eps']:+.2e} | t_entry med {row['t_med']:.3f} s | {row['solve_ms']:.0f} ms/solve ({time.time()-t0:.0f}s)", flush=True)
        pickle.dump(out, open(f"{cf.DATA}/sweep.pkl", "wb"))
    print("saved data/sweep.pkl")
