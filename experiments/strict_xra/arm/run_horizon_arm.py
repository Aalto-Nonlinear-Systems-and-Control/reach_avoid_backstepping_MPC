"""Effect of the horizon at the same period: vanilla MPC and RA-MPC on a subset of the initial set for larger N."""
import pickle, time, sys, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config_arm as cf
NS = [int(a) for a in (sys.argv[1].split(",") if len(sys.argv) > 1 else ["267", "400"])]; N_SUB = int(sys.argv[2]) if len(sys.argv) > 2 else 40
if __name__ == "__main__":
    X0 = cf.initial_set()[:, :N_SUB]; n = X0.shape[1]; par = Parallel(n_jobs=-1, prefer="processes"); M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); omega = M["cfg"]["omega"]
    base_ra = np.array([L["status"] == "reached" for L in M["logs"]])[:n]; base_van = np.array([v["status"] == "reached" for v in M["van"]])[:n]; out = dict(rows=[dict(N=cf.N, ra=base_ra, van=base_van)])
    print(f"N = {cf.N} (main run, same subset): vanilla {base_van.sum()}/{n}, RA-MPC {base_ra.sum()}/{n}", flush=True)
    for N in NS:
        t0 = time.time(); van = par(delayed(cf.mpc_run)(X0[:, i], "target", cf.T, N) for i in range(n)); logs = par(delayed(cf.mpc_run)(X0[:, i], "cert", cf.T, N, omega) for i in range(n))
        ra = np.array([L["status"] == "reached" for L in logs]); vr = np.array([v["status"] == "reached" for v in van]); later = sum(L["status"] == "infeasible" and len(L["J"]) > 0 for L in logs)
        out["rows"].append(dict(N=N, ra=ra, van=vr, later=later, t_med=float(np.median([L["check"]["t_entry"] for L in logs if L["status"] == "reached"])) if ra.any() else np.nan, solve_ms=float(np.mean(np.concatenate([L["solve_ms"] for L in logs if len(L["solve_ms"])])))))
        print(f"N = {N} ({N*cf.T:.2f} s): vanilla {vr.sum()}/{n}, RA-MPC {ra.sum()}/{n} (later infeasible {later}), RA entry med {out['rows'][-1]['t_med']:.3f} s, {out['rows'][-1]['solve_ms']:.0f} ms/solve ({time.time()-t0:.0f}s)", flush=True)
        pickle.dump(out, open(f"{cf.DATA}/horizon.pkl", "wb"))
