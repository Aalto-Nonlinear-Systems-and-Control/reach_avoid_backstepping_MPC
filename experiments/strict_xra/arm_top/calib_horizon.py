"""Horizon calibration on a subset: vanilla and RA-MPC success for several N at T = 10 ms."""
import sys, time, pickle, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config_top as cf
NS = [int(a) for a in sys.argv[1].split(",")]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
X0 = cf.initial_set()[:, :n]; par = Parallel(n_jobs=-1, prefer="processes"); out = {}
for N in NS:
    t0 = time.time(); van = par(delayed(cf.mpc_run)(X0[:, i], "target", cf.T, N) for i in range(n)); ra = par(delayed(cf.mpc_run)(X0[:, i], "cert", cf.T, N) for i in range(n))
    vr = np.array([v["status"] == "reached" for v in van]); rr = np.array([L["status"] == "reached" for L in ra]); out[N] = dict(van=vr, ra=rr, van_st=Counter(v["status"] for v in van), ra_st=Counter(L["status"] for L in ra))
    print(f"N={N} ({N*cf.T:.2f} s): vanilla {vr.sum()}/{n} {dict(out[N]['van_st'])} | RA-MPC {rr.sum()}/{n} {dict(out[N]['ra_st'])} | van&~RA {(vr&~rr).sum()} | solve ms RA {np.mean(np.concatenate([L['solve_ms'] for L in ra if len(L['solve_ms'])])):.0f} ({time.time()-t0:.0f}s)", flush=True)
    pickle.dump(out, open(f"{cf.DATA}/calib_horizon.pkl", "wb"))
