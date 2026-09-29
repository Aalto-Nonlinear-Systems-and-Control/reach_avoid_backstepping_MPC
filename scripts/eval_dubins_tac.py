"""Batch evaluation of the TAC-revision RA-MPC (Dubins). Saves data/tac_dubins_eval.npz.

Env: OMEGAS (comma list), DELTA, N_CAND, PLANT (nominal|zoh), T_MAX.
"""
import os, sys, time
import numpy as np
from joblib import Parallel, delayed
from collections import Counter
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))
import ra_mpc_dubins as rm

OMEGAS = [float(a) for a in os.environ.get("OMEGAS", "10,100,1000").split(",")]
DELTA = float(os.environ.get("DELTA", "0.1"))
N_CAND = int(os.environ.get("N_CAND", "300"))
PLANT = os.environ.get("PLANT", "nominal")
T_MAX = float(os.environ.get("T_MAX", "10"))

np.random.seed(42)
lo = np.array([-2, -2, 2 * np.pi / 3, -1.0]); hi = np.array([2, 2, 4 * np.pi / 3, 1.0])
S = np.random.rand(4, 1000) * (hi - lo)[:, None] + lo[:, None]
cand = np.where((rm.safe_np(*S) >= 0) & (rm.tgt_np(*S) > 0))[0][:N_CAND]

out = []
for om in OMEGAS:
    t0 = time.time()
    res = Parallel(n_jobs=-1, prefer="processes")(
        delayed(rm.run_one)(S[:, i], 0.05, 25, DELTA, om, T_MAX, PLANT) for i in cand)
    for i, r in zip(cand, res):
        r["idx"] = int(i)
    st = Counter(r["status"] for r in res)
    reached = [r for r in res if r["status"] == "reached"]
    viol_runs = sum(r["check"]["n_decrease_viol"] > 0 for r in res)
    dec = np.concatenate([np.diff(r["J"]) + rm.ELL_LOWER for r in reached if len(r["J"]) > 1])
    print(f"[{PLANT}, T_max={T_MAX:g}] omega={om:g}: {dict(st)} | reached {len(reached)}/{len(res)} "
          f"| decrease violations in reached runs {(dec > 1e-6).sum()}/{dec.size} steps, runs with any violation {viol_runs} "
          f"| t_entry<=t_bound {sum(r['check']['bound_ok'] for r in reached)}/{len(reached)} "
          f"| u^c in bounds {np.mean(np.concatenate([r['uc_in_bounds'] for r in res if len(r['uc_in_bounds'])])):.3f} "
          f"| cheaper candidate used in {np.mean(np.concatenate([r['used_candidate'] for r in res if len(r['used_candidate'])])):.3f} of steps "
          f"| solve {np.mean(np.concatenate([r['solve_ms'] for r in res if len(r['solve_ms'])])):.1f} ms  ({time.time()-t0:.0f}s)",
          flush=True)
    out.append(res)

np.savez(os.path.join(ROOT, "data", f"tac_dubins_eval_{PLANT}.npz"), omegas=np.array(OMEGAS),
         delta=DELTA, cand=cand, S=S, t_max=T_MAX, results=np.array(out, dtype=object), allow_pickle=True)
print(f"saved data/tac_dubins_eval_{PLANT}.npz")
