"""Input-bound study for setting R2 (certificate fixed): feasibility, saturation, entry time versus (omega_max, a_max)."""
import sys, time, json, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import ra_mpc_cert as R, ra_mpc_dubins as rm
from eval_mpc import S, cand

K, MU, SC, DF, DT = "k1/W3_l0.06.txt", 0.02, 0.159, 0.2, 0.05

def study(um, N, plant="zoh", T_max=20.0):
    um = np.asarray(um, float); t0 = time.time()
    fv = Parallel(n_jobs=-1, prefer="processes")(delayed(R.vanilla_feasible)(S[:, i], DT, N, um) for i in cand)
    res = Parallel(n_jobs=-1, prefer="processes")(delayed(R.run_one)(S[:, i], K, MU, SC, DT, N, DF, None, T_max, plant, um) for i in cand)
    st = Counter(r["status"] for r in res); reached = [r for r in res if r["status"] == "reached"]
    U = np.concatenate([r["u"] for r in reached if len(r["u"])]) if reached else np.zeros((1, 2))
    sat = (np.abs(U) >= 0.98 * um).mean(0)
    te = np.array([r["check"]["t_entry"] for r in reached]) if reached else np.array([np.nan])
    dec = np.concatenate([np.diff(r["J"]) + rm.ELL_LOWER for r in reached if len(r["J"]) > 1]) if reached else np.array([0.0])
    both = sum(1 for f, r in zip(fv, res) if f and r["status"] == "reached")
    o = dict(u_max=um.tolist(), N=N, vanilla=int(sum(fv)), rampc=len(reached), both=both,
             inf_k0=sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) == 0),
             inf_later=sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) > 0),
             timeout=st.get("timeout", 0), unsafe=st.get("unsafe", 0), dec_viol=int((dec > 1e-6).sum()),
             sat_omega=float(sat[0]), sat_a=float(sat[1]), t_med=float(np.median(te)), t_max=float(np.max(te)),
             runs_with_sat=float(np.mean([np.any(np.abs(r["u"]) >= 0.98 * um) for r in reached])) if reached else 0.0)
    print(f"u_max=({um[0]:g},{um[1]:g}) N={N}: vanilla {o['vanilla']:3d} | RA-MPC {o['rampc']:3d} (k=0 infeasible {o['inf_k0']}, later {o['inf_later']}, timeout {o['timeout']}, unsafe {o['unsafe']}, "
          f"decrease viol {o['dec_viol']}) | saturated steps: omega {sat[0]*100:4.1f}%, a {sat[1]*100:4.1f}%, runs with any saturation {o['runs_with_sat']*100:3.0f}% | "
          f"t_entry med {o['t_med']:.2f} max {o['t_max']:.2f} s ({time.time()-t0:.0f}s)", flush=True)
    return o

if __name__ == "__main__":
    grid = [(5, 5), (2, 2), (2, 1), (1, 2), (1, 1), (1, 0.5), (0.5, 0.5)]
    Ns = [int(a) for a in (sys.argv[1].split(",") if len(sys.argv) > 1 else ["25", "40"])]
    out = []
    for N in Ns:
        for um in grid:
            out.append(study(um, N))
            json.dump(out, open("logs/eval_bounds.json", "w"), indent=1)
    open("logs/eval_bounds.done", "w").write("done")
