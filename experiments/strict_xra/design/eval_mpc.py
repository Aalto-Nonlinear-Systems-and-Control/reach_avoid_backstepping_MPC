"""Closed-loop evaluation of Algorithm 1 for candidate certificate settings (300 initial states)."""
import sys, time, json, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import ra_mpc_cert as R, ra_mpc_dubins as rm

np.random.seed(42); lo = np.array([-2, -2, 2 * np.pi / 3, -1.0]); hi = np.array([2, 2, 4 * np.pi / 3, 1.0])
S = np.random.rand(4, 1000) * (hi - lo)[:, None] + lo[:, None]
cand = np.where((rm.safe_np(*S) >= 0) & (rm.tgt_np(*S) > 0))[0][:300]

def evaluate(name, k1_file, mu, scale, dt, N, dfrac, omega=None, T_max=10.0, plant="nominal", verbose=True):
    t0 = time.time()
    res = Parallel(n_jobs=-1, prefer="processes")(
        delayed(R.run_one)(S[:, i], k1_file, mu, scale, dt, N, dfrac, omega, T_max, plant) for i in cand)
    st = Counter(r["status"] for r in res); reached = [r for r in res if r["status"] == "reached"]
    inf0 = sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) == 0)
    infl = sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) > 0)
    dec = np.concatenate([np.diff(r["J"]) + rm.ELL_LOWER for r in reached if len(r["J"]) > 1]) if reached else np.array([0.0])
    te = np.array([r["check"]["t_entry"] for r in reached]) if reached else np.array([np.nan])
    tb = np.array([r["check"]["t_bound"] for r in reached]) if reached else np.array([np.nan])
    ms = np.mean(np.concatenate([r["solve_ms"] for r in res if len(r["solve_ms"])])) if any(len(r["solve_ms"]) for r in res) else np.nan
    out = dict(name=name, plant=plant, reached=len(reached), infeasible_k0=inf0, infeasible_later=infl,
               timeout=st.get("timeout", 0), unsafe=st.get("unsafe", 0), other=sum(v for k, v in st.items() if k not in ("reached", "infeasible", "timeout", "unsafe")),
               dec_viol_steps=int((dec > 1e-6).sum()), steps=int(dec.size), bound_ok=int(sum(r["check"]["bound_ok"] for r in reached)),
               t_entry_med=float(np.median(te)), t_entry_max=float(np.max(te)), t_bound_med=float(np.median(tb)),
               omega=float(res[0]["omega"]), eps=float(res[0]["eps"]), solve_ms=float(ms), wall=time.time() - t0)
    if verbose:
        print(f"{name:34s} [{plant:7s}] reached {out['reached']:3d}/300 | infeasible k=0 {inf0:3d}, later {infl:2d} | timeout {out['timeout']:2d} unsafe {out['unsafe']:2d} | "
              f"decrease viol {out['dec_viol_steps']}/{out['steps']} | bound ok {out['bound_ok']}/{out['reached']} | t_entry med {out['t_entry_med']:.2f} max {out['t_entry_max']:.2f} s, "
              f"t_bound med {out['t_bound_med']:.3g} s | eps {out['eps']:.1e} omega {out['omega']:.2g} | {ms:.0f} ms/solve ({out['wall']:.0f}s)", flush=True)
    return out, res

if __name__ == "__main__":
    K = "k1/W3_l0.06.txt"
    Tstar = {(0.02, 0.2): 8.39e-3, (0.02, 0.3): 1.40e-2, (0.02, 0.4): 2.18e-2,
             (0.05, 0.2): 5.76e-3, (0.05, 0.3): 9.89e-3, (0.05, 0.4): 1.76e-2}   # uniform second-order bound at s = 1
    Ns = [int(a) for a in (sys.argv[1].split(",") if len(sys.argv) > 1 else ["25", "40"])]
    allout = []
    for N in Ns:
        t0 = time.time()
        fv = Parallel(n_jobs=-1, prefer="processes")(delayed(R.vanilla_feasible)(S[:, i], 0.05, N) for i in cand)
        print(f"vanilla MPC (terminal target), N={N}: first-step feasible {sum(fv)}/300 ({time.time()-t0:.0f}s)", flush=True)
        for (mu, df), ts in Tstar.items():
            s = round(0.95 * ts / 0.05, 3)
            o, _ = evaluate(f"mu={mu} d={df} s={s} N={N}", K, mu, s, 0.05, N, df); o.update(N=N, mu=mu, dfrac=df, s=s); allout.append(o)
    json.dump(allout, open("logs/eval_mpc_certified.json", "w"), indent=1)
