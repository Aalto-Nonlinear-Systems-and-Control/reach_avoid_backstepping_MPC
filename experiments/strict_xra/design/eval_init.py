"""Initial-state study for setting R2: at rest, near rest, and field-relative sampling."""
import sys, time, json, numpy as np
from collections import Counter
from joblib import Parallel, delayed
from scipy.integrate import solve_ivp
from cert import Cert
import ra_mpc_cert as R, ra_mpc_dubins as rm
from eval_mpc import S, cand

K, MU, SC, DF, DT = "k1/W3_l0.06.txt", 0.02, 0.159, 0.2, 0.05
C = Cert(K, mu=MU).scaled(SC)
DELTA = DF * (1 - C.shift)
POS = S[:2, cand]                       # the 300 positions used so far


def ura_one(x0, T_max, um):
    c = Cert(K, mu=MU).scaled(SC)
    if abs(x0[3]) < 1e-6:
        return "singular", np.nan, np.zeros(2)
    def rhs(t, x):
        u = c.uc(x.reshape(4, 1))[:, 0]
        return [x[3] * np.cos(x[2]), x[3] * np.sin(x[2]), u[0], u[1]]
    e1 = lambda t, x: float(rm.tgt_np(*x)); e1.terminal = True; e1.direction = -1
    e2 = lambda t, x: float(rm.safe_np(*x)); e2.terminal = True; e2.direction = -1
    e3 = lambda t, x: abs(x[3]) - 1e-4; e3.terminal = True; e3.direction = -1
    try:
        s = solve_ivp(rhs, (0, T_max), x0, events=[e1, e2, e3], rtol=1e-6, atol=1e-8, max_step=0.1)
    except Exception:
        return "error", np.nan, np.zeros(2)
    u = np.abs(c.uc(s.y)).max(1)
    st = "reached" if s.t_events[0].size else "unsafe" if s.t_events[1].size else "v=0" if s.t_events[2].size else "timeout"
    return st, s.t[-1], u


def study(name, X0, um, N=25, with_ura=True, T_ura=400.0):
    um = np.asarray(um, float); n = X0.shape[1]; t0 = time.time()
    P = C.Psi(np.where(np.abs(X0[3:4]) < 1e-9, np.vstack([X0[:3], 1e-9 * np.ones((1, n))]), X0))
    in0 = int((P >= 0).sum()); inT = int((P >= DELTA).sum())
    fv = Parallel(n_jobs=-1, prefer="processes")(delayed(R.vanilla_feasible)(X0[:, i], DT, N, um) for i in range(n))
    res = Parallel(n_jobs=-1, prefer="processes")(delayed(R.run_one)(X0[:, i], K, MU, SC, DT, N, DF, None, 20.0, "zoh", um) for i in range(n))
    ra = np.array([r["status"] == "reached" for r in res]); bad = sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) > 0)
    te = np.array([r["check"]["t_entry"] for r in res if r["status"] == "reached"])
    line = (f"{name:44s} u_max=({um[0]:g},{um[1]:g}) N={N}: in X_RA {in0:3d}, in terminal set {inT:3d} | vanilla {sum(fv):3d} | RA-MPC {ra.sum():3d} "
            f"(later infeasible {bad}, t_entry med {np.median(te) if te.size else float('nan'):.2f} s)")
    o = dict(name=name, u_max=um.tolist(), N=N, n=n, in_xra=in0, in_terminal=inT, vanilla=int(sum(fv)), rampc=int(ra.sum()), inf_later=bad)
    if with_ura:
        ur = Parallel(n_jobs=-1, prefer="processes")(delayed(ura_one)(X0[:, i], T_ura, um) for i in range(n))
        st = Counter(u[0] for u in ur); ok = np.array([u[0] == "reached" for u in ur])
        viol = np.array([np.any(u[2] > um) for u in ur])
        tt = np.array([u[1] for u in ur if u[0] == "reached"])
        inside = P >= 0
        line += (f" | controller alone: reached {ok.sum():3d} (of those inside X_RA: {int((ok & inside).sum())}/{in0}), {dict(st)}, "
                 f"runs exceeding the bounds {int(viol.sum())}, reach time med {np.median(tt) if tt.size else float('nan'):.1f} s")
        o.update(ura=int(ok.sum()), ura_inside=int((ok & inside).sum()), ura_viol=int(viol.sum()), ura_status=dict(st))
    print(line + f"  ({time.time()-t0:.0f}s)", flush=True)
    return o


if __name__ == "__main__":
    r = np.random.default_rng(7); n = POS.shape[1]; out = []
    th_old = S[2, cand]; th_any = r.uniform(-np.pi, np.pi, n)
    k = np.array([np.broadcast_to(np.asarray(a, float), (n,)) for a in C.k1(*POS)]); kn = np.hypot(k[0], k[1]); kth = np.arctan2(k[1], k[0])
    print(f"field speed |k1| at the 300 positions: min {kn.min():.3f}, median {np.median(kn):.3f}, max {kn.max():.3f}; "
          f"positions with psi >= delta_s/lambda: {int((C.psi(*POS) >= C.shift).sum())}/300", flush=True)
    which = sys.argv[1].split(",") if len(sys.argv) > 1 else ["rest", "near", "field"]
    if "rest" in which:
        for um in ((3, 3), (5, 5)):
            out.append(study("at rest, heading in [2pi/3,4pi/3]", np.vstack([POS, th_old, np.zeros(n)]), um, with_ura=False))
            out.append(study("at rest, any heading", np.vstack([POS, th_any, np.zeros(n)]), um, with_ura=False))
    if "near" in which:
        for v0 in (0.02, 0.05):
            out.append(study(f"near rest v0={v0}, any heading", np.vstack([POS, th_any, v0 * np.ones(n)]), (3, 3)))
    if "field" in which:
        for r0 in (0.02, 0.05, 0.1):
            w = r0 * np.sqrt(r.random(n)); a = r.uniform(0, 2 * np.pi, n)
            yd = k + w * np.array([np.cos(a), np.sin(a)])
            out.append(study(f"field-relative, |ydot - k1| <= {r0}", np.vstack([POS, np.arctan2(yd[1], yd[0]), np.hypot(yd[0], yd[1])]), (3, 3)))
    json.dump(out, open("logs/eval_init.json", "w"), indent=1); open("logs/eval_init.done", "w").write("done")
