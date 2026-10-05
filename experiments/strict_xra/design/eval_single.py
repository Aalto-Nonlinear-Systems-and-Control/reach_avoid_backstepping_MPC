"""One initial set for all three methods (R2): positions uniform in X_S minus X_T, velocity uniform in a disc around k1(y)."""
import sys, time, json, numpy as np
from collections import Counter
from joblib import Parallel, delayed
from scipy.integrate import solve_ivp
from cert import Cert
import ra_mpc_cert as R, ra_mpc_dubins as rm
from eval_init import K, MU, SC, DF, DT, C, DELTA, POS


def build(r0, seed=2026):
    """Same unit-disc draws for every r0, so the sets are nested in r0."""
    r = np.random.default_rng(seed); n = POS.shape[1]
    w = np.sqrt(r.random(n)) * np.exp(1j * r.uniform(0, 2 * np.pi, n))
    k = np.array([np.broadcast_to(np.asarray(a, float), (n,)) for a in C.k1(*POS)])
    yd = k + r0 * np.array([w.real, w.imag])
    return np.vstack([POS, np.arctan2(yd[1], yd[0]), np.hypot(yd[0], yd[1])])


def ura_run(x0, T_max=600.0):
    """u^c alone in continuous time; returns status, time, peak inputs, accumulated stage cost."""
    c = Cert(K, mu=MU).scaled(SC)
    def rhs(t, z):
        x = z[:4]; u = c.uc(x.reshape(4, 1))[:, 0]; dy = x[:2] - rm.TARGET_POS
        return [x[3] * np.cos(x[2]), x[3] * np.sin(x[2]), u[0], u[1], rm.Q_Y * dy @ dy + rm.R_U * u @ u]
    e1 = lambda t, z: float(rm.tgt_np(*z[:4])); e1.terminal = True; e1.direction = -1
    e2 = lambda t, z: float(rm.safe_np(*z[:4])); e2.terminal = True; e2.direction = -1
    e3 = lambda t, z: abs(z[3]) - 1e-4; e3.terminal = True; e3.direction = -1
    try:
        s = solve_ivp(rhs, (0, T_max), np.append(x0, 0.0), events=[e1, e2, e3], rtol=1e-6, atol=1e-8, max_step=0.1)
    except Exception:
        return "error", np.nan, np.zeros(2), np.nan
    u = np.abs(c.uc(s.y[:4])).max(1)
    st = "reached" if s.t_events[0].size else "unsafe" if s.t_events[1].size else "v=0" if s.t_events[2].size else "timeout"
    return st, s.t[-1], u, s.y[4, -1]


def ura_all(X0):
    out = Parallel(n_jobs=-1, prefer="processes")(delayed(ura_run)(X0[:, i]) for i in range(X0.shape[1]))
    return (np.array([o[0] for o in out]), np.array([o[1] for o in out]), np.array([o[2] for o in out]).T, np.array([o[3] for o in out]))


def mpc_all(X0, um, N):
    um = np.asarray(um, float); n = X0.shape[1]
    fv = np.array(Parallel(n_jobs=-1, prefer="processes")(delayed(R.vanilla_feasible)(X0[:, i], DT, N, um) for i in range(n)))
    res = Parallel(n_jobs=-1, prefer="processes")(delayed(R.run_one)(X0[:, i], K, MU, SC, DT, N, DF, None, 20.0, "zoh", um) for i in range(n))
    ra = np.array([r["status"] == "reached" for r in res]); later = sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) > 0)
    te = np.array([r["check"]["t_entry"] if r["status"] == "reached" else np.nan for r in res])
    cost = np.array([DT * float(np.sum(r["ell"])) if r["status"] == "reached" else np.nan for r in res])
    up = np.array([np.abs(r["u"]).max(0) if len(r["u"]) else np.zeros(2) for r in res]).T
    return fv, ra, te, cost, up, later, Counter(r["status"] for r in res)


if __name__ == "__main__":
    med = lambda a: float(np.nanmedian(a)) if np.isfinite(a).any() else float("nan")
    radii = [float(a) for a in sys.argv[1].split(",")] if len(sys.argv) > 1 else [0.01, 0.02, 0.03, 0.04, 0.05]
    mpc_r = [float(a) for a in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.02, 0.03]
    sets = [((3, 3), 25), ((3, 3), 30), ((2, 2), 40)]
    n = POS.shape[1]; rad = np.sqrt(2 * C.mu * np.clip(np.broadcast_to(C.psi(*POS), (n,)) - C.shift, 0, None))
    print(f"X_RA disc radius at the 300 positions: max {rad.max():.4f}, median over certified positions {np.median(rad[rad > 0]):.4f}", flush=True)
    out = []; U = {}
    for r0 in radii:
        t0 = time.time(); X0 = build(r0); P = C.Psi(X0); st, tt, uu, cc = ura_all(X0); U[r0] = (X0, P, st, tt, uu, cc); ok = st == "reached"
        print(f"r0={r0:g}: speed [{X0[3].min():.3f}, {X0[3].max():.3f}] | in X_RA {int((P>=0).sum())}, in terminal set {int((P>=DELTA).sum())} | URA {dict(Counter(st))}, "
              f"of X_RA {int((ok&(P>=0)).sum())}/{int((P>=0).sum())}, time med {med(tt[ok]):.1f} s, peak omega {uu[0][ok].max():.2f}, peak a {uu[1][ok].max():.3f}, "
              f"peak over all runs {uu[0].max():.2f}/{uu[1].max():.2f} ({time.time()-t0:.0f}s)", flush=True)
    for r0 in mpc_r:
        X0, P, st, tt, uu, cc = U[r0] if r0 in U else (lambda X: (X, C.Psi(X), *ura_all(X)))(build(r0)); ok = st == "reached"
        for um, N in sets:
            t0 = time.time(); fv, ra, te, cost, up, later, stc = mpc_all(X0, um, N); both = ok & ra
            viol = ok & ((uu[0] > um[0]) | (uu[1] > um[1]))
            o = dict(r0=r0, u_max=list(um), N=N, in_xra=int((P >= 0).sum()), in_terminal=int((P >= DELTA).sum()),
                     ura=int(ok.sum()), ura_within_bounds=int((ok & ~viol).sum()), vanilla=int(fv.sum()), rampc=int(ra.sum()),
                     ra_not_vanilla=int((ra & ~fv).sum()), vanilla_not_ra=int((fv & ~ra).sum()), ra_not_ura=int((ra & ~ok).sum()), ura_not_ra=int((ok & ~ra).sum()),
                     ra_in_xra=int((ra & (P >= 0)).sum()), ra_in_terminal=int((ra & (P >= DELTA)).sum()), ra_later_infeasible=later, ra_status=dict(stc),
                     common=int(both.sum()), t_ura_med=med(tt[both]), t_ra_med=med(te[both]), cost_ura_med=med(cc[both]), cost_ra_med=med(cost[both]),
                     peak_ura=[float(uu[0][both].max()), float(uu[1][both].max())], peak_ra=[float(up[0][ra].max()), float(up[1][ra].max())])
            out.append(o); print(json.dumps(o), f"({time.time()-t0:.0f}s)", flush=True)
            np.savez(f"logs/single_r{r0:g}_{um[0]}_{um[1]}_N{N}.npz", X0=X0, Psi=P, ura=ok, ura_status=st, t_ura=tt, u_ura=uu, cost_ura=cc, fv=fv, ra=ra, te=te, cost_ra=cost, u_ra=up)
    json.dump(out, open("logs/eval_single.json", "w"), indent=1)
