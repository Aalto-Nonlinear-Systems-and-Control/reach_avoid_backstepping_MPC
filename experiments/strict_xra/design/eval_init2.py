"""Two-group initial set for R2: (I) inside X_RA, velocity drawn from the certified disc; (II) at rest, any heading."""
import sys, time, json, numpy as np
from collections import Counter
from joblib import Parallel, delayed
from cert import Cert
import ra_mpc_cert as R, ra_mpc_dubins as rm
from eval_mpc import S, cand
from eval_init import ura_one, K, MU, SC, DF, DT, C, DELTA, POS


def zoh_uc(X0, T, T_max, nsub=4):
    """Sampled-data u^c alone (zero-order hold, period T), vectorised; status per state."""
    X = X0.copy(); n = X.shape[1]; st = np.zeros(n, int); tt = np.full(n, np.nan); h = T / nsub; um = np.zeros((2, n))
    f = lambda x, u: np.array([x[3] * np.cos(x[2]), x[3] * np.sin(x[2]), u[0], u[1]])
    for k in range(int(T_max / T)):
        act = st == 0
        if not act.any(): break
        u = np.zeros((2, n)); u[:, act] = C.uc(X[:, act]); um = np.maximum(um, np.abs(u) * act)
        for _ in range(nsub):
            k1 = f(X, u); k2 = f(X + h / 2 * k1, u); k3 = f(X + h / 2 * k2, u); k4 = f(X + h * k3, u)
            Xn = X + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); X = np.where(act, Xn, X)
            bad = act & (rm.safe_np(*X) < 0) & (st == 0); st[bad] = 2
            hit = act & (rm.tgt_np(*X) <= 0) & (st == 0); st[hit] = 1; tt[hit] = (k + 1) * T
            act = st == 0
    return st, tt, um


def three_way(name, X0, um, N, T_ura=600.0):
    um = np.asarray(um, float); n = X0.shape[1]; t0 = time.time()
    Xe = X0.copy(); Xe[3, np.abs(Xe[3]) < 1e-9] = 1e-9
    P = C.Psi(Xe); in0 = P >= 0; inT = P >= DELTA
    ur = Parallel(n_jobs=-1, prefer="processes")(delayed(ura_one)(X0[:, i], T_ura, um) for i in range(n))
    uc_ok = np.array([u[0] == "reached" for u in ur]); uc_t = np.array([u[1] if u[0] == "reached" else np.nan for u in ur])
    uc_u = np.array([u[2] for u in ur])
    rest = np.abs(X0[3]) < 1e-6
    zs, zt, zu = (np.zeros(n, int), np.full(n, np.nan), np.zeros((2, n))) if rest.all() else zoh_uc(Xe, DT, T_ura)
    fv = np.array(Parallel(n_jobs=-1, prefer="processes")(delayed(R.vanilla_feasible)(X0[:, i], DT, N, um) for i in range(n)))
    res = Parallel(n_jobs=-1, prefer="processes")(delayed(R.run_one)(X0[:, i], K, MU, SC, DT, N, DF, None, 20.0, "zoh", um) for i in range(n))
    ra = np.array([r["status"] == "reached" for r in res]); later = sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) > 0)
    te = np.array([r["check"]["t_entry"] if r["status"] == "reached" else np.nan for r in res])
    med = lambda a: float(np.nanmedian(a)) if np.isfinite(a).any() else float("nan")
    o = dict(name=name, u_max=um.tolist(), N=N, n=n, in_xra=int(in0.sum()), in_terminal=int(inT.sum()),
             uc=int(uc_ok.sum()), uc_in_xra=int((uc_ok & in0).sum()), uc_status=dict(Counter(u[0] for u in ur)), uc_t_med=med(uc_t),
             uc_umax=[float(np.nanmax(uc_u[uc_ok, 0])) if uc_ok.any() else 0.0, float(np.nanmax(uc_u[uc_ok, 1])) if uc_ok.any() else 0.0],
             zoh=int((zs == 1).sum()), zoh_in_terminal=int(((zs == 1) & inT).sum()), zoh_unsafe=int((zs == 2).sum()), zoh_t_med=med(zt),
             vanilla=int(fv.sum()), vanilla_in_xra=int((fv & in0).sum()),
             rampc=int(ra.sum()), rampc_in_xra=int((ra & in0).sum()), rampc_in_terminal=int((ra & inT).sum()), rampc_later=later, rampc_t_med=med(te),
             rampc_only=int((ra & ~fv).sum()), vanilla_only=int((fv & ~ra).sum()),
             rampc_psi_ok=int((ra & (C.psi(*X0[:2]) >= C.shift)).sum()), pos_psi_ok=int((C.psi(*X0[:2]) >= C.shift).sum()))
    print(json.dumps(o), f"({time.time()-t0:.0f}s)", flush=True)
    np.savez(f"logs/init2_{name}_{um[0]:g}_{um[1]:g}_N{N}.npz", X0=X0, Psi=P, uc_ok=uc_ok, uc_t=uc_t, zs=zs, zt=zt, fv=fv, ra=ra, te=te)
    return o


def groups(seed=11):
    """S_RA: the positions that admit a certified velocity, velocity uniform in the disc of X_RA; S_0: all positions at rest."""
    r = np.random.default_rng(seed); n = POS.shape[1]
    psi = np.broadcast_to(C.psi(*POS), (n,)); ok = psi > C.shift
    k = np.array([np.broadcast_to(np.asarray(a, float), (n,)) for a in C.k1(*POS)])
    rad = np.sqrt(2 * C.mu * np.clip(psi - C.shift, 0, None)); a = r.uniform(0, 2 * np.pi, n)
    yd = k + rad * np.sqrt(r.random(n)) * np.array([np.cos(a), np.sin(a)])
    XI = np.vstack([POS, np.arctan2(yd[1], yd[0]), np.hypot(yd[0], yd[1])])[:, ok]
    XII = np.vstack([POS, r.uniform(-np.pi, np.pi, n), np.zeros(n)])
    return XI, XII


def mpc_only(name, X0, um, N):
    um = np.asarray(um, float); n = X0.shape[1]; t0 = time.time()
    Xe = X0.copy(); Xe[3, np.abs(Xe[3]) < 1e-9] = 1e-9
    P = C.Psi(Xe); in0 = P >= 0; inT = P >= DELTA; pT = (C.psi(*X0[:2]) - C.shift) >= DELTA
    fv = np.array(Parallel(n_jobs=-1, prefer="processes")(delayed(R.vanilla_feasible)(X0[:, i], DT, N, um) for i in range(n)))
    res = Parallel(n_jobs=-1, prefer="processes")(delayed(R.run_one)(X0[:, i], K, MU, SC, DT, N, DF, None, 20.0, "zoh", um) for i in range(n))
    ra = np.array([r["status"] == "reached" for r in res]); later = sum(1 for r in res if r["status"] == "infeasible" and len(r["J"]) > 0)
    te = np.array([r["check"]["t_entry"] if r["status"] == "reached" else np.nan for r in res])
    o = dict(name=name, u_max=um.tolist(), N=N, n=n, in_xra=int(in0.sum()), in_terminal=int(inT.sum()), vanilla=int(fv.sum()), rampc=int(ra.sum()),
             rampc_in_terminal=int((ra & inT).sum()), rampc_pos_terminal=f"{int((ra & pT).sum())}/{int(pT.sum())}", rampc_later=later,
             rampc_only=int((ra & ~fv).sum()), vanilla_only=int((fv & ~ra).sum()), rampc_t_med=float(np.nanmedian(te)), rampc_t_max=float(np.nanmax(te)))
    print(json.dumps(o), f"({time.time()-t0:.0f}s)", flush=True)
    np.savez(f"logs/init3_{name}_{um[0]:g}_{um[1]:g}_N{N}.npz", X0=X0, Psi=P, fv=fv, ra=ra, te=te)
    return o


if __name__ == "__main__":
    XI, XII = groups()
    print(f"group I: {XI.shape[1]} states, speed range [{XI[3].min():.4f}, {XI[3].max():.4f}]; group II: {XII.shape[1]} states at rest", flush=True)
    out = []
    if len(sys.argv) > 1 and sys.argv[1] == "sweep":
        for um, N in (((3, 3), 30), ((3, 3), 40), ((2, 2), 25), ((2, 2), 30)):
            out.append(mpc_only("I", XI, um, N)); out.append(mpc_only("II", XII, um, N))
        json.dump(out, open("logs/eval_init3.json", "w"), indent=1)
    else:
        for um, N in (((3, 3), 25), ((2, 2), 40)):
            out.append(three_way("I", XI, um, N)); out.append(three_way("II", XII, um, N))
        json.dump(out, open("logs/eval_init2.json", "w"), indent=1)
