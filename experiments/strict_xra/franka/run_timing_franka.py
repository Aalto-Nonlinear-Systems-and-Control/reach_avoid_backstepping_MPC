"""Per-step solve times of RA-MPC and vanilla MPC on the Panda-like arm, one run at a time in a single process (as
run_timing_top.py for the two-link arm), excluding the one-off JIT compilation; terminal controller evaluation time.
Writes data/red_fwd_sq/timing.pkl."""
import pickle, time, numpy as np
import config_franka as cf
M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0 = M["X0"]
ok_r = np.array([L["status"] == "reached" for L in M["logs"]]); ok_v = np.array([v["status"] == "reached" for v in M["van"]])
idx_r = np.where(ok_r)[0][:: max(1, ok_r.sum() // 5)][:5]; idx_v = np.where(ok_v)[0][:: max(1, ok_v.sum() // 5)][:5]; ra_ms, va_ms, ur = [], [], []
t0 = time.time(); cf.mpc_run(X0[:, idx_r[0]], "cert"); cf.mpc_run(X0[:, idx_v[0]], "target"); print(f"warm-up incl. JIT {time.time() - t0:.0f}s", flush=True)
for i in idx_r:
    L = cf.mpc_run(X0[:, i], "cert"); ra_ms.append(np.asarray(L["solve_ms"])); x = X0[:, i]; t1 = time.perf_counter()
    for _ in range(2000): cf.C.F["uc"](x)
    ur.append(1e3 * (time.perf_counter() - t1) / 2000)
for i in idx_v:
    V = cf.mpc_run(X0[:, i], "target"); va_ms.append(np.asarray(V["solve_ms"]))
f = lambda a: dict(mean=float(np.mean(np.concatenate(a))), median=float(np.median(np.concatenate(a))), first=float(np.mean([b[0] for b in a])), later=float(np.mean(np.concatenate([b[1:] for b in a]))), n=int(sum(len(b) for b in a)))
out = dict(idx_r=idx_r, idx_v=idx_v, ra=f(ra_ms), vanilla=f(va_ms), ura_ms=float(np.mean(ur)), machine="Apple M5, 32 GB, single process"); print(out); pickle.dump(out, open(f"{cf.DATA}/timing.pkl", "wb"))
