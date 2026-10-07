"""Per-step solve times measured one run at a time."""
import pickle, time, numpy as np
import config_arm as cf
M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0 = M["X0"]
ok_r = np.array([L["status"] == "reached" for L in M["logs"]]); ok_v = np.array([v["status"] == "reached" for v in M["van"]])
idx_r = np.where(ok_r)[0][:: max(1, ok_r.sum() // 5)][:5]; idx_v = np.where(ok_v)[0][:5]; ra_ms, va_ms, ur = [], [], []
for i in idx_r:
    L = cf.mpc_run(X0[:, i], "cert"); ra_ms.append(L["solve_ms"]); x = X0[:, [i]]; t0 = time.perf_counter()
    for _ in range(2000): cf.C.uc(x)
    ur.append(1e3 * (time.perf_counter() - t0) / 2000)
for i in idx_v:
    V = cf.mpc_run(X0[:, i], "target"); va_ms.append(V["solve_ms"])
f = lambda a: dict(mean=float(np.mean(np.concatenate(a))), median=float(np.median(np.concatenate(a))), first=float(np.mean([b[0] for b in a])), later=float(np.mean(np.concatenate([b[1:] for b in a]))), n=int(sum(len(b) for b in a))) if a else None
out = dict(idx_r=idx_r, idx_v=idx_v, ra=f(ra_ms), vanilla=f(va_ms), ura_ms=float(np.mean(ur))); print(out); pickle.dump(out, open(f"{cf.DATA}/timing.pkl", "wb"))
