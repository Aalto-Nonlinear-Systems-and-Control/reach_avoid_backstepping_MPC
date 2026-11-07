"""Per-step solve times measured one run at a time (no parallel load). Writes data/timing.pkl."""
import pickle, time, numpy as np
import config as cf

M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0 = M["X0"]
ok = np.array([u["status"] == "reached" for u in M["ura"]]) & M["fv"].astype(bool) & np.array([L["status"] == "reached" for L in M["logs"]])
idx = np.where(ok)[0][:: max(1, ok.sum() // 12)][:12]
ra, va, ur = [], [], []
for i in idx:
    L = cf.ra_run(X0[:, i]); ra.append(L["solve_ms"]); V = cf.vanilla_rollout(X0[:, i]); va.append(V["solve_ms"])
    c = cf.C; x = X0[:, [i]]; t0 = time.perf_counter()
    for _ in range(2000):
        c.uc(x)
    ur.append(1e3 * (time.perf_counter() - t0) / 2000)
f = lambda a: dict(mean=float(np.mean(np.concatenate(a))), median=float(np.median(np.concatenate(a))), first=float(np.mean([b[0] for b in a])),
                   later=float(np.mean(np.concatenate([b[1:] for b in a]))), n=int(sum(len(b) for b in a)))
out = dict(idx=idx, ra=f(ra), vanilla=f(va), ura_ms=float(np.mean(ur))); print(out)
pickle.dump(out, open(f"{cf.DATA}/timing.pkl", "wb"))
