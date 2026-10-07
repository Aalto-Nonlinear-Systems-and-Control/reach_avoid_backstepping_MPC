"""Lemma 2 for the arm's terminal controller (held torque): one-step increase against the period, sampled-data closed
loop from X_RA^delta and X_RA, bounds. Writes data/lemma2.pkl."""
import pickle, time, numpy as np
import config_arm as cf
from config_arm import C, fk, phi_np, gravity

TG = sorted({0.0005, 0.001, cf.T, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.1}); TG_CL = sorted({0.001, cf.T, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05})


def arc(x, u, T):
    n = max(10, int(round(T / 0.0002))); h = T / n; x = x.copy(); m = x.shape[1]
    touched = np.zeros(m, bool); pmin = C.Psi(x); pd0 = C.Psidot(x, u); M2 = np.zeros(m)
    for _ in range(n):
        a1 = cf.f_np(x, u); a2 = cf.f_np(x + h / 2 * a1, u); a3 = cf.f_np(x + h / 2 * a2, u); a4 = cf.f_np(x + h * a3, u); x = x + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        touched |= phi_np(*fk(x[:2])) <= 0; pd = C.Psidot(x, u)
        M2 = np.where(touched, M2, np.maximum(M2, np.abs(pd - pd0) / h)); pd0 = pd; pmin = np.where(touched, pmin, np.minimum(pmin, C.Psi(x)))
    return x, touched, pmin, M2


if __name__ == "__main__":
    t0 = time.time(); xd, _ = C.sample(cf.DELTA, n_y=80000, seed=3); x0, _ = C.sample(0.0, n_y=80000, seed=4)
    keep = phi_np(*fk(xd[:2])) > 0; xd = xd[:, keep]; keep0 = phi_np(*fk(x0[:2])) > 0; x0 = x0[:, keep0]
    ud = C.uc(xd); P = C.Psi(xd); eta = float(C.Psidot(xd, ud).min()); rows = []
    for T in TG:
        xn, touched, pmin, M2 = arc(xd, ud, T); k = ~touched; d = (C.Psi(xn) - P)[k]
        rows.append(dict(T=T, n=int(k.sum()), dmin=float(d.min()), dq01=float(np.quantile(d, 0.01)), frac_inc=float((d > 0).mean()), frac_inv=float((C.Psi(xn)[k] >= cf.DELTA).mean()),
                         frac_in_xra=float((pmin[k] >= 0).mean()), M2=float(M2[k].max()), lb2=float(eta * T - 0.5 * M2[k].max() * T * T)))
        print(f"T={1e3*T:<5g} ms arcs {k.sum():5d} | min dPsi {d.min():+.3e} | increasing {100*(d>0).mean():6.2f}% | end in terminal set {100*rows[-1]['frac_inv']:6.2f}% | arc in X_RA {100*rows[-1]['frac_in_xra']:6.2f}% | M2 {rows[-1]['M2']:.2f} | 2nd-order bound {rows[-1]['lb2']:+.2e}", flush=True)
    M2ref = [r for r in rows if r["T"] == cf.T][0]["M2"]
    out = dict(TG=TG, TG_CL=TG_CL, n_delta=xd.shape[1], n_0=x0.shape[1], eta=eta, lam_delta=C.lam * cf.DELTA, onestep=rows, M2_ref=M2ref, Tstar2=2 * eta / M2ref,
               umax_delta=np.abs(ud).max(1), umax_0=np.abs(C.uc(x0)).max(1))
    print(f"eta {eta:.4e} (lambda delta {C.lam*cf.DELTA:.4e}), M2 at T={cf.T}: {M2ref:.3f} -> second-order T* = {2e3*eta/M2ref:.2f} ms", flush=True)
    r = np.random.default_rng(0); Xd = xd[:, r.choice(xd.shape[1], 400, replace=False)]; Xo = x0[:, r.choice(x0.shape[1], 400, replace=False)]; cl = []
    for T in TG_CL:
        row = dict(T=T)
        for nm, X in (("delta", Xd), ("zero", Xo)):
            st, tt, dm = cf.uc_zoh(X, dt=T, T_max=40.0)
            row[nm] = dict(reached=int((st == 1).sum()), unsafe=int((st == 2).sum()), timeout=int((st == 0).sum()), n=X.shape[1], monotone=int(((dm > 0) | ~np.isfinite(dm)).sum()),
                           dmin=float(dm[np.isfinite(dm)].min()) if np.isfinite(dm).any() else np.nan, t_med=float(np.nanmedian(tt)))
        cl.append(row)
        print(f"T={1e3*T:<5g} ms from X_RA^delta: reached {row['delta']['reached']}/400, Psi increasing at all samples {row['delta']['monotone']}/400 | from X_RA: reached {row['zero']['reached']}/400, increasing {row['zero']['monotone']}/400 ({time.time()-t0:.0f}s)", flush=True)
    out["closed_loop"] = cl
    st, tt, dm = cf.uc_zoh(Xd[:, :60], dt=cf.T, T_max=40.0)
    # trajectories of Psi at the samples for the figure
    X = Xd[:, :40].copy(); hist = [C.Psi(X)]; nsub = 4; h = cf.T / nsub; act = np.ones(40, bool)
    for k in range(int(40.0 / cf.T)):
        u = C.uc(X)
        for _ in range(nsub):
            a1 = cf.f_np(X, u); a2 = cf.f_np(X + h / 2 * a1, u); a3 = cf.f_np(X + h / 2 * a2, u); a4 = cf.f_np(X + h * a3, u); X = np.where(act, X + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4), X)
        act &= phi_np(*fk(X[:2])) > 0; hist.append(np.where(act, C.Psi(X), np.nan))
        if not act.any(): break
    out["traj"] = dict(psi=np.array(hist), dt=cf.T)
    pickle.dump(out, open(f"{cf.DATA}/lemma2.pkl", "wb")); print(f"saved data/lemma2.pkl ({time.time()-t0:.0f}s)")
