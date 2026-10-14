"""Lemma 2 for the Panda terminal controller (held torque): one-step increase of Psi on the terminal set against the
period, second-order bound, and the sampled-data closed loop from X_RA^delta and X_RA (MuJoCo). Writes data/<TASK>/lemma2.pkl."""
import pickle, time, numpy as np
import config_franka as cf
from config_franka import C
from cert_franka import sample, zoh, phi_xy

TG = sorted({0.002, 0.005, 0.0075, cf.T, 0.0125, 0.015, 0.02, 0.03, 0.05}); TG_CL = sorted({0.005, cf.T, 0.0125, 0.015, 0.02, 0.03})

if __name__ == "__main__":
    t0 = time.time(); xd, _ = sample(C, cf.DELTA, n=8000, seed=3); x0, _ = sample(C, 0.0, n=8000, seed=4)
    ud = C.uc(xd); P = C.Psi(xd); eta = float(C.Psidot(xd, ud).min()); rows = []
    for T in TG:
        n = max(2, int(round(T / 0.0005))); h = T / n; X = xd.copy(); touched = np.zeros(X.shape[1], bool); pd0 = C.Psidot(X, ud); M2 = np.zeros(X.shape[1]); pmin = P.copy()
        for _ in range(n):
            X = zoh(X, ud, h, n=1); touched |= phi_xy(C.y(X)) <= 0; pd = C.Psidot(X, ud); M2 = np.where(touched, M2, np.maximum(M2, np.abs(pd - pd0) / h)); pd0 = pd; pmin = np.where(touched, pmin, np.minimum(pmin, C.Psi(X)))
        k = ~touched; d = (C.Psi(X) - P)[k]
        rows.append(dict(T=T, n=int(k.sum()), dmin=float(d.min()), frac_inc=float((d > 0).mean()), frac_inv=float((C.Psi(X)[k] >= cf.DELTA).mean()), frac_in_xra=float((pmin[k] >= 0).mean()), M2=float(M2[k].max()), lb2=float(eta * T - 0.5 * M2[k].max() * T * T)))
        print(f"T={1e3*T:<5g} ms arcs {k.sum():5d} | min dPsi {d.min():+.3e} | increasing {100*(d>0).mean():6.2f}% | end in terminal set {100*rows[-1]['frac_inv']:6.2f}% | M2 {rows[-1]['M2']:.3f} | 2nd-order bound {rows[-1]['lb2']:+.2e} ({time.time()-t0:.0f}s)", flush=True)
    M2ref = [r for r in rows if r["T"] == cf.T][0]["M2"]; out = dict(TG=TG, TG_CL=TG_CL, n_delta=xd.shape[1], eta=eta, lam_delta=C.lam * cf.DELTA, onestep=rows, M2_ref=M2ref, Tstar2=2 * eta / M2ref, umax_delta=np.abs(ud).max(1), umax_0=np.abs(C.uc(x0)).max(1))
    print(f"eta {eta:.4e} (lambda delta {C.lam*cf.DELTA:.4e}), M2 at T: {M2ref:.3f} -> second-order T* = {2e3*eta/M2ref:.2f} ms", flush=True)
    r = np.random.default_rng(0); Xd = xd[:, r.choice(xd.shape[1], 150, replace=False)]; Xo = x0[:, r.choice(x0.shape[1], 150, replace=False)]; cl = []
    for T in TG_CL:
        row = dict(T=T)
        for nm, X in (("delta", Xd), ("zero", Xo)):
            st, tt, dm = cf.uc_zoh(X, dt=T, T_max=8.0)
            row[nm] = dict(reached=int((st == 1).sum()), unsafe=int((st == 2).sum()), timeout=int((st == 0).sum()), n=X.shape[1], monotone=int(((dm > 0) | ~np.isfinite(dm)).sum()), t_med=float(np.nanmedian(tt)))
        cl.append(row); print(f"T={1e3*T:<5g} ms from X_RA^delta: reached {row['delta']['reached']}/150, increasing {row['delta']['monotone']}/150 | from X_RA: reached {row['zero']['reached']}/150, increasing {row['zero']['monotone']}/150 ({time.time()-t0:.0f}s)", flush=True)
    out["closed_loop"] = cl; pickle.dump(out, open(f"{cf.DATA}/lemma2.pkl", "wb")); print(f"saved lemma2.pkl ({time.time()-t0:.0f}s)")
