"""Lemma 2 on the final certificate: one-step increase of Psi under zero-order hold, sampled-data closed loop,
period sweep, and the sampling-period bounds. Writes data/lemma2.pkl."""
import pickle, time, numpy as np
import config as cf
from config import C, rm

TG = [0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
TG_CL = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]


def arc(x, u, T):
    """One ZOH arc of length T: end state, entered X_T, min Psi and min safe value on the arc, sup |d Psidot/dt|."""
    n = max(40, int(round(T / 0.0025))); h = T / n; x = x.copy(); m = x.shape[1]
    touched = np.zeros(m, bool); pmin = C.Psi(x); smin = rm.safe_np(*x) + 0 * pmin; pd0 = C.Psidot(x, u); M2 = np.zeros(m)
    for _ in range(n):
        a1 = cf.f_np(x, u); a2 = cf.f_np(x + h / 2 * a1, u); a3 = cf.f_np(x + h / 2 * a2, u); a4 = cf.f_np(x + h * a3, u)
        x = x + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        touched |= rm.tgt_np(*x) <= 0
        pd = C.Psidot(x, u); M2 = np.where(touched, M2, np.maximum(M2, np.abs(pd - pd0) / h)); pd0 = pd
        pmin = np.where(touched, pmin, np.minimum(pmin, C.Psi(x))); smin = np.where(touched, smin, np.minimum(smin, rm.safe_np(*x)))
    return x, touched, pmin, smin, M2


def paper_constants(x, h=1e-6):
    """Constants of the closed-form bound in the proof of Lemma 2, by central differences on the samples."""
    e = np.eye(4)
    dPsi = np.array([(C.Psi(x + h * e[:, [i]]) - C.Psi(x - h * e[:, [i]])) / (2 * h) for i in range(4)])
    Mg = np.abs(dPsi[2:]).max(1)                                   # L_g1 Psi = dPsi/dtheta, L_g2 Psi = dPsi/dv
    du = np.array([(C.uc(x + h * e[:, [i]]) - C.uc(x - h * e[:, [i]])) / (2 * h) for i in range(4)])   # (4, 2, n)
    lu = np.sqrt((du**2).sum(0)).max(1)
    u = C.uc(x); Mv = np.sqrt(x[3]**2 + (u**2).sum(0)).max()
    return Mg, lu, Mv


if __name__ == "__main__":
    t0 = time.time()
    xd, _ = C.sample(cf.DELTA, n_y=60000, seed=3); x0, _ = C.sample(0.0, n_y=60000, seed=4)
    ud = C.uc(xd); pd = C.Psidot(xd, ud); P = C.Psi(xd)
    eta = float(pd.min()); out = dict(TG=TG, TG_CL=TG_CL, n_delta=xd.shape[1], n_0=x0.shape[1], eta=eta, lam_delta=C.lam * cf.DELTA,
                                      umax_delta=np.abs(ud).max(1), umax_0=np.abs(C.uc(x0)).max(1), speed_delta=[float(xd[3].min()), float(xd[3].max())])
    rows = []
    for T in TG:
        xn, touched, pmin, smin, M2 = arc(xd, ud, T); k = ~touched; d = (C.Psi(xn) - P)[k]
        rows.append(dict(T=T, n=int(k.sum()), dmin=float(d.min()), dq01=float(np.quantile(d, 0.01)), dmed=float(np.median(d)), frac_inc=float((d > 0).mean()),
                         frac_inv=float((C.Psi(xn)[k] >= cf.DELTA).mean()), frac_in_xra=float((pmin[k] >= 0).mean()), safe_min=float(smin[k].min()),
                         M2=float(M2[k].max()), lb2=float(eta * T - 0.5 * M2[k].max() * T * T)))
        print(f"T={T:<6g} arcs {k.sum():5d} | min dPsi {d.min():+.3e} (1% {np.quantile(d,0.01):+.2e}) | increasing {100*(d>0).mean():6.2f}% | "
              f"end in terminal set {100*rows[-1]['frac_inv']:6.2f}% | arc inside X_RA {100*rows[-1]['frac_in_xra']:6.2f}% | M2 {rows[-1]['M2']:.3f} | 2nd-order bound {rows[-1]['lb2']:+.2e}", flush=True)
    out["onestep"] = rows
    # eta and M2 on a larger sample (sup and inf estimates grow with the sample), used for T* = 2 eta / M2
    xb, _ = C.sample(cf.DELTA, n_y=500000, seed=11); ub = C.uc(xb); eta_b = float(C.Psidot(xb, ub).min())
    _, tb, _, _, M2b = arc(xb, ub, 0.05); out["n_big"] = xb.shape[1]
    eta = min(eta, eta_b); out["eta"] = eta
    M2ref = max([r for r in rows if r["T"] == 0.05][0]["M2"], float(M2b[~tb].max())); out["M2_ref"] = M2ref; out["Tstar2"] = 2 * eta / M2ref
    d_b = (C.Psi(cf.zoh_np(xb, ub, 0.05)) - C.Psi(xb))[~tb]; out["eps_big"] = float(d_b.min())
    ell_b = rm.Q_Y * ((xb[0] - rm.TARGET_POS[0])**2 + (xb[1] - rm.TARGET_POS[1])**2) + rm.R_U * (ub**2).sum(0)
    out["vf_check"] = dict(dpsi=d_b[::8], ell=ell_b[~tb][::8], sup_ell=float(ell_b.max()), n=int((~tb).sum()))
    print(f"large sample: {xb.shape[1]} states, eta {eta_b:.4e}, M2 {M2b[~tb].max():.4f}, min one-step increase at T=0.05 {d_b.min():.3e}, sup stage cost {ell_b.max():.3f}", flush=True)
    # self-consistent T*: largest T with eta*T - M2(T) T^2/2 > 0, M2 taken over arcs of that length
    Tf = np.geomspace(0.005, 1.0, 60); ok = []
    for T in Tf:
        _, touched, _, _, M2 = arc(xd[:, ::4], ud[:, ::4], T); ok.append(eta - 0.5 * M2[~touched].max() * T > 0)
    out["Tstar2_self"] = float(Tf[np.where(ok)[0].max()]) if any(ok) else 0.0
    Mg, lu, Mv = paper_constants(xd); a = float((Mg * lu).sum() * Mv)
    out.update(Mg=Mg, lu=lu, Mv=float(Mv), Tstar_paper=float(np.log(1 + C.lam * cf.DELTA / a)))
    print(f"eta = {eta:.4e} (lambda*delta = {C.lam*cf.DELTA:.4e}), M2 at T=0.05: {M2ref:.4f} -> T* = 2 eta/M2 = {out['Tstar2']*1e3:.1f} ms (self-consistent {out['Tstar2_self']*1e3:.1f} ms)", flush=True)
    print(f"paper bound: M_g = {np.round(Mg,3)}, l_uc = {np.round(lu,2)}, M_v = {Mv:.3f} -> T* = {out['Tstar_paper']:.3e} s", flush=True)
    # sampled-data closed loop
    r = np.random.default_rng(0); Xd = xd[:, r.choice(xd.shape[1], 600, replace=False)]; Xo = x0[:, r.choice(x0.shape[1], 600, replace=False)]
    cl = []
    for T in TG_CL:
        row = dict(T=T)
        for nm, X in (("delta", Xd), ("zero", Xo)):
            st, tt, dm = cf.ura_zoh(X, dt=T, T_max=400.0)
            row[nm] = dict(reached=int((st == 1).sum()), unsafe=int((st == 2).sum()), timeout=int((st == 0).sum()), n=X.shape[1],
                           monotone=int(((dm > 0) | ~np.isfinite(dm)).sum()), dmin=float(dm[np.isfinite(dm)].min()), t_med=float(np.nanmedian(tt)), t_max=float(np.nanmax(tt)) if np.isfinite(tt).any() else np.nan)
        cl.append(row)
        print(f"T={T:<5g} from X_RA^delta: reached {row['delta']['reached']}/600, Psi increasing at all samples {row['delta']['monotone']}/600 (min step {row['delta']['dmin']:+.2e}) | "
              f"from X_RA: reached {row['zero']['reached']}/600, increasing {row['zero']['monotone']}/600 ({time.time()-t0:.0f}s)", flush=True)
    out["closed_loop"] = cl
    st, tt, dm, hist = cf.ura_zoh(Xd[:, :60], dt=cf.T, T_max=400.0, track=True)
    out["traj"] = dict(psi=hist, status=st, t_reach=tt, dt=cf.T)
    pickle.dump(out, open(f"{cf.DATA}/lemma2.pkl", "wb")); print(f"saved data/lemma2.pkl ({time.time()-t0:.0f}s)")
