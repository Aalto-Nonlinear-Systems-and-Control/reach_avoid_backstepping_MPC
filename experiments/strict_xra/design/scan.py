"""Scan (weight, lambda, mu, delta): terminal-set size, singularity, input bounds, eta, M2, T*, simulated limits."""
import sys, itertools, numpy as np
from cert import Cert
TG = [0.005, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.1]

def metrics(c, dfrac, n_y=150000, seed=0, T_ref=0.05):
    pm = 1.0 - c.shift; delta = dfrac * pm
    x, info = c.sample(delta, n_y=n_y, seed=seed)
    x0, info0 = c.sample(0.0, n_y=n_y, seed=seed)
    share = info["w"].sum() / info0["w"].sum()
    vmin = float((info["kn"] - info["rad"]).min())          # smallest speed in the set (<=0: contains v = 0)
    ok = x[3] > 0.02; x = x[:, ok]
    u = c.uc(x); inb = (np.abs(u) <= 5.0).all(0)
    p0 = c.Psi(x); pd = c.Psidot(x, u)
    _, M2, pmin = c.zoh(x, u, T_ref, track=True)
    eta = float(pd.min()); M2max = float(M2.max())
    res = dict(share=share, vmin=vmin, frac_slow=float(1 - ok.mean()), inb=float(inb.mean()), umax=np.abs(u).max(1),
               eta=eta, eta1=float(np.quantile(pd, 0.01)), M2=M2max, M2_99=float(np.quantile(M2, 0.99)),
               Tb=2 * eta / M2max, Tp=float((2 * pd / np.maximum(M2, 1e-9)).min()), n=x.shape[1])
    inc, inv, eps = 0.0, 0.0, {}
    for T in TG:
        xn, _, pmn = c.zoh(x, u, T, track=True)
        d = c.Psi(xn) - p0
        eps[T] = float(d.min())
        if np.all(d > 0): inc = T
        if np.all(c.Psi(xn) >= delta) and np.all(pmn >= 0): inv = T
    res.update(T_inc=inc, T_inv=inv, eps=eps, frac_dec_ref=float((c.Psi(c.zoh(x, u, T_ref)) - p0 < 0).mean()))
    return res

if __name__ == "__main__":
    tags = sys.argv[1].split(",") if len(sys.argv) > 1 else ["W3"]
    lams = sys.argv[2].split(",") if len(sys.argv) > 2 else ["0.005", "0.02", "0.05"]
    mus = [float(a) for a in (sys.argv[3].split(",") if len(sys.argv) > 3 else ["0.05", "0.1", "0.2"])]
    dfs = [float(a) for a in (sys.argv[4].split(",") if len(sys.argv) > 4 else ["0.1", "0.2", "0.3", "0.5"])]
    for tag, l in itertools.product(tags, lams):
        c0 = Cert(f"k1/{tag}_l{l}.txt", mu=0.1); g = c0.margin_on_grid()
        print(f"\n=== {tag}  lambda={c0.lam:.4g}  delta_s/lambda={c0.shift:.4f}  min k1-margin={g['min_margin']:+.1e}  max|k1_i|={g['kmax']:.3f}  output-area share={g['area_share']*100:.0f}%")
        print("  mu   d/max | share  vmin   |u|<=5  max|w| max|a| |   eta    M2(99%)      T*=2eta/M2  min 2Pd/M2loc | T_inc  T_inv | dec@0.05  eps@0.05")
        for mu in mus:
            c = Cert(f"k1/{tag}_l{l}.txt", mu=mu)
            for df in dfs:
                r = metrics(c, df)
                print(f"  {mu:<4g} {df:<5g} | {r['share']*100:4.0f}%  {r['vmin']:+.2f}  {r['inb']*100:5.1f}%  {r['umax'][0]:5.1f}  {r['umax'][1]:5.1f} | "
                      f"{r['eta']:+.4f}  {r['M2']:6.1f}({r['M2_99']:5.1f})  {r['Tb']:.2e}    {r['Tp']:.2e}  | {r['T_inc']:<5g}  {r['T_inv']:<5g} | {r['frac_dec_ref']*100:5.2f}%   {r['eps'][0.05]:+.1e}", flush=True)
