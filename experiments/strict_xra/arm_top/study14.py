"""Study: k1 from program (14) solved directly (input bounds enforced at sampled states) against (12) plus time scaling.
At a time scale s the target controller has lambda = s lambda0 and mu = s^2 mu0, as the scaled design A = (12) at lambda0
scaled by s. Design B solves (14) at the same lambda and mu: the ECGBF condition of (12) as an SOS constraint
(solve_k1_sampled14.m), |k1_i| <= kfac s on {psi >= 0} (kfac = 0: no bound), and |u^c_i(x_j)| <= U_i at states x_j of
the terminal set X_RA^delta. u^c(x_j) is affine in the coefficients of k1 (Theorem 1), so these are linear constraints.
The terminal set depends on k1, hence an exchange loop: solve, sample the new design's terminal set afresh, add the
states closest to or beyond the bounds, solve again, until a fresh sample of the design's own terminal set respects the
bounds. Then the period check of A (largest period with a one-step increase of Psi on the terminal set) and the grid
check of the ECGBF margin. Usage: TASK=blue_fwd_sq LAM=0.1 python study14.py <s> <kfac>
Writes k1/<task>_s14_s<s>_k<kfac>.txt and data/study14/<task>_s<s>_k<kfac>.json; existing designs are not touched."""
import os, sys, json, time, subprocess, numpy as np, sympy as sp, scipy.io as sio
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_top import CertArm, fk, jac, jacdot_qdot, mass, coriolis_qdot, _solve2, _mul2, phi_np, SPEC, TASK, y1, y2
SC = json.load(open(os.path.join(HERE, "data", f"scale_{TASK}.json"))); MU0, DFRAC, T = SC["MU0"], SC["DFRAC"], SC["T"]; U = np.array(SC["U_MAX"], float)
LAM0 = float(os.environ.get("LAM", "0.1")); C0 = CertArm(os.path.join(HERE, "k1", f"{TASK}_d8_l{os.environ.get('LAM', '0.1')}.txt"), mu=MU0)
s, kfac = float(sys.argv[1]), float(sys.argv[2]); lam, mu = s * LAM0, s * s * MU0; SIG, DEG, M_ADD, MAXIT = 4.0, 8, 250, 14
MONS = [(i, j) for i in range(DEG + 1) for j in range(DEG + 1 - i)]; NB = len(MONS); P_, Q_ = np.array(MONS).T
OUT = os.path.join(HERE, "data", "study14"); os.makedirs(OUT, exist_ok=True); TAG = f"{TASK}_s{s:g}_k{kfac:g}"; TMP = os.path.join(HERE, "build", "s14"); os.makedirs(TMP, exist_ok=True)
z1, z2 = sp.symbols("z1 z2")


def theta_of(k1_syms):
    out = []
    for k in k1_syms:
        P = sp.Poly(sp.expand(k.subs({y1: SIG * z1, y2: SIG * z2}, simultaneous=True) / SIG), z1, z2)
        assert all(sum(m) <= DEG for m in P.monoms()); out += [float(P.coeff_monomial(z1**i * z2**j)) for i, j in MONS]
    return np.array(out)


def k1_of(th):
    return [sp.expand(SIG * sum(sp.Float(c, 17) * (y1 / SIG)**i * (y2 / SIG)**j for c, (i, j) in zip(th[o * NB:(o + 1) * NB], MONS))) for o in range(2)]


def rows(x):
    """u^c(x) = a + Gu theta for the controller at (lam, mu); a (2, n), Gu (2, n, 2 NB)."""
    q, dq = x[:2], x[2:]; y = fk(q); J = jac(q); yd = _mul2(J, dq); n = x.shape[1]
    g = np.array([np.broadcast_to(np.asarray(v, float), (n,)) for v in C0.gpsi(*y)])
    W = np.einsum("ijn,jkn->ikn", mass(q), np.linalg.inv(J.transpose(2, 0, 1)).transpose(1, 2, 0))
    a = np.einsum("ijn,jn->in", W, mu * g - jacdot_qdot(q, dq) + 0.5 * lam * yd) + coriolis_qdot(q, dq)
    t1, t2 = y[0] / SIG, y[1] / SIG
    m = t1[:, None]**P_ * t2[:, None]**Q_; m1 = P_ * t1[:, None]**np.maximum(P_ - 1, 0) * t2[:, None]**Q_; m2 = Q_ * t1[:, None]**P_ * t2[:, None]**np.maximum(Q_ - 1, 0)
    Z = np.zeros_like(m); Phi = np.stack([np.hstack([SIG * m, Z]), np.hstack([Z, SIG * m]), np.hstack([m1, Z]), np.hstack([m2, Z]), np.hstack([Z, m1]), np.hstack([Z, m2])])
    F = np.zeros((2, 6, n)); F[0, 0] = F[1, 1] = -0.5 * lam; F[0, 2], F[0, 3], F[1, 4], F[1, 5] = yd[0], yd[1], yd[0], yd[1]
    return a, np.einsum("icn,cnk->ink", np.einsum("ijn,jcn->icn", W, F), Phi)


def ratio(c, x):
    return (np.abs(c.uc(x)) / U[:, None]).max(0)


def T_eps(c, x, lo=1e-5, hi=1.0):
    tau = c.uc(x); P0 = c.Psi(x)
    def ok(Tp):
        xn = c.zoh(x, tau, Tp, n=max(10, int(Tp / 0.0005))); d = c.Psi(xn) - P0; k = phi_np(*fk(xn[:2])) > 0; return (d[k].min() > 0) if k.any() else True
    if not ok(lo): return 0.0
    for _ in range(16):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def terminal(c, n, seed):
    x, _ = c.sample(DFRAC * (1 - c.shift), n_y=200000, seed=seed); x = x[:, phi_np(*fk(x[:2])) > 0]; return x[:, :n]


t_start = time.time(); A = C0.scaled(s); thA = theta_of(A.k1_sym)
xa = terminal(A, 30000, 1); aa, Ga = rows(xa[:, :2000]); err = np.abs(aa + Ga @ thA - A.uc(xa[:, :2000])).max()
rA = ratio(A, xa); log = dict(task=TASK, s=s, kfac=kfac, lam=lam, mu=mu, U=U.tolist(), affine_map_check=float(err), A=dict(max_ratio=float(rA.max()), share_over=float((rA > 1).mean())), it=[])
print(f"{TAG}: lambda {lam:.4f}, mu {mu:.4f} | affine map check {err:.1e} | design A at this scale: max |u_i|/U_i on its terminal set {rA.max():.3f} ({100 * (rA > 1).mean():.1f} % over)", flush=True)
XC = xa[:, np.argsort(-rA)[:M_ADD]]; th, dl, c = None, None, None
for it in range(MAXIT):
    a, G = rows(XC); Gm = np.vstack([G[0], -G[0], G[1], -G[1]]); bm = np.r_[U[0] - a[0], U[0] + a[0], U[1] - a[1], U[1] + a[1]]
    fin, fout = os.path.join(TMP, f"{TAG}_in.mat"), os.path.join(TMP, f"{TAG}_out.mat")
    sio.savemat(fin, dict(psi_str=SPEC["psi"], pmax=SPEC["pmax"], target=np.array(SPEC["target"], float), RT=SPEC["RT"], lam=lam, kmax=kfac * s, sigma=SIG, deg=DEG,
                          mons=np.array(MONS, float), G=Gm, b=bm, out_mat=fout))
    t0 = time.time(); r = subprocess.run(["matlab", "-batch", f"IN_MAT='{fin}'; run('{os.path.join(HERE, 'solve_k1_sampled14.m')}')"], capture_output=True, text=True)
    sol = [l for l in r.stdout.splitlines() if l.startswith("__SOL__")]
    if not sol: print(r.stdout[-1500:], r.stderr[-800:]); log["status"] = "matlab error"; break
    info = [float(v) for v in sol[0].split(",")[1:7]]; R = sio.loadmat(fout); th, dl = R["theta"].ravel(), float(R["dl"].ravel()[0])
    if not np.all(np.isfinite(th)) or info[0] or info[1] or info[3] < 0.5:
        print(f"  iter {it}: SOS program infeasible or failed (pinf {info[0]:g}, dinf {info[1]:g}, feasratio {info[3]:.2f})", flush=True); log["it"].append(dict(info=info, n_cons=int(XC.shape[1]))); log["status"] = "infeasible"; th = None; break
    c = CertArm(mu=mu, psi=C0.psi_sym, k1=k1_of(th), lam=lam, delta_s=max(dl, 0.0)); xb = terminal(c, 30000, 100 + it); rb = ratio(c, xb)
    if xb.shape[1] == 0: print(f"  iter {it}: terminal set of the new design is empty", flush=True); log["status"] = "terminal set empty"; th = None; break
    log["it"].append(dict(n_cons=int(XC.shape[1]), delta=dl, info=info, max_ratio=float(rb.max()), share_over=float((rb > 1).mean()), n_terminal=int(xb.shape[1]), sec=time.time() - t0))
    print(f"  iter {it}: {XC.shape[1]} sampled states | delta {dl:.2e}, feasratio {info[3]:.2f} | fresh terminal set ({xb.shape[1]} states): max ratio {rb.max():.4f}, over {100 * (rb > 1).mean():.2f} % ({time.time() - t0:.0f}s)", flush=True)
    if rb.max() <= 1.0: log["status"] = "converged"; break
    XC = np.hstack([XC, xb[:, np.argsort(-rb)[:M_ADD]]])
else:
    log["status"] = "max iterations"
if th is not None and c is not None:
    xv = terminal(c, 60000, 999); rv = ratio(c, xv); g = c.margin_on_grid(n=601); te = T_eps(c, terminal(c, 8000, 7))
    kf = os.path.join(HERE, "k1", f"{TASK}_s14_s{s:g}_k{kfac:g}.txt")
    open(kf, "w").write(f"{c.k1_sym[0]}\n{c.k1_sym[1]}\n{lam!r}\n{max(dl, 0.0)!r}\n{sp.sympify(SPEC['psi'])}\n{SPEC['pmax']!r}\n")
    log.update(k1_file=os.path.basename(kf), verify=dict(n=int(xv.shape[1]), max_ratio=float(rv.max()), share_over=float((rv > 1).mean()), grid_min_margin=g["min_margin"], kmax=g["kmax"],
               area_share=g["area_share"], T_eps_ms=1e3 * te, period_ok_at_T=bool(te >= T), period_ok_with_08=bool(0.8 * te >= T)))
    print(f"  verify on {xv.shape[1]} fresh terminal states: max ratio {rv.max():.4f} ({100 * (rv > 1).mean():.2f} % over) | grid min ECGBF margin {g['min_margin']:+.2e}, max|k1_i| {g['kmax']:.3f} | "
          f"T_eps {1e3 * te:.2f} ms (A at this scale {SC['T_eps1_ms'] / s:.2f} ms) | status {log['status']} | total {time.time() - t_start:.0f}s", flush=True)
log["A"].update(T_eps_ms=SC["T_eps1_ms"] / s, kmax=s * C0.margin_on_grid(n=601)["kmax"]); log["total_s"] = time.time() - t_start
json.dump(log, open(os.path.join(OUT, f"{TAG}.json"), "w"), indent=1)
