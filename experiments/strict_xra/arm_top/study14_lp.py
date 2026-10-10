"""Study C: program (14) solved without SOS, as a linear program in the coefficients theta of k1 (all constraints sampled).
At a time scale s, lambda = s lambda0 and mu = s^2 mu0 are fixed (as for the scaled design A and the SOS design B of
study14.py). Variables theta (k1 in monomials of y/4 up to degree D) and gamma; maximise gamma subject to
    grad psi(y_i) . k1(y_i) - lambda psi(y_i) >= gamma    at samples of the closure of X_S minus X_T (grid and boundary)
    |k1_o(y_l)| <= s K0                                   at samples of {psi >= 0}
    |u^c_o(x_j)| <= U_o                                   at states of the design's own terminal set X_RA^delta
(u^c is affine in theta at a fixed state, Theorem 1). Multiresolution exchange: grids of spacing 0.04 ... 0.0025, at
each level the LP is solved, the level's grid and a fresh terminal sample are checked and the violating points added,
until the level passes; final checks on spacings 0.00125 and 0.000625 without adding. Then the period check.
Usage: TASK=blue_fwd_sq LAM=0.1 python study14_lp.py <s> [D]. Writes k1/<tag>.txt and data/study14lp/<tag>.json (tag = <task>_lp_s<s>_D<D>[_<OBJ>])."""
import os, sys, json, time, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import linprog
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_top import CertArm, fk, jac, jacdot_qdot, mass, coriolis_qdot, _mul2, phi_np, SPEC, TASK, BOX, y1, y2
SC = json.load(open(os.path.join(HERE, "data", f"scale_{TASK}.json"))); MU0, DFRAC, T = SC["MU0"], SC["DFRAC"], SC["T"]; U = np.array(SC["U_MAX"], float)
LAM0 = float(os.environ.get("LAM", "0.1")); C0 = CertArm(os.path.join(HERE, "k1", f"{TASK}_d8_l{os.environ.get('LAM', '0.1')}.txt"), mu=MU0)
s = float(sys.argv[1]); DEG = int(sys.argv[2]) if len(sys.argv) > 2 else 8; lam, mu, K = s * LAM0, s * s * MU0, s * 1.0; SIG, M_ADD, CB = 4.0, 250, 200.0
MONS = [(i, j) for i in range(DEG + 1) for j in range(DEG + 1 - i)]; NB = len(MONS); P_, Q_ = np.array(MONS).T
OUT = os.path.join(HERE, "data", "study14lp"); os.makedirs(OUT, exist_ok=True); TAG = f"{TASK}_lp_s{s:g}_D{DEG}" + ("" if os.environ.get("OBJ", "gamma") == "gamma" else "_" + os.environ["OBJ"])
LEVELS, CHECKS = (0.04, 0.02, 0.01, 0.005, 0.0025), (0.00125, 0.000625)
OBJ, RHO = os.environ.get("OBJ", "gamma"), float(os.environ.get("RHO", "0.25"))   # gamma: max margin; umax, jac, hess, coef: max margin, then with gamma >= RHO gamma* minimise the largest torque ratio, Jacobian entry, Hessian entry of k1, or coefficient
psi_n = sp.lambdify((y1, y2), C0.psi_sym, "numpy"); gpsi_n = sp.lambdify((y1, y2), [sp.diff(C0.psi_sym, v) for v in (y1, y2)], "numpy")


def mono(A, B):
    t1, t2 = A / SIG, B / SIG; return t1[:, None]**P_ * t2[:, None]**Q_


def grid_points(h, chunk=400):
    """Grid of spacing h over the box, split into row blocks; yields (A, B, psi, phi) of the points with psi >= 0."""
    g1 = np.arange(BOX[0][0], BOX[0][1] + h, h); g2 = np.arange(BOX[1][0], BOX[1][1] + h, h)
    for i in range(0, len(g2), chunk):
        A, B = np.meshgrid(g1, g2[i:i + chunk]); A, B = A.ravel(), B.ravel(); p = psi_n(A, B); keep = p >= 0
        yield A[keep], B[keep], p[keep], phi_np(A[keep], B[keep])


def boundary_points(h):
    """Points on the zero levels of psi (corridor edge) and of phi (target circle) inside the closure of X_S minus X_T."""
    g1 = np.arange(BOX[0][0], BOX[0][1] + h, h); g2 = np.arange(BOX[1][0], BOX[1][1] + h, h); G1, G2 = np.meshgrid(g1, g2)
    pts = [np.vstack(cs.allsegs[0]) for cs in [plt.contour(G1, G2, psi_n(G1, G2), levels=[0])] if cs.allsegs[0]]; plt.close("all")
    a = np.linspace(0, 2 * np.pi, int(2 * np.pi * SPEC["RT"] / h) + 8); c = np.array(SPEC["target"]) + SPEC["RT"] * np.c_[np.cos(a), np.sin(a)]
    P = np.vstack(pts + [c[psi_n(c[:, 0], c[:, 1]) >= 0]]); return P[:, 0], P[:, 1]


def margin(th, A, B):
    Mb = mono(A, B); g = np.array([np.broadcast_to(np.asarray(v, float), A.shape) for v in gpsi_n(A, B)]); k = SIG * np.stack([Mb @ th[:NB], Mb @ th[NB:]])
    return (g * k).sum(0) - lam * psi_n(A, B), np.abs(k).max(0)


def rows_u(x):
    q, dq = x[:2], x[2:]; y = fk(q); J = jac(q); yd = _mul2(J, dq); n = x.shape[1]
    g = np.array([np.broadcast_to(np.asarray(v, float), (n,)) for v in gpsi_n(*y)])
    W = np.einsum("ijn,jkn->ikn", mass(q), np.linalg.inv(J.transpose(2, 0, 1)).transpose(1, 2, 0))
    a = np.einsum("ijn,jn->in", W, mu * g - jacdot_qdot(q, dq) + 0.5 * lam * yd) + coriolis_qdot(q, dq)
    t1, t2 = y[0] / SIG, y[1] / SIG
    m = t1[:, None]**P_ * t2[:, None]**Q_; m1 = P_ * t1[:, None]**np.maximum(P_ - 1, 0) * t2[:, None]**Q_; m2 = Q_ * t1[:, None]**P_ * t2[:, None]**np.maximum(Q_ - 1, 0)
    Z = np.zeros_like(m); Phi = np.stack([np.hstack([SIG * m, Z]), np.hstack([Z, SIG * m]), np.hstack([m1, Z]), np.hstack([m2, Z]), np.hstack([Z, m1]), np.hstack([Z, m2])])
    F = np.zeros((2, 6, n)); F[0, 0] = F[1, 1] = -0.5 * lam; F[0, 2], F[0, 3], F[1, 4], F[1, 5] = yd[0], yd[1], yd[0], yd[1]
    return a, np.einsum("icn,cnk->ink", np.einsum("ijn,jcn->icn", W, F), Phi)


def solve(LA, LB, KA, KB, XU):
    Mb = mono(LA, LB); g = np.array([np.broadcast_to(np.asarray(v, float), LA.shape) for v in gpsi_n(LA, LB)]); z1 = np.ones((len(LA), 1))
    R = [np.hstack([-SIG * g[0][:, None] * Mb, -SIG * g[1][:, None] * Mb, z1])]; b = [-lam * psi_n(LA, LB)]
    Mk = SIG * mono(KA, KB); Zk = np.zeros_like(Mk); zk = np.zeros((len(KA), 1))
    R += [np.hstack([Mk, Zk, zk]), np.hstack([-Mk, Zk, zk]), np.hstack([Zk, Mk, zk]), np.hstack([Zk, -Mk, zk])]; b += [K * np.ones(4 * len(KA))]
    if XU.shape[1]:
        a, G = rows_u(XU); zu = np.zeros((XU.shape[1], 1))
        for o in range(2): R += [np.hstack([G[o] / U[o], zu]), np.hstack([-G[o] / U[o], zu])]; b += [1 - a[o] / U[o], 1 + a[o] / U[o]]
    t0 = time.time(); R, b = np.vstack(R), np.concatenate(b)
    r = linprog(np.r_[np.zeros(2 * NB), -1.0], A_ub=R, b_ub=b, bounds=[(-CB, CB)] * (2 * NB) + [(None, 1.0)], method="highs")
    if r.x is None: return None, None, time.time() - t0
    if OBJ == "gamma": return r.x[:-1], r.x[-1], time.time() - t0
    g0 = RHO * r.x[-1]; nr = R.shape[0]                       # stage 2: gamma fixed at RHO gamma*, new variable t in the last column
    R2 = np.hstack([R[:, :-1], np.zeros((nr, 1))]); b2 = b - R[:, -1] * g0
    if OBJ == "l2":                                           # minimum Euclidean norm of theta (second-order cone, Mosek)
        import mosek.fusion as mf
        with mf.Model() as Mo:
            th_ = Mo.variable(2 * NB, mf.Domain.inRange(-CB, CB)); t_ = Mo.variable(1, mf.Domain.greaterThan(0.0))
            Mo.constraint(mf.Expr.mul(mf.Matrix.dense(R[:, :-1]), th_), mf.Domain.lessThan(b2)); Mo.constraint(mf.Expr.vstack(t_, th_), mf.Domain.inQCone())
            Mo.objective(mf.ObjectiveSense.Minimize, t_); Mo.solve()
            if Mo.getPrimalSolutionStatus() not in (mf.SolutionStatus.Optimal,): return r.x[:-1], r.x[-1], time.time() - t0
            return np.array(th_.level()), g0, time.time() - t0
    if OBJ == "umax":
        a, G = rows_u(XU); ex = [np.hstack([G[o] / U[o], -np.ones((XU.shape[1], 1))]) for o in range(2)] + [np.hstack([-G[o] / U[o], -np.ones((XU.shape[1], 1))]) for o in range(2)]
        eb = [-a[0] / U[0], -a[1] / U[1], a[0] / U[0], a[1] / U[1]]
    elif OBJ == "jac":
        t1, t2 = LA / SIG, LB / SIG; m1 = P_ * t1[:, None]**np.maximum(P_ - 1, 0) * t2[:, None]**Q_; m2 = Q_ * t1[:, None]**P_ * t2[:, None]**np.maximum(Q_ - 1, 0); Z = np.zeros_like(m1); o1 = -np.ones((len(LA), 1))
        ex = [np.hstack([sg * blk, o1]) for blk in (np.hstack([m1, Z]), np.hstack([m2, Z]), np.hstack([Z, m1]), np.hstack([Z, m2])) for sg in (1, -1)]; eb = [np.zeros(len(LA))] * 8
    elif OBJ == "hess":
        t1, t2 = LA / SIG, LB / SIG; d = lambda a, b: (P_ * (P_ - 1) if a == 2 else P_ if a == 1 else 1) * (Q_ * (Q_ - 1) if b == 2 else Q_ if b == 1 else 1) * t1[:, None]**np.maximum(P_ - a, 0) * t2[:, None]**np.maximum(Q_ - b, 0)
        H = [d(2, 0), d(1, 1), d(0, 2)]; Z = np.zeros_like(H[0]); o1 = -np.ones((len(LA), 1))
        ex = [np.hstack([sg * blk, o1]) for h_ in H for blk in (np.hstack([h_, Z]), np.hstack([Z, h_])) for sg in (1, -1)]; eb = [np.zeros(len(LA))] * 12
    else:                                                    # coef: minimise the largest |theta_k| (coefficient size)
        I = np.eye(2 * NB); o1 = -np.ones((2 * NB, 1)); ex = [np.hstack([I, o1]), np.hstack([-I, o1])]; eb = [np.zeros(2 * NB)] * 2
    r2 = linprog(np.r_[np.zeros(2 * NB), 1.0], A_ub=np.vstack([R2] + ex), b_ub=np.concatenate([b2] + eb), bounds=[(-CB, CB)] * (2 * NB) + [(0, None)], method="highs")
    if r2.x is None: return r.x[:-1], r.x[-1], time.time() - t0
    return r2.x[:-1], g0, time.time() - t0


def cert_of(th):
    k1 = [sp.expand(SIG * sum(sp.Float(c, 17) * (y1 / SIG)**i * (y2 / SIG)**j for c, (i, j) in zip(th[o * NB:(o + 1) * NB], MONS))) for o in range(2)]
    return CertArm(mu=mu, psi=C0.psi_sym, k1=k1, lam=lam, delta_s=1e-9)


def terminal(c, n, seed):
    x, _ = c.sample(DFRAC * (1 - c.shift), n_y=200000, seed=seed); x = x[:, phi_np(*fk(x[:2])) > 0]; return x[:, :n]


def T_eps(c, x, lo=1e-5, hi=1.0):
    tau = c.uc(x); P0 = c.Psi(x)
    def ok(Tp):
        xn = c.zoh(x, tau, Tp, n=max(10, int(Tp / 0.0005))); d = c.Psi(xn) - P0; k = phi_np(*fk(xn[:2])) > 0; return (d[k].min() > 0) if k.any() else True
    if not ok(lo): return 0.0
    for _ in range(16):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


t_start = time.time(); A0 = C0.scaled(s); xa = terminal(A0, 30000, 1); ra = (np.abs(A0.uc(xa)) / U[:, None]).max(0)
LA, LB = [], []
for A_, B_, p_, f_ in grid_points(0.1):
    LA.append(A_[f_ >= 0]); LB.append(B_[f_ >= 0])
LA, LB = np.concatenate(LA), np.concatenate(LB); KA, KB = LA.copy(), LB.copy(); XU = xa[:, np.argsort(-ra)[:M_ADD]]
log = dict(task=TASK, s=s, deg=DEG, lam=lam, mu=mu, K=K, U=U.tolist(), A=dict(max_ratio=float(ra.max())), levels=[]); n_lp, t_lp, th = 0, 0.0, None
print(f"{TAG}: lambda {lam:.4f}, mu {mu:.4f}, K {K:.3f}, {2 * NB} coefficients | design A at this scale: max ratio {ra.max():.3f}", flush=True)
for h in LEVELS:
    bA, bB = boundary_points(h)
    for rnd in range(25):
        th, gam, dt = solve(LA, LB, KA, KB, XU); n_lp += 1; t_lp += dt
        if th is None: print(f"  h {h:g}: LP infeasible", flush=True); break
        worst, nbad, kmx, Ladd, Kadd = np.inf, 0, 0.0, [], []
        for A_, B_, p_, f_ in list(grid_points(h)) + [(bA, bB, psi_n(bA, bB), phi_np(bA, bB))]:
            if A_.size == 0: continue
            mg, kn = margin(th, A_, B_); sel = f_ >= -1e-12; bad = sel & (mg < 0.5 * gam if gam > 0 else mg < 0)
            if sel.any(): worst = min(worst, mg[sel].min())
            nbad += bad.sum(); kmx = max(kmx, kn.max()); Ladd.append(np.c_[A_[bad], B_[bad]]); kb = kn > K * 1.005; Kadd.append(np.c_[A_[kb], B_[kb]])
        Ladd, Kadd = np.vstack(Ladd), np.vstack(Kadd); c = cert_of(th); xb = terminal(c, 30000, 100 + n_lp); rb = (np.abs(c.uc(xb)) / U[:, None]).max(0) if xb.shape[1] else np.array([np.inf])
        ok = worst >= 0 and kmx <= K * 1.005 and rb.max() <= 1.0
        print(f"  h {h:g} round {rnd}: {len(LA)} Lie, {len(KA)} K, {XU.shape[1]} input samples | gamma {gam:+.3e} | grid min margin {worst:+.3e} ({nbad} below), max|k1| {kmx:.3f} | terminal ratio {rb.max():.3f} | LP {dt:.2f}s", flush=True)
        if ok: break
        if len(Ladd): j = np.argsort(margin(th, Ladd[:, 0], Ladd[:, 1])[0])[:3000]; LA, LB = np.r_[LA, Ladd[j, 0]], np.r_[LB, Ladd[j, 1]]
        if len(Kadd): j = np.argsort(-margin(th, Kadd[:, 0], Kadd[:, 1])[1])[:2000]; KA, KB = np.r_[KA, Kadd[j, 0]], np.r_[KB, Kadd[j, 1]]
        if rb.max() > 1.0: XU = np.hstack([XU, xb[:, np.argsort(-rb)[:M_ADD]]])
    log["levels"].append(dict(h=h, rounds=rnd + 1, gamma=float(gam) if th is not None else None, min_margin=float(worst) if th is not None else None, ok=bool(th is not None and ok)))
    if th is None: break
if th is not None:
    chk = {}
    for h in CHECKS:
        bA, bB = boundary_points(h); worst, kmx, at = np.inf, 0.0, None
        for A_, B_, p_, f_ in list(grid_points(h, chunk=200)) + [(bA, bB, psi_n(bA, bB), phi_np(bA, bB))]:
            if A_.size == 0: continue
            mg, kn = margin(th, A_, B_); sel = f_ >= -1e-12
            if sel.any() and mg[sel].min() < worst: worst = mg[sel].min(); i = np.where(sel)[0][np.argmin(mg[sel])]; at = (float(A_[i]), float(B_[i]), float(p_[i]))
            kmx = max(kmx, kn.max())
        chk[str(h)] = dict(min_margin=float(worst), at=at, kmax=float(kmx)); print(f"  check on spacing {h:g}: min margin {worst:+.3e} at y = ({at[0]:.3f}, {at[1]:.3f}), psi {at[2]:.4f}, max|k1| {kmx:.3f}", flush=True)
    c = cert_of(th); xv = terminal(c, 60000, 999); rv = (np.abs(c.uc(xv)) / U[:, None]).max(0); te = T_eps(c, terminal(c, 8000, 7))
    kf = os.path.join(HERE, "k1", f"{TAG}.txt"); open(kf, "w").write(f"{c.k1_sym[0]}\n{c.k1_sym[1]}\n{lam!r}\n{1e-9!r}\n{sp.sympify(SPEC['psi'])}\n{SPEC['pmax']!r}\n")
    valid = all(v["min_margin"] >= 0 and v["kmax"] <= K * 1.01 for v in chk.values()) and rv.max() <= 1.0
    log.update(k1_file=os.path.basename(kf), checks=chk, n_lp=n_lp, lp_seconds=t_lp, n_samples=dict(lie=int(len(LA)), k=int(len(KA)), input=int(XU.shape[1])),
               verify=dict(n=int(xv.shape[1]), max_ratio=float(rv.max()), T_eps_ms=1e3 * te, period_ok_with_08=bool(0.8 * te >= T), valid=bool(valid)), status="done")
    print(f"  verify: {xv.shape[1]} fresh terminal states, max ratio {rv.max():.4f} | T_eps {1e3 * te:.2f} ms (A at this scale {SC['T_eps1_ms'] / s:.2f}) | valid {valid}, period ok {0.8 * te >= T} | "
          f"{n_lp} LPs ({t_lp:.1f} s in HiGHS) | total {time.time() - t_start:.0f}s", flush=True)
else:
    log["status"] = "infeasible"
log["total_s"] = time.time() - t_start; json.dump(log, open(os.path.join(OUT, f"{TAG}.json"), "w"), indent=1)
