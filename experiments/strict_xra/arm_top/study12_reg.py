"""Study: the S-procedure bound |k1_i| <= K in (12) against a coefficient-norm penalty in the objective (lambda fixed).
For each penalty design (OBJ inf: delta + gamma max|theta_k|; OBJ one: delta + gamma sum|theta_k|; no K bound) the
realised K_P = max |k1_i| on {psi >= 0} is read off, and the hard-bound program (12) is solved again with K = K_P at the
same lambda. The two designs are then compared: slack delta, field difference on the safe set, grid margin, area share
of the terminal set, period threshold T_eps(1), and the same quantities at the normalised speed max|k1| = 1 (exact
scaling by s = 1 / K_P). Usage: TASK=red_fwd_sq LAM=0.1 python study12_reg.py. Writes data/study12_reg/<task>.json."""
import os, sys, json, time, subprocess, numpy as np, sympy as sp, scipy.io as sio
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_top import CertArm, fk, phi_np, SPEC, TASK, BOX, y1, y2
SC = json.load(open(os.path.join(HERE, "data", f"scale_{TASK}.json"))); MU0, DFRAC = SC["MU0"], SC["DFRAC"]; LAM = float(os.environ.get("LAM", "0.1"))
C_A = CertArm(os.path.join(HERE, "k1", f"{TASK}_d8_l{os.environ.get('LAM', '0.1')}.txt"), mu=MU0); SIG, DEG = 4.0, 8
MONS = [(i, j) for i in range(DEG + 1) for j in range(DEG + 1 - i)]; NB = len(MONS); z1, z2 = sp.symbols("z1 z2")
OUT = os.path.join(HERE, "data", "study12_reg"); os.makedirs(OUT, exist_ok=True); TMP = os.path.join(HERE, "build", "s12"); os.makedirs(TMP, exist_ok=True)
G1, G2 = np.meshgrid(np.linspace(*BOX[0], 401), np.linspace(*BOX[1], 401)); SAFE = (C_A.psi(G1, G2) >= 0) & (phi_np(G1, G2) >= 0)


def theta_of(k1_syms):
    out = []
    for k in k1_syms:
        P = sp.Poly(sp.expand(k.subs({y1: SIG * z1, y2: SIG * z2}, simultaneous=True) / SIG), z1, z2); out += [float(P.coeff_monomial(z1**i * z2**j)) for i, j in MONS]
    return np.array(out)


def k1_of(th):
    return [sp.expand(SIG * sum(sp.Float(c, 17) * (y1 / SIG)**i * (y2 / SIG)**j for c, (i, j) in zip(th[o * NB:(o + 1) * NB], MONS))) for o in range(2)]


def solve(obj, gamma, kmax, tag):
    fin, fout = os.path.join(TMP, f"{tag}_in.mat"), os.path.join(TMP, f"{tag}_out.mat")
    sio.savemat(fin, dict(psi_str=SPEC["psi"], pmax=SPEC["pmax"], target=np.array(SPEC["target"], float), RT=SPEC["RT"], lam=LAM, kmax=float(kmax), gamma=float(gamma), obj=obj, sigma=SIG, deg=DEG, mons=np.array(MONS, float), out_mat=fout))
    t0 = time.time(); r = subprocess.run(["matlab", "-batch", f"IN_MAT='{fin}'; run('{os.path.join(HERE, 'solve_k1_reg12.m')}')"], capture_output=True, text=True)
    sol = [l for l in r.stdout.splitlines() if l.startswith("__SOL__")]
    if not sol: print(r.stdout[-1200:], r.stderr[-600:]); return None
    R = sio.loadmat(fout); info = [float(v) for v in sol[0].split(",")[1:7]]
    return dict(theta=R["theta"].ravel(), delta=float(R["dl"].ravel()[0]), info=info, sec=time.time() - t0)


def field(c):
    k = c.k1(G1, G2); return np.stack([np.broadcast_to(np.asarray(a, float), G1.shape) for a in k])


def T_eps(c, x, lo=1e-4, hi=3.0):
    tau = c.uc(x); P0 = c.Psi(x)
    def ok(Tp):
        xn = c.zoh(x, tau, Tp, n=max(10, int(Tp / 0.0005))); d = c.Psi(xn) - P0; k = phi_np(*fk(xn[:2])) > 0; return (d[k].min() > 0) if k.any() else True
    if not ok(lo): return 0.0
    for _ in range(16):
        mid = np.sqrt(lo * hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def evaluate(c, name):
    g = c.margin_on_grid(n=601); K = g["kmax"]; cn = c.scaled(1.0 / K)                         # normalised to max|k1| = 1
    x, _ = cn.sample(DFRAC * (1 - cn.shift), n_y=200000, seed=11); x = x[:, phi_np(*fk(x[:2])) > 0][:, :8000]
    te = T_eps(cn, x); return dict(name=name, K=K, delta=c.ds, shift=c.shift, grid_margin=g["min_margin"], area_share=g["area_share"], lam_norm=cn.lam, T_eps1_norm_ms=1e3 * te)


def compare(ca, cb):
    fa, fb = field(ca), field(cb); d = np.sqrt(((fa - fb)**2).sum(0))[SAFE]; na = np.sqrt((fa**2).sum(0))[SAFE]
    cos = ((fa * fb).sum(0) / np.maximum(np.sqrt((fa**2).sum(0) * (fb**2).sum(0)), 1e-12))[SAFE]
    return dict(max_diff_rel=float(d.max() / na.max()), mean_diff_rel=float(d.mean() / na.mean()), min_cos=float(cos.min()), mean_cos=float(cos.mean()))


t0 = time.time(); thA = theta_of(C_A.k1_sym); print(f"{TASK}, lambda {LAM}: design A (hard bound K = 1): max|theta| {np.abs(thA).max():.3f}, sum|theta| {np.abs(thA).sum():.2f}, delta/lambda {C_A.shift:.2e}", flush=True)
EA = evaluate(C_A, "A hard K=1"); print(f"  A: K {EA['K']:.3f}, grid margin {EA['grid_margin']:+.2e}, area {EA['area_share']:.3f}, T_eps(1) at unit speed {EA['T_eps1_norm_ms']:.1f} ms", flush=True)
rows = [dict(design=EA, matched=None, vs_matched=None, vs_A=None)]
specs = [("delta", 0.0, 0.0)] + [("inf", g, 0.0) for g in (1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1)] + [("one", g, 0.0) for g in (1e-5, 1e-4, 1e-3, 1e-2)]
for obj, gamma, kmax in specs:
    tag = f"{TASK}_{obj}_g{gamma:g}"; P = solve(obj, gamma, kmax, tag)
    if P is None or not np.all(np.isfinite(P["theta"])): print(f"{tag}: solver failed", flush=True); continue
    cP = CertArm(mu=MU0, psi=C_A.psi_sym, k1=k1_of(P["theta"]), lam=LAM, delta_s=max(P["delta"], 0.0)); EP = evaluate(cP, tag)
    if EP["K"] > 50: print(f"{tag}: max|k1| {EP['K']:.3g}, delta {P['delta']:.2e} | unbounded, no matched design", flush=True); rows.append(dict(design=EP, matched=None, vs_matched=None, vs_A=compare(C_A, cP))); continue
    H = solve("delta", 0.0, EP["K"], tag + "_matched"); cH = CertArm(mu=MU0, psi=C_A.psi_sym, k1=k1_of(H["theta"]), lam=LAM, delta_s=max(H["delta"], 0.0)); EH = evaluate(cH, tag + " matched hard K")
    vs = compare(cH, cP); vA = compare(C_A, cP.scaled(1.0 / EP["K"]))
    print(f"{tag}: K_P {EP['K']:.3f}, delta {EP['delta']:.2e} (hard K={EP['K']:.3f}: {EH['delta']:.2e}) | field vs matched: max diff {100*vs['max_diff_rel']:.1f} %, min cos {vs['min_cos']:.3f} | "
          f"area {EP['area_share']:.3f} vs {EH['area_share']:.3f} | T_eps(1) at unit speed {EP['T_eps1_norm_ms']:.1f} vs {EH['T_eps1_norm_ms']:.1f} ms | vs A (unit speed): max diff {100*vA['max_diff_rel']:.1f} %, min cos {vA['min_cos']:.3f} | {P['sec']:.0f}+{H['sec']:.0f}s", flush=True)
    rows.append(dict(design=EP, matched=EH, vs_matched=vs, vs_A=vA))
json.dump(dict(task=TASK, lam=LAM, A_theta_max=float(np.abs(thA).max()), rows=rows, total_s=time.time() - t0), open(os.path.join(OUT, f"{TASK}.json"), "w"), indent=1)
print(f"total {time.time() - t0:.0f}s")
