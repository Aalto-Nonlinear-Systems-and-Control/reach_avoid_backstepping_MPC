"""Constants of Lemma 2 / Theorem 2 (TAC revision) for the Dubins car.

All sup/Lipschitz quantities are estimated by dense sampling over the closure of
X_RA = {Psi >= 0} intersected with a state box (trajectories of Lemma 2 stay in
X_RA, so constants on X_RA suffice).  Sampled maxima are inflated by `infl`.
"""
import numpy as np
import sympy as sp


def dubins_constants(n_samples=400_000, seed=0, infl=1.1, lam=None,
                     lo=(-2.2, -2.2, 2 * np.pi / 3, 0.1), hi=(2.2, 2.2, 4 * np.pi / 3, 1.0)):
    from sop_bounded_control_dubins_car_result import u_opt, certificate_opt

    x1, x2, th, v = sp.symbols("x1 x2 th v")
    X = [x1, x2, th, v]
    Psi = sp.sympify(certificate_opt)
    u = [sp.sympify(e) for e in u_opt]
    f = sp.Matrix([v * sp.cos(th), v * sp.sin(th), 0, 0])
    G = sp.Matrix([[0, 0], [0, 0], [1, 0], [0, 1]])

    gradPsi = sp.Matrix([sp.diff(Psi, xi) for xi in X])
    LgPsi = [(gradPsi.T * G[:, j])[0] for j in range(2)]
    LfPsi = (gradPsi.T * f)[0]
    grad_u = [[sp.diff(uj, xi) for xi in X] for uj in u]
    Jf = f.jacobian(X)
    F = sp.lambdify(X, [Psi, LfPsi, *LgPsi, *u, *sum(grad_u, []), *list(Jf)], "numpy")

    rng = np.random.default_rng(seed)
    lo = np.asarray(lo, float); hi = np.asarray(hi, float)
    S = lo[:, None] + (hi - lo)[:, None] * rng.random((4, n_samples))
    vals = F(*S)
    bc = lambda a: np.broadcast_to(np.asarray(a, float), S.shape[1:])
    vals = [bc(a) for a in vals]
    psi, lf, lg1, lg2, u1, u2 = vals[:6]
    gu = np.array(vals[6:14]).reshape(2, 4, -1)
    jf = np.array(vals[14:30]).reshape(4, 4, -1)
    tgt = (S[1] ** 2 + 4 * (S[0] + 1.7) ** 2 - 0.4)  # <=0 inside target
    ra = psi >= 0
    out = ra & (tgt > 0)

    c = {}
    c["n_RA"] = int(ra.sum())
    c["M_Psi"] = float(psi[ra].max()) * infl
    c["M_g"] = [float(np.abs(lg1[ra]).max()) * infl, float(np.abs(lg2[ra]).max()) * infl]
    c["l_uc"] = [float(np.linalg.norm(gu[j][:, ra], axis=0).max()) * infl for j in range(2)]
    c["M_uc"] = [float(np.abs(u1[ra]).max()), float(np.abs(u2[ra]).max())]
    c["l_f"] = float(max(np.linalg.norm(jf[:, :, i], 2) for i in np.where(ra)[0][:: max(1, ra.sum() // 20000)])) * infl
    c["l_g"] = [0.0, 0.0]
    c["l_sys"] = c["l_f"] + sum(a * b for a, b in zip(c["l_g"], c["M_uc"]))
    vel = np.sqrt((S[3] * np.cos(S[2])) ** 2 + (S[3] * np.sin(S[2])) ** 2 + u1 ** 2 + u2 ** 2)
    c["M_v"] = float(vel[ra].max()) * infl
    psidot = lf + lg1 * u1 + lg2 * u2
    m = out & (psi > 1e-3 * c["M_Psi"])
    c["lambda_est"] = float((psidot[m] / psi[m]).min())
    c["lambda"] = float(lam) if lam is not None else c["lambda_est"]
    c["psi_min_on_target"] = float(psi[(tgt <= 0) & (np.abs(S[3]) <= 1)].min()) if ((tgt <= 0)).any() else float("nan")
    return c


def beta(c, T):
    a = sum(m * l for m, l in zip(c["M_g"], c["l_uc"]))
    gamma = a * (np.exp(c["l_sys"] * T) - 1) / c["l_sys"]
    return gamma * c["M_v"]


def T_star(c, delta):
    """Largest T with beta(T) < lambda * delta."""
    a = sum(m * l for m, l in zip(c["M_g"], c["l_uc"])) * c["M_v"]
    return np.log(1 + c["l_sys"] * c["lambda"] * delta / a) / c["l_sys"]


def omega_min(c, T, delta, sup_ell):
    return sup_ell / ((np.exp(c["lambda"] * T) - 1) * (delta - beta(c, T) / c["lambda"]))
