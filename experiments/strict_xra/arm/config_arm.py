"""Final setting of the two-link arm (option A: the torque is held over the period, no inner-loop compensation)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_arm import CertArm, fk, ik, phi_np, gravity, f_np, jac, TARGET_POS
import ra_mpc_arm as ra

K1_FILE = os.path.join(HERE, "k1", "ARM_l0.04.txt"); MU_BASE, DFRAC = 0.1, 0.5
T, N, U_MAX = 0.0015, 133, (120.0, 40.0)
T_HOR = T * N; T_MAX_MPC, T_MAX_URA = 8.0, 60.0
SAFE_MARGIN, NSUB = 0.005, 1
_sj = os.path.join(HERE, "data", "scale_arm.json")
SCALE = float(os.environ.get("SCALE", json.load(open(_sj))["s"] if os.path.exists(_sj) else 2.5))
URA_SCALE = float(os.environ.get("URA_SCALE", 5.0))
TAG = os.environ.get("TAG", ""); DATA, FIGS = os.path.join(HERE, "data" + TAG), os.path.join(HERE, "figs" + TAG)
os.makedirs(DATA, exist_ok=True); os.makedirs(FIGS, exist_ok=True)
SEED_POS, N_INIT = 5, 150

C = CertArm(K1_FILE, mu=MU_BASE).scaled(SCALE); M_PSI = 1.0 - C.shift; DELTA = DFRAC * M_PSI


def positions(n=N_INIT, seed=SEED_POS):
    """Band positions outside the target, uniform in the output space."""
    r = np.random.default_rng(seed); Y = np.array([r.uniform(-4, 8.5, 60000), r.uniform(-3, 3, 60000)])
    keep = (np.broadcast_to(C.psi(*Y), (Y.shape[1],)) > 0) & (phi_np(*Y) > 0)
    return Y[:, keep][:, :n]


def initial_set():
    """At rest, elbow-down."""
    Y = positions(); q = ik(Y); return np.vstack([q, np.zeros_like(q)])


def ura_run(x0, T_max=T_MAX_URA, scale=URA_SCALE, dense=False):
    from scipy.integrate import solve_ivp
    c = CertArm(K1_FILE, mu=MU_BASE).scaled(scale)
    def rhs(t, z):
        x = z[:4]; u = c.uc(x.reshape(4, 1))[:, 0]; dy = fk(x[:2]) - TARGET_POS
        return np.concatenate([f_np(x.reshape(4, 1), u)[:, 0], [ra.Q_Y * dy @ dy + ra.R_U * u @ u]])
    e1 = lambda t, z: float(phi_np(*fk(z[:2]))); e1.terminal = True; e1.direction = -1
    e2 = lambda t, z: float(c.psi(*fk(z[:2]))); e2.terminal = True; e2.direction = -1
    s = solve_ivp(rhs, (0, T_max), np.append(np.asarray(x0, float), 0.0), events=[e1, e2], rtol=1e-7, atol=1e-9, max_step=0.01)
    u = c.uc(s.y[:4]); st = "reached" if s.t_events[0].size else "unsafe" if s.t_events[1].size else "timeout"
    out = dict(status=st, t=float(s.t[-1]), upeak=np.abs(u).max(1), cost=float(s.y[4, -1]), psi0=float(c.Psi(s.y[:4, :1])[0]),
               t_viol=float(np.sum(np.diff(s.t) * (np.abs(u[:, 1:]) > np.asarray(U_MAX)[:, None]).any(0))))
    if dense: out.update(tt=s.t, x=s.y[:4].T, u=u.T, y=fk(s.y[:2]).T)
    return out


def uc_zoh(X0, dt=T, T_max=T_MAX_URA):
    """Sampled-data terminal controller alone (held torque), vectorised. status 1 reached, 2 unsafe, 0 timeout."""
    X = X0.copy(); n = X.shape[1]; st = np.zeros(n, int); tt = np.full(n, np.nan); dmin = np.full(n, np.inf); P = C.Psi(X)
    nsub = max(2, int(round(dt / 0.0005))); h = dt / nsub
    for k in range(int(round(T_max / dt))):
        act = st == 0
        if not act.any(): break
        u = np.zeros((2, n)); u[:, act] = C.uc(X[:, act])
        for _ in range(nsub):
            a1 = f_np(X, u); a2 = f_np(X + h / 2 * a1, u); a3 = f_np(X + h / 2 * a2, u); a4 = f_np(X + h * a3, u)
            X = np.where(act, X + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4), X); y = fk(X[:2])
            hit = act & (st == 0) & (phi_np(*y) <= 0); st[hit] = 1; tt[hit] = (k + 1) * dt
            st[act & (st == 0) & (np.broadcast_to(C.psi(*y), (n,)) < 0)] = 2; act = st == 0
        Pn = C.Psi(X); dmin = np.where(act, np.minimum(dmin, Pn - P), dmin); P = Pn
    return st, tt, dmin


def mpc_run(x0, terminal="cert", dt=T, N=N, omega=None, T_max=T_MAX_MPC):
    return ra.run_one(np.asarray(x0, float), K1_FILE, MU_BASE, SCALE, dt, N, DFRAC, omega, T_max, terminal, U_MAX, SAFE_MARGIN, max(1, int(round(dt / 0.002))))


def safety_until_target(x_samples, u_applied, dt, nsub=8):
    m = np.inf
    for k in range(len(u_applied)):
        x = np.asarray(x_samples[k], float).copy()
        for _ in range(nsub):
            x = ra.rk4_np(x, u_applied[k], dt / nsub); y = fk(x[:2])
            if phi_np(y[0], y[1]) <= 0: return m, True
            m = min(m, float(C.psi(y[0], y[1])))
    return m, False
