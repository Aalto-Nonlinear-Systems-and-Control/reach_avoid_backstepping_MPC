"""Final Dubins setting for the TAC paper: certificate R2, one initial set for all methods."""
import os, sys
import numpy as np
import casadi as ca

HERE = os.path.dirname(os.path.abspath(__file__))
DESIGN = os.path.abspath(os.path.join(HERE, "..", "design"))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for _p in (DESIGN, os.path.join(ROOT, "python")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from cert import Cert                      # noqa: E402
import ra_mpc_dubins as rm                 # noqa: E402
import ra_mpc_cert as R                    # noqa: E402

# ── certificate (REPORT.md Section 11): psi = alpha*h/661.409, k1 of design/k1/W3_l0.06.txt, time scale s ──
K1_FILE = os.path.join(DESIGN, "k1", "W3_l0.06.txt")
MU_BASE, DFRAC = 0.02, 0.2
# time scale of the terminal controller: SCALE env, else data/scale.json (solve_scale.py), else the earlier 0.159
import json as _json
_sj = os.path.join(HERE, "data", "scale.json")
SCALE = float(os.environ.get("SCALE", _json.load(open(_sj))["s"] if os.path.exists(_sj) else 0.159))
# URA baseline (the synthesis without regard to the bounds or the period) at time scale URA_SCALE; the initial velocities
# are drawn around the field of scale SET_SCALE with the radius scaled accordingly. Defaults: the paper setting (s = 5).
# Environment overrides (SET_SCALE, URA_SCALE, SCALE, TAG) allow a variant to be run into data<TAG>/ and figs<TAG>/;
# the earlier variants (terminal controller at s = 0.159, data_s0159/; s = 5, data_s5/) are on branch tac-revision-archive
# and in data_s5/ (terminal 0.159, set and URA at 5).
SET_SCALE = float(os.environ.get("SET_SCALE", 5.0))
URA_SCALE = float(os.environ.get("URA_SCALE", 5.0))
TAG = os.environ.get("TAG", "")
# ── MPC ──
T, N, U_MAX = 0.05, 25, (3.0, 3.0)
T_HOR = T * N
# model and path constraint: RK4 sub-steps of H_SUB seconds; the safety constraint is imposed at every sub-step node
# with the back-off SAFE_MARGIN (in units of rm.safe_ca, about 1.2 mm at the inner boundary) for inter-sample safety
H_SUB, SAFE_MARGIN = 0.0125, 0.05
MPC_OPTS = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 300, "ipopt.tol": 1e-5, "ipopt.constr_viol_tol": 1e-8,
            "ipopt.acceptable_tol": 1e-4, "ipopt.acceptable_constr_viol_tol": 1e-7, "ipopt.acceptable_iter": 5}
T_MAX_MPC, T_MAX_URA = 20.0, 600.0
# ── initial set: positions uniform in X_S minus X_T, velocity uniform in the disc of radius R0 around k1(y) ──
R0_S1, SEED_POS, SEED_VEL, N_INIT = 0.03 / 0.159, 42, 2026, 300
R0 = R0_S1 * SET_SCALE                      # the largest velocity radius of X_RA at the set's scale (0.189 at s = 1)
DATA, FIGS = os.path.join(HERE, "data" + TAG), os.path.join(HERE, "figs" + TAG)
os.makedirs(DATA, exist_ok=True); os.makedirs(FIGS, exist_ok=True)


def cert():
    return Cert(K1_FILE, mu=MU_BASE).scaled(SCALE)


C = cert()
M_PSI = 1.0 - C.shift
DELTA = DFRAC * M_PSI


def positions():
    np.random.seed(SEED_POS)
    lo = np.array([-2, -2, 2 * np.pi / 3, -1.0]); hi = np.array([2, 2, 4 * np.pi / 3, 1.0])
    S = np.random.rand(4, 1000) * (hi - lo)[:, None] + lo[:, None]
    cand = np.where((rm.safe_np(*S) >= 0) & (rm.tgt_np(*S) > 0))[0][:N_INIT]
    return S[:2, cand]


def initial_set(r0=R0, set_scale=SET_SCALE):
    P = positions(); n = P.shape[1]; r = np.random.default_rng(SEED_VEL)
    w = np.sqrt(r.random(n)) * np.exp(1j * r.uniform(0, 2 * np.pi, n))
    cs = C if abs(set_scale - SCALE) < 1e-12 else Cert(K1_FILE, mu=MU_BASE).scaled(set_scale)
    k = np.array([np.broadcast_to(np.asarray(a, float), (n,)) for a in cs.k1(*P)])
    yd = k + r0 * np.array([w.real, w.imag])
    return np.vstack([P, np.arctan2(yd[1], yd[0]), np.hypot(yd[0], yd[1])])


def f_np(x, u):
    return np.array([x[3] * np.cos(x[2]), x[3] * np.sin(x[2]), u[0] + 0 * x[0], u[1] + 0 * x[0]])


def zoh_np(x, u, T, h=0.005):
    """Continuous plant under a frozen input, fine RK4 sub-steps (works on (4,) and (4, n) arrays)."""
    n = max(1, int(round(T / h))); h = T / n
    for _ in range(n):
        a1 = f_np(x, u); a2 = f_np(x + h / 2 * a1, u); a3 = f_np(x + h / 2 * a2, u); a4 = f_np(x + h * a3, u)
        x = x + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
    return x


# ── vanilla MPC: same model, stage cost, path constraints and solver options as RA-MPC; terminal constraint = target set ──
def nsub_of(dt):
    return max(1, int(round(dt / H_SUB)))


_MODELS = {}


def model(dt):
    """F: one sampling period (nsub RK4 sub-steps); Fmid: the states at the interior sub-step nodes."""
    ns = nsub_of(dt)
    if (dt, ns) not in _MODELS:
        f1 = rm._rk4_ca(dt / ns); xc, uc = ca.SX.sym("x", 4), ca.SX.sym("u", 2); xn = xc; mids = []
        for i in range(ns):
            xn = f1(xn, uc)
            if i < ns - 1:
                mids.append(xn)
        _MODELS[(dt, ns)] = (ca.Function("F", [xc, uc], [xn]), ca.Function("Fmid", [xc, uc], [ca.horzcat(*mids)]) if mids else None, ns)
    return _MODELS[(dt, ns)]


def step_np(x, u, dt):
    ns = nsub_of(dt)
    for _ in range(ns):
        x = rm.rk4_np(x, u, dt / ns)
    return x


def vanilla_solve(x0, dt=T, N=N, u_max=U_MAX, U_warm=None, Qf=80.0, _retry=False, _shift=True):
    F, Fmid, ns = model(dt); opti = ca.Opti(); X = opti.variable(4, N + 1); U = opti.variable(2, N)
    x0 = np.asarray(x0, float); opti.subject_to(X[:, 0] == x0); J = 0; tgt = ca.DM(rm.TARGET_POS)
    for k in range(N):
        opti.subject_to(X[:, k + 1] == F(X[:, k], U[:, k]))
        if k > 0:
            opti.subject_to(rm.safe_ca(X[:, k]) >= SAFE_MARGIN)
        if Fmid is not None:
            Xm = Fmid(X[:, k], U[:, k]); m = SAFE_MARGIN if k > 0 else min(SAFE_MARGIN, float(rm.safe_np(*x0)))
            for j in range(ns - 1):
                opti.subject_to(rm.safe_ca(Xm[:, j]) >= m)
        dy = X[:2, k] - tgt; J += rm.Q_Y * ca.dot(dy, dy) + rm.R_U * ca.dot(U[:, k], U[:, k])
    opti.subject_to(X[1, N]**2 + (2 * (X[0, N] + 1.7))**2 - 0.4 <= 0)
    opti.subject_to(opti.bounded(-u_max[0], U[0, :], u_max[0])); opti.subject_to(opti.bounded(-u_max[1], U[1, :], u_max[1]))
    dyN = X[:2, N] - tgt; J += Qf * ca.dot(dyN, dyN); opti.minimize(J)
    if U_warm is not None:
        Uw = np.hstack([U_warm[:, 1:], U_warm[:, -1:]]) if _shift else np.asarray(U_warm); Xw = np.zeros((4, N + 1)); Xw[:, 0] = x0
        for k in range(N):
            Xw[:, k + 1] = step_np(Xw[:, k], Uw[:, k], dt)
        opti.set_initial(U, Uw); opti.set_initial(X, Xw)
    else:
        tgt_full = np.array([*rm.TARGET_POS, x0[2], 0.0])
        opti.set_initial(X, np.array([(1 - k / N) * x0 + k / N * tgt_full for k in range(N + 1)]).T); opti.set_initial(U, 0)
    opti.solver("ipopt", MPC_OPTS)
    try:
        sol = opti.solve(); return True, np.atleast_2d(sol.value(U)).reshape(2, N)
    except Exception:
        if U_warm is None and not _retry:
            # same fallback as RA-MPC at k = 0: a steer-to-target rollout as initial guess
            return vanilla_solve(x0, dt, N, u_max, R.steer_guess(x0, N, dt, u_max, step=lambda x, u: step_np(x, u, dt)), Qf, _retry=True, _shift=False)
        return False, None


def vanilla_rollout(x0, dt=T, N=N, u_max=U_MAX, T_max=T_MAX_MPC):
    """Receding-horizon vanilla MPC on the continuous plant; returns t, u, x, status, solve times [ms]."""
    import time
    x = np.asarray(x0, float).copy(); ts, us, xs, ms = [0.0], [], [x.copy()], []; Uw = None; status = "timeout"
    for _ in range(int(round(T_max / dt))):
        if float(rm.tgt_np(*x)) <= 0:
            status = "reached"; break
        if float(rm.safe_np(*x)) < -1e-6:
            status = "unsafe"; break
        t0 = time.perf_counter(); ok, U = vanilla_solve(x, dt, N, u_max, Uw); ms.append(1e3 * (time.perf_counter() - t0))
        if not ok:
            status = "infeasible"; break
        Uw = U; x = zoh_np(x, U[:, 0], dt); us.append(U[:, 0].copy()); ts.append(ts[-1] + dt); xs.append(x.copy())
    else:
        if float(rm.tgt_np(*x)) <= 0:
            status = "reached"
    return dict(t=np.array(ts), u=np.array(us).reshape(-1, 2), x=np.array(xs), status=status, solve_ms=np.array(ms))


def ra_run(x0, dt=T, N=N, u_max=U_MAX, omega=None, T_max=T_MAX_MPC):
    """Algorithm 1 (RA-MPC) on the continuous plant with the final certificate."""
    return R.run_one(np.asarray(x0, float), K1_FILE, MU_BASE, SCALE, dt, N, DFRAC, omega, T_max, "zoh", np.asarray(u_max, float),
                     nsub=nsub_of(dt), safe_sub=True, safe_margin=SAFE_MARGIN)


def ring(x):
    """h(y) of the safe set X_S = {h(y) >= 0} (ring between two superellipses)."""
    r = x[0]**4 + x[1]**4
    return -(r - 16) * (r - 4)


def safety_until_target(x_samples, u_applied, dt, nsub=20):
    """Smallest h along the continuous closed loop up to the first entry into X_T."""
    m = np.inf
    for k in range(len(u_applied)):
        x = np.asarray(x_samples[k], float).copy()
        for _ in range(nsub):
            x = zoh_np(x, u_applied[k], dt / nsub, h=min(0.005, dt / nsub))
            if float(rm.tgt_np(*x)) <= 0:
                return m, True
            m = min(m, float(ring(x)))
    return m, False


def ura_run(x0, T_max=T_MAX_URA, dense=False, scale=SCALE):
    """u^c alone in continuous time: status, time, peak inputs, accumulated stage cost, smallest safe-set value."""
    from scipy.integrate import solve_ivp
    c = Cert(K1_FILE, mu=MU_BASE).scaled(scale)
    def rhs(t, z):
        x = z[:4]; u = c.uc(x.reshape(4, 1))[:, 0]; dy = x[:2] - rm.TARGET_POS
        return [x[3] * np.cos(x[2]), x[3] * np.sin(x[2]), u[0], u[1], rm.Q_Y * dy @ dy + rm.R_U * u @ u]
    e1 = lambda t, z: float(rm.tgt_np(*z[:4])); e1.terminal = True; e1.direction = -1
    e2 = lambda t, z: float(rm.safe_np(*z[:4])); e2.terminal = True; e2.direction = -1
    e3 = lambda t, z: abs(z[3]) - 1e-4; e3.terminal = True; e3.direction = -1
    try:
        s = solve_ivp(rhs, (0, T_max), np.append(np.asarray(x0, float), 0.0), events=[e1, e2, e3], rtol=1e-7, atol=1e-9, max_step=0.05)
    except Exception:
        return dict(status="error", t=np.nan, upeak=np.zeros(2), cost=np.nan)
    u = c.uc(s.y[:4])
    st = "reached" if s.t_events[0].size else "unsafe" if s.t_events[1].size else "v=0" if s.t_events[2].size else "timeout"
    out = dict(status=st, t=float(s.t[-1]), upeak=np.abs(u).max(1), cost=float(s.y[4, -1]), psi0=float(c.Psi(s.y[:4, :1])[0]),
               t_viol=float(np.sum(np.diff(s.t) * (np.abs(u[:, 1:]) > np.asarray(U_MAX)[:, None]).any(0))))
    if dense:
        out.update(tt=s.t, x=s.y[:4].T, u=u.T, psi=c.Psi(s.y[:4]))
    return out


def ura_zoh(X0, dt=T, T_max=T_MAX_URA, track=False):
    """Sampled-data u^c (zero-order hold, period dt), vectorised over the columns of X0.
    status: 1 reached, 2 left the safe set, 0 timeout. Also the smallest one-step increase of Psi at the samples."""
    X = X0.copy(); n = X.shape[1]; st = np.zeros(n, int); tt = np.full(n, np.nan)
    nsub = max(4, int(round(dt / 0.005))); h = dt / nsub; dmin = np.full(n, np.inf); P = C.Psi(X); hist = [P.copy()] if track else None
    for k in range(int(round(T_max / dt))):
        act = st == 0
        if not act.any():
            break
        u = np.zeros((2, n)); u[:, act] = C.uc(X[:, act])
        for _ in range(nsub):
            a1 = f_np(X, u); a2 = f_np(X + h / 2 * a1, u); a3 = f_np(X + h / 2 * a2, u); a4 = f_np(X + h * a3, u)
            X = np.where(act, X + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4), X)
            hit = act & (st == 0) & (rm.tgt_np(*X) <= 0); st[hit] = 1; tt[hit] = (k + 1) * dt
            st[act & (st == 0) & (rm.safe_np(*X) < 0)] = 2
            act = st == 0
        Pn = C.Psi(X); dmin = np.where(act, np.minimum(dmin, Pn - P), dmin); P = Pn
        if track:
            hist.append(np.where(act, P, np.nan))
    return (st, tt, dmin, np.array(hist)) if track else (st, tt, dmin)
