"""Algorithm 1 of the paper for the Dubins car, with the certificate of cert.Cert, and the vanilla MPC.
  * terminal set {Psi >= delta} and cost sum_i l(x_i, u_i) + omega (M_Psi - Psi(x_N))
  * warm start = shifted candidate {u*_1..u*_{N-1}, u^c(x*_N)}
  * applied sequence = cheaper of {solver output, feasible candidate}
  * adaptive horizon: N <- max(N-1, 1) if the candidate rollout touches X_T
  * logs J*_k, candidate cost and horizon, and checks Corollary 1.
"""
import numpy as np
import sympy as sp
import casadi as ca

# ── problem data ─────────────────────────────────────────────────────────────
U_MAX = np.array([3.0, 3.0])
TARGET_POS = np.array([-1.7, 0.0])
Q_Y, R_U = 5.0, 0.05
ELL_LOWER = Q_Y * 0.1  # ||h(x)-y*||^2 > 0.1 outside the target ellipse

x1, x2, th, v = sp.symbols("x1 x2 th v")
_X = [x1, x2, th, v]
_y1, _y2 = x1, x2
_safe_sym = 1e-3 * (-((_y2) ** 2 + (2 * (_y1 + 1.7)) ** 2 - 0.4) + 300) * (
    -(_y1**4 + _y2**4 - 16) * (_y1**4 + _y2**4 - 4))
_tgt_sym = _y2**2 + (2 * (_y1 + 1.7)) ** 2 - 0.4

safe_np = sp.lambdify(_X, _safe_sym, "numpy")
tgt_np = sp.lambdify(_X, _tgt_sym, "numpy")


def _s2c(e, m):
    if isinstance(e, sp.Symbol):
        return m[e]
    if e.is_Number:
        return float(e)
    if isinstance(e, sp.Add):
        return sum(_s2c(a, m) for a in e.args)
    if isinstance(e, sp.Mul):
        r = 1.0
        for a in e.args:
            r = r * _s2c(a, m)
        return r
    if isinstance(e, sp.Pow):
        return _s2c(e.args[0], m) ** _s2c(e.args[1], m)
    if isinstance(e, sp.sin):
        return ca.sin(_s2c(e.args[0], m))
    if isinstance(e, sp.cos):
        return ca.cos(_s2c(e.args[0], m))
    raise ValueError(type(e))


_xs = ca.SX.sym("x", 4)
_m = dict(zip(_X, [_xs[i] for i in range(4)]))
safe_ca = ca.Function("safe", [_xs], [_s2c(sp.expand(_safe_sym), _m)])


def rhs_np(x, u):
    return np.array([x[3] * np.cos(x[2]), x[3] * np.sin(x[2]), u[0], u[1]])


def rk4_np(x, u, dt):
    k1 = rhs_np(x, u); k2 = rhs_np(x + dt / 2 * k1, u)
    k3 = rhs_np(x + dt / 2 * k2, u); k4 = rhs_np(x + dt * k3, u)
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def _rk4_ca(dt):
    xc, uc = ca.SX.sym("x", 4), ca.SX.sym("u", 2)
    f = ca.Function("f", [xc, uc], [ca.vertcat(xc[3] * ca.cos(xc[2]), xc[3] * ca.sin(xc[2]), uc[0], uc[1])])
    k1 = f(xc, uc); k2 = f(xc + dt / 2 * k1, uc); k3 = f(xc + dt / 2 * k2, uc); k4 = f(xc + dt * k3, uc)
    return ca.Function("F", [xc, uc], [xc + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)])


def stage_cost(x, u):
    dy = np.asarray(x[:2]) - TARGET_POS
    return Q_Y * float(dy @ dy) + R_U * float(np.asarray(u) @ np.asarray(u))


class DubinsRAMPC:
    def __init__(self, dt=0.05, N=25, delta=0.1, omega=100.0, M_psi=11.63, safe_margin=1e-3):
        # safe_margin: small tightening of the safety constraint that absorbs the RK4
        # discretisation error w.r.t. the continuous plant
        self.dt, self.N0, self.delta, self.omega, self.M_psi = dt, N, delta, omega, M_psi
        self.safe_margin = safe_margin
        self.F = _rk4_ca(dt)
        self.opts = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 500,
                     "ipopt.tol": 1e-4, "ipopt.acceptable_tol": 1e-3}

    def plant_step(self, x, u, plant):
        """nominal: the MPC's own RK4 step; zoh: continuous plant, fine RK4 sub-steps."""
        if plant == "nominal":
            return rk4_np(x, u, self.dt)
        nsub = max(1, int(round(self.dt / 0.005)))
        for _ in range(nsub):
            x = rk4_np(x, u, self.dt / nsub)
        return x

    def run(self, x0, T_max=10.0, keep_cheaper=True, plant="nominal"):
        """Closed loop of Algorithm 1."""
        import time
        x = np.asarray(x0, float); N = self.N0; U_warm = None
        log = dict(t=[0.0], x=[x.copy()], u=[], J=[], J_cand=[], N=[], ell=[], solve_ms=[],
                   uc_in_bounds=[], used_candidate=[], status="timeout")
        for k in range(int(round(T_max / self.dt))):
            if float(tgt_np(*x)) <= 0:
                log["status"] = "reached"; break
            if float(safe_np(*x)) < -1e-6:
                log["status"] = "unsafe"; break
            _t0 = time.perf_counter()
            ok, Us, Xs = self.solve(x, N, U_warm)
            log["solve_ms"].append(1e3 * (time.perf_counter() - _t0))
            if not ok:
                if U_warm is not None and U_warm.shape[1] >= N:
                    Xw = self.rollout(x, U_warm[:, :N])
                    if self.feasible(Xw, U_warm[:, :N]):
                        ok, Us, Xs = True, U_warm[:, :N], Xw
                if not ok:
                    log["status"] = "infeasible"; break
            Js = self.cost(Xs, Us); used = False
            if keep_cheaper and U_warm is not None and U_warm.shape[1] >= N:
                Xw = self.rollout(x, U_warm[:, :N])
                if self.feasible(Xw, U_warm[:, :N]):
                    Jw = self.cost(Xw, U_warm[:, :N])
                    if Jw < Js:
                        Us, Xs, Js, used = U_warm[:, :N], Xw, Jw, True
            log["J"].append(Js); log["N"].append(N); log["used_candidate"].append(used)
            log["u"].append(Us[:, 0].copy()); log["ell"].append(stage_cost(x, Us[:, 0]))
            Uc, Xc, touch, inb = self.candidate(Us, Xs)
            log["uc_in_bounds"].append(inb)
            if touch:
                N = max(N - 1, 1)
            U_warm = Uc
            x = self.plant_step(x, Us[:, 0], plant)
            if not np.all(np.isfinite(x)):
                log["status"] = "nonfinite"; break
            Xw = self.rollout(x, Uc[:, :N])
            log["J_cand"].append(self.cost(Xw, Uc[:, :N]) if self.feasible(Xw, Uc[:, :N]) else np.nan)
            log["t"].append(log["t"][-1] + self.dt); log["x"].append(x.copy())
        else:
            if float(tgt_np(*x)) <= 0:
                log["status"] = "reached"
        for key in ("t", "J", "J_cand", "N", "ell"):
            log[key] = np.array(log[key])
        log["x"] = np.array(log["x"]); log["u"] = np.array(log["u"]).reshape(-1, 2)
        return log


def corollary_check(log, dt, ell_lower=ELL_LOWER, tol=1e-6):
    """Decrease J_{k+1} <= J_k - ell_lower and entry time <= (floor(J0/ell)+1) dt."""
    J = log["J"]
    dec = np.diff(J) + ell_lower if len(J) > 1 else np.array([])
    t_bound = (np.floor(J[0] / ell_lower) + 1) * dt if len(J) else np.nan
    t_entry = log["t"][-1] if log["status"] == "reached" else np.nan
    return dict(n_steps=len(J), n_decrease_viol=int((dec > tol).sum()),
                max_decrease_viol=float(dec.max()) if len(dec) else 0.0,
                t_entry=t_entry, t_bound=t_bound,
                bound_ok=bool(t_entry <= t_bound + 1e-9) if log["status"] == "reached" else False)


def terminal_stats(cert, delta, dt, n_y=120000):
    """Sampled one-step increase of Psi and sup of the stage cost on {Psi >= delta} \\ X_T."""
    x, _ = cert.sample(delta, n_y=n_y, seed=1)
    u = cert.uc(x); xn = cert.zoh(x, u, dt)
    eps = float((cert.Psi(xn) - cert.Psi(x)).min())
    ell = Q_Y * ((x[0] - TARGET_POS[0])**2 + (x[1] - TARGET_POS[1])**2) + R_U * (u**2).sum(0)
    return eps, float(ell.max()), float(np.abs(u).max())


class CertMPC(DubinsRAMPC):
    def __init__(self, cert, dt=0.05, N=25, dfrac=0.3, omega=None, safe_margin=1e-3, u_max=None, nsub=None, safe_sub=False):
        self.c = cert
        if nsub is not None:
            self.NSUB = int(nsub)
        # safe_sub: enforce the safety constraint at the RK4 sub-step nodes as well (inter-sample safety)
        self.safe_sub = safe_sub
        self.u_max = np.asarray(U_MAX if u_max is None else u_max, float)
        M_psi = 1.0 - cert.shift
        delta = dfrac * M_psi
        self.eps, self.sup_ell, self.uc_max = terminal_stats(cert, delta, dt)
        if omega is None:
            omega = 1.2 * self.sup_ell / self.eps if self.eps > 0 else 1e4
        super().__init__(dt=dt, N=N, delta=delta, omega=omega, M_psi=M_psi, safe_margin=safe_margin)
        expr, X = cert.psi_state_sym()
        xs = ca.SX.sym("x", 4)
        self.psi_ca = ca.Function("psi", [xs], [_s2c(sp.expand(expr), dict(zip(X, [xs[i] for i in range(4)])))])
        self.cscale = 1.0 / max(1.0, omega)
        self._build_model()
        # tight tolerances: with a large omega the cost is sensitive to dynamics-constraint residuals
        self.opts = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 300, "ipopt.tol": 1e-5,
                     "ipopt.constr_viol_tol": 1e-8, "ipopt.acceptable_tol": 1e-4, "ipopt.acceptable_constr_viol_tol": 1e-7,
                     "ipopt.acceptable_iter": 5}

    NSUB = 4   # RK4 sub-steps per sampling interval in the MPC model (removes the integrator mismatch)

    def _build_model(self):
        f1 = _rk4_ca(self.dt / self.NSUB)
        xc, uc = ca.SX.sym("x", 4), ca.SX.sym("u", 2); xn = xc; mids = []
        for i in range(self.NSUB):
            xn = f1(xn, uc)
            if i < self.NSUB - 1:
                mids.append(xn)
        self.F = ca.Function("F", [xc, uc], [xn])
        self.Fmid = ca.Function("Fmid", [xc, uc], [ca.horzcat(*mids)]) if mids else None

    def _mids(self, x, u):
        out = []
        for _ in range(self.NSUB - 1):
            x = rk4_np(x, u, self.dt / self.NSUB); out.append(x)
        return np.array(out).T

    def _step(self, x, u):
        for _ in range(self.NSUB):
            x = rk4_np(x, u, self.dt / self.NSUB)
        return x

    def rollout(self, x0, U):
        X = [np.asarray(x0, float)]
        for k in range(U.shape[1]):
            X.append(self._step(X[-1], U[:, k]))
        return np.array(X).T

    def plant_step(self, x, u, plant):
        if plant == "nominal":
            return self._step(x, u)
        return super().plant_step(x, u, plant)

    def _psi(self, X):
        X = np.atleast_2d(np.asarray(X, float)); X = X.reshape(4, -1)
        return self.c.Psi(X)

    def cost(self, X, U):
        J = sum(stage_cost(X[:, k], U[:, k]) for k in range(U.shape[1]))
        return J + self.omega * (self.M_psi - float(self._psi(X[:, -1])[0]))

    def feasible(self, X, U, tol=1e-6):
        ok = (np.all(np.abs(U) <= self.u_max[:, None] + tol)
              and np.all(safe_np(*X[:, 1:]) >= self.safe_margin - tol)
              and float(self._psi(X[:, -1])[0]) >= self.delta - tol)
        if ok and self.safe_sub and self.NSUB > 1:
            for k in range(U.shape[1]):
                m = self.safe_margin if k > 0 else min(self.safe_margin, float(safe_np(*X[:, 0])))
                if np.any(safe_np(*self._mids(X[:, k], U[:, k])) < m - tol):
                    return False
        return bool(ok)

    def solve(self, x0, N, U_warm=None):
        opti = ca.Opti()
        X = opti.variable(4, N + 1); U = opti.variable(2, N)
        x0 = np.asarray(x0, float)
        opti.subject_to(X[:, 0] == x0)
        J = 0; tgt = ca.DM(TARGET_POS)
        for k in range(N):
            opti.subject_to(X[:, k + 1] == self.F(X[:, k], U[:, k]))
            if not self.safe_sub:
                opti.subject_to(safe_ca(X[:, k]) >= self.safe_margin)
            else:
                # the current state is given; later samples and all sub-step nodes are constrained
                if k > 0:
                    opti.subject_to(safe_ca(X[:, k]) >= self.safe_margin)
                if self.Fmid is not None:
                    Xm = self.Fmid(X[:, k], U[:, k]); m = self.safe_margin if k > 0 else min(self.safe_margin, float(safe_np(*x0)))
                    for j in range(self.NSUB - 1):
                        opti.subject_to(safe_ca(Xm[:, j]) >= m)
            dy = X[:2, k] - tgt
            J += Q_Y * ca.dot(dy, dy) + R_U * ca.dot(U[:, k], U[:, k])
        opti.subject_to(safe_ca(X[:, N]) >= self.safe_margin)
        opti.subject_to(self.psi_ca(X[:, N]) >= self.delta)
        opti.subject_to(opti.bounded(-self.u_max[0], U[0, :], self.u_max[0]))
        opti.subject_to(opti.bounded(-self.u_max[1], U[1, :], self.u_max[1]))
        J += self.omega * (self.M_psi - self.psi_ca(X[:, N]))
        opti.minimize(self.cscale * J)
        if U_warm is not None:
            Uw = np.asarray(U_warm)[:, :N]
            opti.set_initial(U, Uw); opti.set_initial(X, self.rollout(x0, Uw))
        else:
            tgt_full = np.array([*TARGET_POS, x0[2], 0.2])
            opti.set_initial(X, np.array([(1 - k / N) * x0 + k / N * tgt_full for k in range(N + 1)]).T)
            opti.set_initial(U, 0)
        opti.solver("ipopt", self.opts)
        try:
            sol = opti.solve()
            Us = np.atleast_2d(sol.value(U)).reshape(2, N)
            return True, Us, self.rollout(x0, Us)      # true nominal trajectory, not the solver iterate
        except Exception:
            if U_warm is None and not getattr(self, "_retry", False):
                self._retry = True
                try:
                    for Ug in self._guesses(x0, N):
                        ok, Us, Xs = self.solve(x0, N, U_warm=Ug)
                        if ok:
                            return ok, Us, Xs
                finally:
                    self._retry = False
            return False, None, None

    def _guesses(self, x0, N):
        """Initial guesses tried at k = 0 after the cold start fails: the terminal-target (vanilla) solution, a
        steer-to-target rollout, and a continuation on the time scale of the terminal set (slower sets first)."""
        ok, Uv = vanilla_solve(self, x0, N)
        if ok:
            yield Uv
        yield steer_guess(x0, N, self.dt, self.u_max, step=self._step)
        if getattr(self, "_no_continuation", False):
            return
        Uw = None
        for f in (0.25, 0.5):
            m = CertMPC(self.c.scaled(f), dt=self.dt, N=N, dfrac=self.delta / self.M_psi, omega=self.omega, safe_margin=self.safe_margin,
                        u_max=self.u_max, nsub=self.NSUB, safe_sub=self.safe_sub)
            m._no_continuation = True
            ok, Uw, _ = m.solve(x0, N, U_warm=Uw)
            if not ok:
                return
        yield Uw

    def candidate(self, U_star, X_star):
        xN = X_star[:, -1].reshape(4, 1)
        uN = self.c.uc(xN)[:, 0] if xN[3, 0] > 1e-6 else np.zeros(2)
        Uc = np.hstack([U_star[:, 1:], uN[:, None]])
        Xc = self.rollout(X_star[:, 1], Uc)
        touch = np.any(tgt_np(*Xc) <= 0)
        xs = Xc[:, -2].copy()
        for _ in range(10):
            xs = rk4_np(xs, uN, self.dt / 10)
            touch = touch or float(tgt_np(*xs)) <= 0
        return Uc, Xc, bool(touch), bool(np.all(np.abs(uN) <= self.u_max))


def steer_guess(x0, N, dt, u_max, v_cruise=1.0, step=None):
    """Open-loop guess: turn the heading towards the target at the rate limit and track a cruise speed."""
    x = np.asarray(x0, float).copy(); U = np.zeros((2, N)); step = step or (lambda x, u: rk4_np(x, u, dt))
    for k in range(N):
        d = TARGET_POS - x[:2]; th_d = np.arctan2(d[1], d[0]); e = (th_d - x[2] + np.pi) % (2 * np.pi) - np.pi
        if x[3] < 0:
            e = (e + np.pi + np.pi) % (2 * np.pi) - np.pi
        U[:, k] = [np.clip(4.0 * e, -u_max[0], u_max[0]), np.clip(3.0 * (v_cruise - x[3]), -u_max[1], u_max[1])]
        x = step(x, U[:, k])
    return U


def vanilla_solve(mpc, x0, N, Qf=80.0):
    """Terminal-target MPC of the paper's baseline (same dynamics, cost and safety constraint); fallback initial guess of RA-MPC."""
    um = getattr(mpc, "u_max", U_MAX)
    opti = ca.Opti(); X = opti.variable(4, N + 1); U = opti.variable(2, N)
    x0 = np.asarray(x0, float); opti.subject_to(X[:, 0] == x0); J = 0; tgt = ca.DM(TARGET_POS)
    for k in range(N):
        opti.subject_to(X[:, k + 1] == mpc.F(X[:, k], U[:, k])); opti.subject_to(safe_ca(X[:, k]) >= 0)
        dy = X[:2, k] - tgt; J += Q_Y * ca.dot(dy, dy) + R_U * ca.dot(U[:, k], U[:, k])
    opti.subject_to(safe_ca(X[:, N]) >= 0)
    opti.subject_to(X[1, N]**2 + (2 * (X[0, N] + 1.7))**2 - 0.4 <= 0)
    opti.subject_to(opti.bounded(-um[0], U[0, :], um[0])); opti.subject_to(opti.bounded(-um[1], U[1, :], um[1]))
    dyN = X[:2, N] - tgt; J += Qf * ca.dot(dyN, dyN); opti.minimize(J)
    tgt_full = np.array([*TARGET_POS, x0[2], 0.0])
    opti.set_initial(X, np.array([(1 - k / N) * x0 + k / N * tgt_full for k in range(N + 1)]).T); opti.set_initial(U, 0)
    opti.solver("ipopt", mpc.opts)
    try:
        sol = opti.solve(); return True, np.atleast_2d(sol.value(U)).reshape(2, N)
    except Exception:
        return False, None


def run_one(x0, k1_file, mu, scale, dt, N, dfrac, omega, T_max, plant, u_max=None, nsub=None, safe_sub=False, safe_margin=1e-3):
    from cert import Cert
    c = Cert(k1_file, mu=mu)
    if scale != 1.0:
        c = c.scaled(scale)
    mpc = CertMPC(c, dt=dt, N=N, dfrac=dfrac, omega=omega, u_max=u_max, nsub=nsub, safe_sub=safe_sub, safe_margin=safe_margin)
    log = mpc.run(x0, T_max=T_max, plant=plant)
    log["uc_max_terminal"] = mpc.uc_max
    log["check"] = corollary_check(log, dt)
    log["omega"] = mpc.omega; log["eps"] = mpc.eps
    return log
