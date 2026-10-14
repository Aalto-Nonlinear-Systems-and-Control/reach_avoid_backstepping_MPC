"""Algorithm 1 (RA-MPC) and the vanilla MPC for the top-view arm (copy of ../arm/ra_mpc_arm.py on cert_top), with a certificate of cert_top.CertArm.
Same structure as design/ra_mpc_cert.py for the Dubins car: RK4 sub-step model, state constraint at every sub-step
node with a back-off, terminal set {Psi >= delta}, terminal cost omega (M_Psi - Psi), shifted candidate ending with
u^c, keep-the-cheaper rule, adaptive horizon, fallback initial guesses at k = 0."""
import os, sys, time
import numpy as np, sympy as sp, casadi as ca
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cert_top import (CertArm, fk, ik, phi_np, gravity, f_np, mass, coriolis_qdot, jac, M1, M2, L1, L2, LC1, LC2, I1, I2, GR, TARGET_POS, PHI, RT, y1, y2)

Q_Y, R_U, QF_VAN = 5.0, 0.001, 80.0
ELL_LOWER = Q_Y * RT**2         # ||y - y*||^2 >= RT^2 outside the target disc
X1, X2, X3, X4 = sp.symbols("x1 x2 x3 x4")
H1 = L1 * sp.cos(X1) + L2 * sp.cos(X1 + X2); H2 = L1 * sp.sin(X1) + L2 * sp.sin(X1 + X2)


def stage_cost(x, u):
    dy = fk(np.asarray(x, float)[:2]) - TARGET_POS
    return Q_Y * float(dy @ dy) + R_U * float(np.asarray(u) @ np.asarray(u))


def tgt_np(x):
    y = fk(np.asarray(x)[:2]); return phi_np(y[0], y[1])


def s2c(e, mp):
    """sympy -> casadi (polynomials, sin, cos)."""
    if isinstance(e, sp.Symbol): return mp[e]
    if e.is_Number: return float(e)
    if isinstance(e, sp.Add): return sum(s2c(a, mp) for a in e.args)
    if isinstance(e, sp.Mul):
        r = 1.0
        for a in e.args: r = r * s2c(a, mp)
        return r
    if isinstance(e, sp.Pow): return s2c(e.args[0], mp) ** float(e.args[1])
    if isinstance(e, sp.sin): return ca.sin(s2c(e.args[0], mp))
    if isinstance(e, sp.cos): return ca.cos(s2c(e.args[0], mp))
    raise ValueError(type(e))


def arm_rhs_ca(x, u):
    q1, q2, d1, d2 = x[0], x[1], x[2], x[3]; c2, s2 = ca.cos(q2), ca.sin(q2)
    M11 = I1 + I2 + M1 * LC1**2 + M2 * (L1**2 + LC2**2 + 2 * L1 * LC2 * c2); M12 = M2 * (LC2**2 + L1 * LC2 * c2) + I2; M22 = M2 * LC2**2 + I2
    h = M2 * L1 * LC2 * s2
    Cq = ca.vertcat(-h * d2 * d1 - h * (d1 + d2) * d2, h * d1 * d1)
    G = ca.vertcat((M1 * GR * LC1 + M2 * GR * L1) * ca.cos(q1) + M2 * GR * LC2 * ca.cos(q1 + q2), M2 * GR * LC2 * ca.cos(q1 + q2))
    r = u - Cq - G; det = M11 * M22 - M12 * M12
    qdd = ca.vertcat((M22 * r[0] - M12 * r[1]) / det, (-M12 * r[0] + M11 * r[1]) / det)
    return ca.vertcat(d1, d2, qdd[0], qdd[1])


def rk4_np(x, u, h):
    k1 = f_np(x, u); k2 = f_np(x + h / 2 * k1, u); k3 = f_np(x + h / 2 * k2, u); k4 = f_np(x + h * k3, u)
    return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def steer_guess(x0, N, dt, u_max, step):
    """Open-loop guess: gravity compensation plus a saturated PD towards the target's elbow-down configuration."""
    x = np.asarray(x0, float).copy(); U = np.zeros((2, N)); qt = ik(TARGET_POS.reshape(2, 1))[:, 0]
    for k in range(N):
        e = (qt - x[:2] + np.pi) % (2 * np.pi) - np.pi
        U[:, k] = np.clip(gravity(x[:2]) + 60.0 * e - 20.0 * x[2:], -u_max, u_max); x = step(x, U[:, k])
    return U


def uc_guess(cert, x0, N, u_max, step):
    """Open-loop guess: the terminal controller u^c, clipped to the bounds, rolled out over the horizon."""
    x = np.asarray(x0, float).copy(); U = np.zeros((2, N))
    for k in range(N):
        U[:, k] = np.clip(cert.uc(x.reshape(4, 1))[:, 0], -u_max, u_max); x = step(x, U[:, k])
    return U


class ArmMPC:
    """terminal='cert' (RA-MPC, Algorithm 1) or 'target' (vanilla MPC)."""
    def __init__(self, cert, dt=0.002, N=100, dfrac=0.5, omega=None, safe_margin=0.005, u_max=(120.0, 40.0), nsub=1, terminal="cert", opts=None):
        self.c, self.dt, self.N0, self.terminal = cert, dt, N, terminal
        self.u_max = np.asarray(u_max, float); self.NSUB = int(nsub); self.safe_margin = safe_margin
        self.M_psi = 1.0 - cert.shift; self.delta = dfrac * self.M_psi
        if terminal == "cert":
            self.eps, self.sup_ell, self.uc_max = self.terminal_stats()
            self.omega = omega if omega is not None else (1.2 * self.sup_ell / self.eps if self.eps > 0 else 1e4)
        else:
            self.eps, self.sup_ell, self.uc_max, self.omega = np.nan, np.nan, np.zeros(2), 0.0
        # casadi functions
        xs = ca.SX.sym("x", 4); mp = {X1: xs[0], X2: xs[1], X3: xs[2], X4: xs[3]}
        ys = ca.SX.sym("y", 2); my = {y1: ys[0], y2: ys[1]}            # psi and k1 composed with the kinematics inside casadi (no symbolic expansion)
        psi_f = ca.Function("psi_y", [ys], [s2c(sp.expand(self.c.psi_sym), my)]); k1_f = ca.Function("k1_y", [ys], [ca.vertcat(*[s2c(sp.expand(k), my) for k in self.c.k1_sym])])
        hx = ca.vertcat(s2c(H1, mp), s2c(H2, mp)); q1, q2, d1, d2 = xs[0], xs[1], xs[2], xs[3]
        yd = ca.vertcat((-L1 * ca.sin(q1) - L2 * ca.sin(q1 + q2)) * d1 - L2 * ca.sin(q1 + q2) * d2, (L1 * ca.cos(q1) + L2 * ca.cos(q1 + q2)) * d1 + L2 * ca.cos(q1 + q2) * d2)
        self.h_ca = ca.Function("h", [xs], [hx]); self.safe_ca = ca.Function("safe", [xs], [psi_f(hx)])
        self.tgt_ca = ca.Function("tgt", [xs], [(hx[0] - TARGET_POS[0])**2 + (hx[1] - TARGET_POS[1])**2 - RT**2])
        ev = yd - k1_f(hx); self.Psi_ca = ca.Function("Psi", [xs], [psi_f(hx) - self.c.shift - ca.dot(ev, ev) / (2 * self.c.mu)])
        uc_ = ca.SX.sym("u", 2); h = dt / self.NSUB
        def rk4(x, u):
            k1 = arm_rhs_ca(x, u); k2 = arm_rhs_ca(x + h / 2 * k1, u); k3 = arm_rhs_ca(x + h / 2 * k2, u); k4 = arm_rhs_ca(x + h * k3, u)
            return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        xn = xs; mids = []
        for i in range(self.NSUB):
            xn = rk4(xn, uc_)
            if i < self.NSUB - 1: mids.append(xn)
        self.F = ca.Function("F", [xs, uc_], [xn]); self.Fmid = ca.Function("Fmid", [xs, uc_], [ca.horzcat(*mids)]) if mids else None
        self.cscale = 1.0 / max(1.0, self.omega)
        self.opts = opts or {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 300, "ipopt.tol": 1e-5, "ipopt.constr_viol_tol": 1e-8,
                             "ipopt.acceptable_tol": 1e-4, "ipopt.acceptable_constr_viol_tol": 1e-7, "ipopt.acceptable_iter": 5}

    def terminal_stats(self, n_y=120000):
        x, _ = self.c.sample(self.delta, n_y=n_y, seed=1); u = self.c.uc(x); xn = self.c.zoh(x, u, self.dt, n=max(10, int(self.dt / 0.0002)))
        keep = phi_np(*fk(xn[:2])) > 0; eps = float((self.c.Psi(xn) - self.c.Psi(x))[keep].min())
        y = fk(x[:2]); ell = Q_Y * ((y[0] - TARGET_POS[0])**2 + (y[1] - TARGET_POS[1])**2) + R_U * (u**2).sum(0)
        return eps, float(ell.max()), np.abs(u).max(1)

    # ── numpy model ──
    def _step(self, x, u):
        for _ in range(self.NSUB): x = rk4_np(x, u, self.dt / self.NSUB)
        return x

    def rollout(self, x0, U):
        X = [np.asarray(x0, float)]
        for k in range(U.shape[1]): X.append(self._step(X[-1], U[:, k]))
        return np.array(X).T

    def _psi(self, x):
        return float(self.c.Psi(np.asarray(x, float).reshape(4, 1))[0])

    def safe_np(self, x):
        y = fk(np.asarray(x)[:2]); return float(self.c.psi(y[0], y[1]))

    def cost(self, X, U):
        J = sum(stage_cost(X[:, k], U[:, k]) for k in range(U.shape[1]))
        if self.terminal == "cert": return J + self.omega * (self.M_psi - self._psi(X[:, -1]))
        dy = fk(X[:2, -1]) - TARGET_POS; return J + QF_VAN * float(dy @ dy)

    def feasible(self, X, U, tol=1e-6):
        if not np.all(np.abs(U) <= self.u_max[:, None] + tol): return False
        for k in range(1, X.shape[1]):
            if self.safe_np(X[:, k]) < self.safe_margin - tol: return False
        if self.terminal == "cert": return self._psi(X[:, -1]) >= self.delta - tol
        return tgt_np(X[:, -1]) <= tol

    # ── OCP: one parametric NLP per horizon length, built once and cached ──
    def _nlp(self, N):
        if not hasattr(self, "_solvers"): self._solvers = {}
        if N in self._solvers: return self._solvers[N]
        opti = ca.Opti(); X = opti.variable(4, N + 1); U = opti.variable(2, N); x0p = opti.parameter(4); m0p = opti.parameter()
        opti.subject_to(X[:, 0] == x0p); J = 0; tgt = ca.DM(TARGET_POS)
        for k in range(N):
            opti.subject_to(X[:, k + 1] == self.F(X[:, k], U[:, k]))
            if k > 0: opti.subject_to(self.safe_ca(X[:, k]) >= self.safe_margin)
            if self.Fmid is not None:
                Xm = self.Fmid(X[:, k], U[:, k])
                for j in range(self.NSUB - 1): opti.subject_to(self.safe_ca(Xm[:, j]) >= (self.safe_margin if k > 0 else m0p))
            dy = self.h_ca(X[:, k]) - tgt; J += Q_Y * ca.dot(dy, dy) + R_U * ca.dot(U[:, k], U[:, k])
        opti.subject_to(self.safe_ca(X[:, N]) >= self.safe_margin)
        if self.terminal == "cert":
            opti.subject_to(self.Psi_ca(X[:, N]) >= self.delta); J += self.omega * (self.M_psi - self.Psi_ca(X[:, N]))
        else:
            opti.subject_to(self.tgt_ca(X[:, N]) <= 0); dyN = self.h_ca(X[:, N]) - tgt; J += QF_VAN * ca.dot(dyN, dyN)
        opti.subject_to(opti.bounded(-self.u_max[0], U[0, :], self.u_max[0])); opti.subject_to(opti.bounded(-self.u_max[1], U[1, :], self.u_max[1]))
        opti.minimize(self.cscale * J); opti.solver("ipopt", self.opts)
        self._solvers[N] = (opti, X, U, x0p, m0p); return self._solvers[N]

    def solve(self, x0, N, U_warm=None):
        opti, X, U, x0p, m0p = self._nlp(N); x0 = np.asarray(x0, float)
        opti.set_value(x0p, x0); opti.set_value(m0p, min(self.safe_margin, self.safe_np(x0)))
        if U_warm is not None:
            Uw = np.asarray(U_warm)[:, :N]; opti.set_initial(U, Uw); opti.set_initial(X, self.rollout(x0, Uw))
        else:
            qt = ik(TARGET_POS.reshape(2, 1))[:, 0]; xt = np.array([qt[0], qt[1], 0.0, 0.0])
            Xg = np.array([(1 - k / N) * x0 + k / N * xt for k in range(N + 1)]).T; opti.set_initial(X, Xg)
            opti.set_initial(U, np.array([gravity(Xg[:2, k]) for k in range(N)]).T)
        try:
            sol = opti.solve(); Us = np.atleast_2d(sol.value(U)).reshape(2, N); return True, Us, self.rollout(x0, Us)
        except Exception:
            if U_warm is None and not getattr(self, "_retry", False):
                self._retry = True
                try:
                    for Ug in self._guesses(x0, N):
                        ok, Us, Xs = self.solve(x0, N, U_warm=Ug)
                        if ok: return ok, Us, Xs
                finally:
                    self._retry = False
            return False, None, None

    def _guesses(self, x0, N):
        yield uc_guess(self.c, x0, N, self.u_max, self._step)
        yield steer_guess(x0, N, self.dt, self.u_max, self._step)
        if self.terminal == "cert":
            if not hasattr(self, "_van"): self._van = ArmMPC(self.c, self.dt, N, self.delta / self.M_psi, 0.0, self.safe_margin, self.u_max, self.NSUB, "target"); self._van._retry = True
            ok, Uv, _ = self._van.solve(x0, N)
            if ok: yield Uv
            if not getattr(self, "_no_continuation", False):
                Uw = None
                if not hasattr(self, "_cont"):
                    self._cont = {}
                    for f in (0.25, 0.5):
                        self._cont[f] = ArmMPC(self.c.scaled(f), self.dt, N, self.delta / self.M_psi, self.omega, self.safe_margin, self.u_max, self.NSUB, "cert"); self._cont[f]._no_continuation = True
                for f in (0.25, 0.5):
                    ok, Uw, _ = self._cont[f].solve(x0, N, U_warm=Uw)
                    if not ok: return
                yield Uw

    def candidate(self, U_star, X_star):
        xN = X_star[:, -1].reshape(4, 1); uN = self.c.uc(xN)[:, 0]
        Uc = np.hstack([U_star[:, 1:], uN[:, None]]); Xc = self.rollout(X_star[:, 1], Uc)
        touch = any(tgt_np(Xc[:, k]) <= 0 for k in range(Xc.shape[1]))
        xs = Xc[:, -2].copy()
        for _ in range(5):
            xs = rk4_np(xs, uN, self.dt / 5); touch = touch or tgt_np(xs) <= 0
        return Uc, Xc, bool(touch), bool(np.all(np.abs(uN) <= self.u_max))

    def plant_step(self, x, u):
        nsub = max(2, int(round(self.dt / 0.0005)))
        for _ in range(nsub): x = rk4_np(x, u, self.dt / nsub)
        return x

    def run(self, x0, T_max=8.0, keep_cheaper=True):
        """Closed loop. RA-MPC: Algorithm 1; vanilla: receding horizon with a shifted warm start."""
        x = np.asarray(x0, float); N = self.N0; U_warm = None
        log = dict(t=[0.0], x=[x.copy()], u=[], J=[], N=[], ell=[], solve_ms=[], uc_in_bounds=[], used_candidate=[], status="timeout")
        for k in range(int(round(T_max / self.dt))):
            if tgt_np(x) <= 0: log["status"] = "reached"; break
            if self.safe_np(x) < -1e-6: log["status"] = "unsafe"; break
            t0 = time.perf_counter(); ok, Us, Xs = self.solve(x, N, U_warm); log["solve_ms"].append(1e3 * (time.perf_counter() - t0))
            if not ok:
                if U_warm is not None and U_warm.shape[1] >= N:
                    Xw = self.rollout(x, U_warm[:, :N])
                    if self.feasible(Xw, U_warm[:, :N]): ok, Us, Xs = True, U_warm[:, :N], Xw
                if not ok: log["status"] = "infeasible"; break
            Js = self.cost(Xs, Us); used = False
            if self.terminal == "cert" and keep_cheaper and U_warm is not None and U_warm.shape[1] >= N:
                Xw = self.rollout(x, U_warm[:, :N])
                if self.feasible(Xw, U_warm[:, :N]):
                    Jw = self.cost(Xw, U_warm[:, :N])
                    if Jw < Js: Us, Xs, Js, used = U_warm[:, :N], Xw, Jw, True
            log["J"].append(Js); log["N"].append(N); log["used_candidate"].append(used); log["u"].append(Us[:, 0].copy()); log["ell"].append(stage_cost(x, Us[:, 0]))
            if self.terminal == "cert":
                Uc, Xc, touch, inb = self.candidate(Us, Xs); log["uc_in_bounds"].append(inb)
                if touch: N = max(N - 1, 1)
                U_warm = Uc
            else:
                U_warm = np.hstack([Us[:, 1:], Us[:, -1:]])
            x = self.plant_step(x, Us[:, 0])
            if not np.all(np.isfinite(x)): log["status"] = "nonfinite"; break
            log["t"].append(log["t"][-1] + self.dt); log["x"].append(x.copy())
        else:
            if tgt_np(x) <= 0: log["status"] = "reached"
        for key in ("t", "J", "N", "ell"): log[key] = np.array(log[key])
        log["x"] = np.array(log["x"]); log["u"] = np.array(log["u"]).reshape(-1, 2)
        J = log["J"]; log["check"] = dict(n_steps=len(J), t_entry=log["t"][-1] if log["status"] == "reached" else np.nan,
                                         t_bound=(np.floor(J[0] / ELL_LOWER) + 1) * self.dt if len(J) else np.nan,
                                         n_decrease_viol=int(((np.diff(J) + ELL_LOWER) > 1e-6).sum()) if len(J) > 1 else 0)
        log["check"]["bound_ok"] = bool(log["check"]["t_entry"] <= log["check"]["t_bound"] + 1e-9) if log["status"] == "reached" else False
        log["omega"], log["eps"] = self.omega, self.eps
        return log


def run_one(x0, k1_file, mu, scale, dt, N, dfrac, omega, T_max, terminal="cert", u_max=(120.0, 40.0), safe_margin=0.005, nsub=1):
    c = CertArm(k1_file, mu=mu).scaled(scale)
    m = ArmMPC(c, dt=dt, N=N, dfrac=dfrac, omega=omega, safe_margin=safe_margin, u_max=u_max, nsub=nsub, terminal=terminal)
    return m.run(x0, T_max=T_max)
