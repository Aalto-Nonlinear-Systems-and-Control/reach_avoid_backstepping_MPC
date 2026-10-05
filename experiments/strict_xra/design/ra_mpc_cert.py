"""Algorithm 1 (TAC revision) for the Dubins car with an arbitrary certificate (cert.Cert)."""
import os, sys
import numpy as np, sympy as sp, casadi as ca
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "python"))
import ra_mpc_dubins as rm


def terminal_stats(cert, delta, dt, n_y=120000):
    """Sampled one-step increase of Psi and sup of the stage cost on {Psi >= delta} \\ X_T."""
    x, _ = cert.sample(delta, n_y=n_y, seed=1)
    u = cert.uc(x); xn = cert.zoh(x, u, dt)
    eps = float((cert.Psi(xn) - cert.Psi(x)).min())
    ell = rm.Q_Y * ((x[0] - rm.TARGET_POS[0])**2 + (x[1] - rm.TARGET_POS[1])**2) + rm.R_U * (u**2).sum(0)
    return eps, float(ell.max()), float(np.abs(u).max())


class CertMPC(rm.DubinsRAMPC):
    def __init__(self, cert, dt=0.05, N=25, dfrac=0.3, omega=None, safe_margin=1e-3, u_max=None):
        self.c = cert
        self.u_max = np.asarray(rm.U_MAX if u_max is None else u_max, float)
        M_psi = 1.0 - cert.shift
        delta = dfrac * M_psi
        self.eps, self.sup_ell, self.uc_max = terminal_stats(cert, delta, dt)
        if omega is None:
            omega = 1.2 * self.sup_ell / self.eps if self.eps > 0 else 1e4
        super().__init__(dt=dt, N=N, delta=delta, omega=omega, M_psi=M_psi, safe_margin=safe_margin)
        expr, X = cert.psi_state_sym()
        xs = ca.SX.sym("x", 4)
        self.psi_ca = ca.Function("psi", [xs], [rm._s2c(sp.expand(expr), dict(zip(X, [xs[i] for i in range(4)])))])
        self.cscale = 1.0 / max(1.0, omega)
        self._build_model()
        # tight tolerances: with a large omega the cost is sensitive to dynamics-constraint residuals
        self.opts = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 300, "ipopt.tol": 1e-5,
                     "ipopt.constr_viol_tol": 1e-8, "ipopt.acceptable_tol": 1e-4, "ipopt.acceptable_constr_viol_tol": 1e-7,
                     "ipopt.acceptable_iter": 5}

    NSUB = 4   # RK4 sub-steps per sampling interval in the MPC model (removes the integrator mismatch)

    def _build_model(self):
        f1 = rm._rk4_ca(self.dt / self.NSUB)
        xc, uc = ca.SX.sym("x", 4), ca.SX.sym("u", 2); xn = xc
        for _ in range(self.NSUB):
            xn = f1(xn, uc)
        self.F = ca.Function("F", [xc, uc], [xn])

    def _step(self, x, u):
        for _ in range(self.NSUB):
            x = rm.rk4_np(x, u, self.dt / self.NSUB)
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
        J = sum(rm.stage_cost(X[:, k], U[:, k]) for k in range(U.shape[1]))
        return J + self.omega * (self.M_psi - float(self._psi(X[:, -1])[0]))

    def feasible(self, X, U, tol=1e-6):
        return (np.all(np.abs(U) <= self.u_max[:, None] + tol)
                and np.all(rm.safe_np(*X[:, 1:]) >= self.safe_margin - tol)
                and float(self._psi(X[:, -1])[0]) >= self.delta - tol)

    def solve(self, x0, N, U_warm=None):
        opti = ca.Opti()
        X = opti.variable(4, N + 1); U = opti.variable(2, N)
        x0 = np.asarray(x0, float)
        opti.subject_to(X[:, 0] == x0)
        J = 0; tgt = ca.DM(rm.TARGET_POS)
        for k in range(N):
            opti.subject_to(X[:, k + 1] == self.F(X[:, k], U[:, k]))
            opti.subject_to(rm.safe_ca(X[:, k]) >= self.safe_margin)
            dy = X[:2, k] - tgt
            J += rm.Q_Y * ca.dot(dy, dy) + rm.R_U * ca.dot(U[:, k], U[:, k])
        opti.subject_to(rm.safe_ca(X[:, N]) >= self.safe_margin)
        opti.subject_to(self.psi_ca(X[:, N]) >= self.delta)
        opti.subject_to(opti.bounded(-self.u_max[0], U[0, :], self.u_max[0]))
        opti.subject_to(opti.bounded(-self.u_max[1], U[1, :], self.u_max[1]))
        J += self.omega * (self.M_psi - self.psi_ca(X[:, N]))
        opti.minimize(self.cscale * J)
        if U_warm is not None:
            Uw = np.asarray(U_warm)[:, :N]
            opti.set_initial(U, Uw); opti.set_initial(X, self.rollout(x0, Uw))
        else:
            tgt_full = np.array([*rm.TARGET_POS, x0[2], 0.2])
            opti.set_initial(X, np.array([(1 - k / N) * x0 + k / N * tgt_full for k in range(N + 1)]).T)
            opti.set_initial(U, 0)
        opti.solver("ipopt", self.opts)
        try:
            sol = opti.solve()
            Us = np.atleast_2d(sol.value(U)).reshape(2, N)
            return True, Us, self.rollout(x0, Us)      # true nominal trajectory, not the solver iterate
        except Exception:
            if U_warm is None and not getattr(self, "_retry", False):
                # second initial guess at k = 0: solution of the terminal-target (vanilla) problem
                ok, Uv = vanilla_solve(self, x0, N)
                if ok:
                    self._retry = True
                    try:
                        return self.solve(x0, N, U_warm=Uv)
                    finally:
                        self._retry = False
            return False, None, None

    def candidate(self, U_star, X_star):
        xN = X_star[:, -1].reshape(4, 1)
        uN = self.c.uc(xN)[:, 0] if xN[3, 0] > 1e-6 else np.zeros(2)
        Uc = np.hstack([U_star[:, 1:], uN[:, None]])
        Xc = self.rollout(X_star[:, 1], Uc)
        touch = np.any(rm.tgt_np(*Xc) <= 0)
        xs = Xc[:, -2].copy()
        for _ in range(10):
            xs = rm.rk4_np(xs, uN, self.dt / 10)
            touch = touch or float(rm.tgt_np(*xs)) <= 0
        return Uc, Xc, bool(touch), bool(np.all(np.abs(uN) <= self.u_max))


def vanilla_solve(mpc, x0, N, Qf=80.0):
    um = getattr(mpc, "u_max", rm.U_MAX)
    """Terminal-target MPC of the paper's baseline (same dynamics, cost and safety constraint)."""
    opti = ca.Opti(); X = opti.variable(4, N + 1); U = opti.variable(2, N)
    x0 = np.asarray(x0, float); opti.subject_to(X[:, 0] == x0); J = 0; tgt = ca.DM(rm.TARGET_POS)
    for k in range(N):
        opti.subject_to(X[:, k + 1] == mpc.F(X[:, k], U[:, k])); opti.subject_to(rm.safe_ca(X[:, k]) >= 0)
        dy = X[:2, k] - tgt; J += rm.Q_Y * ca.dot(dy, dy) + rm.R_U * ca.dot(U[:, k], U[:, k])
    opti.subject_to(rm.safe_ca(X[:, N]) >= 0)
    opti.subject_to(X[1, N]**2 + (2 * (X[0, N] + 1.7))**2 - 0.4 <= 0)
    opti.subject_to(opti.bounded(-um[0], U[0, :], um[0])); opti.subject_to(opti.bounded(-um[1], U[1, :], um[1]))
    dyN = X[:2, N] - tgt; J += Qf * ca.dot(dyN, dyN); opti.minimize(J)
    tgt_full = np.array([*rm.TARGET_POS, x0[2], 0.0])
    opti.set_initial(X, np.array([(1 - k / N) * x0 + k / N * tgt_full for k in range(N + 1)]).T); opti.set_initial(U, 0)
    opti.solver("ipopt", mpc.opts)
    try:
        sol = opti.solve(); return True, np.atleast_2d(sol.value(U)).reshape(2, N)
    except Exception:
        return False, None


def vanilla_feasible(x0, dt, N, u_max=None):
    class _M: pass
    m = _M(); m.F = rm._rk4_ca(dt); m.u_max = np.asarray(rm.U_MAX if u_max is None else u_max, float); m.opts = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 500, "ipopt.tol": 1e-4, "ipopt.acceptable_tol": 1e-3}
    return vanilla_solve(m, x0, N)[0]


def run_one(x0, k1_file, mu, scale, dt, N, dfrac, omega, T_max, plant, u_max=None):
    from cert import Cert
    c = Cert(k1_file, mu=mu)
    if scale != 1.0:
        c = c.scaled(scale)
    mpc = CertMPC(c, dt=dt, N=N, dfrac=dfrac, omega=omega, u_max=u_max)
    log = mpc.run(x0, T_max=T_max, plant=plant)
    log["uc_max_terminal"] = mpc.uc_max
    log["check"] = rm.corollary_check(log, dt)
    log["omega"] = mpc.omega; log["eps"] = mpc.eps
    return log
