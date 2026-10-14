"""Reach-avoid MPC for the Dubins car, TAC revision (Algorithm 1 of the paper).

Changes w.r.t. the LCSS notebooks:
  * terminal set {Psi >= delta}                       (Theorem 2)
  * cost  sum_i l(x_i,u_i) + omega (M_Psi - Psi(x_N)) (Remark 4; no Q_f term)
  * warm start = shifted candidate {u*_1..u*_{N-1}, u^c(x*_N)}  (Alg. 1)
  * applied sequence = cheaper of {solver output, feasible candidate}
  * adaptive horizon: N <- max(N-1, 1) if the candidate rollout touches X_T
  * logs J*_k, candidate cost, horizon, and checks Corollary 1.
"""
import os, sys
import numpy as np
import sympy as sp
import casadi as ca

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "controllers"))
from sop_bounded_control_dubins_car_result import u_opt, certificate_opt  # noqa: E402

# ── problem data (unchanged from the LCSS notebooks) ─────────────────────────
U_MAX = np.array([5.0, 5.0])
TARGET_POS = np.array([-1.7, 0.0])
Q_Y, R_U = 5.0, 0.05
ELL_LOWER = Q_Y * 0.1  # ||h(x)-y*||^2 > 0.1 outside the target ellipse

x1, x2, th, v = sp.symbols("x1 x2 th v")
_X = [x1, x2, th, v]
_psi_sym = sp.sympify(certificate_opt)
_y1, _y2 = x1, x2
_safe_sym = 1e-3 * (-((_y2) ** 2 + (2 * (_y1 + 1.7)) ** 2 - 0.4) + 300) * (
    -(_y1**4 + _y2**4 - 16) * (_y1**4 + _y2**4 - 4))
_tgt_sym = _y2**2 + (2 * (_y1 + 1.7)) ** 2 - 0.4

psi_np = sp.lambdify(_X, _psi_sym, "numpy")
safe_np = sp.lambdify(_X, _safe_sym, "numpy")
tgt_np = sp.lambdify(_X, _tgt_sym, "numpy")
_uc_np = sp.lambdify(_X, [sp.sympify(e) for e in u_opt], "numpy")


def uc_np(x):
    return np.array(_uc_np(*x), dtype=float).flatten()


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
psi_ca = ca.Function("psi", [_xs], [_s2c(_psi_sym, _m)])
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
        # discretisation error w.r.t. the continuous plant (see TAC_EXPERIMENT_UPDATE.md)
        self.dt, self.N0, self.delta, self.omega, self.M_psi = dt, N, delta, omega, M_psi
        self.safe_margin = safe_margin
        self.F = _rk4_ca(dt)
        self.opts = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 500,
                     "ipopt.tol": 1e-4, "ipopt.acceptable_tol": 1e-3}

    # cost of an input sequence from x0 under the nominal (RK4) model
    def rollout(self, x0, U):
        X = [np.asarray(x0, float)]
        for k in range(U.shape[1]):
            X.append(rk4_np(X[-1], U[:, k], self.dt))
        return np.array(X).T

    def cost(self, X, U):
        J = sum(stage_cost(X[:, k], U[:, k]) for k in range(U.shape[1]))
        return J + self.omega * (self.M_psi - float(psi_np(*X[:, -1])))

    def feasible(self, X, U, tol=1e-6):
        return (np.all(np.abs(U) <= U_MAX[:, None] + tol)
                and np.all(safe_np(*X[:, 1:]) >= self.safe_margin - tol)
                and float(psi_np(*X[:, -1])) >= self.delta - tol)

    def solve(self, x0, N, U_warm=None):
        opti = ca.Opti()
        X = opti.variable(4, N + 1); U = opti.variable(2, N)
        x0 = np.asarray(x0, float)
        opti.subject_to(X[:, 0] == x0)
        J = 0
        tgt = ca.DM(TARGET_POS)
        for k in range(N):
            opti.subject_to(X[:, k + 1] == self.F(X[:, k], U[:, k]))
            opti.subject_to(safe_ca(X[:, k]) >= self.safe_margin)
            dy = X[:2, k] - tgt
            J += Q_Y * ca.dot(dy, dy) + R_U * ca.dot(U[:, k], U[:, k])
        opti.subject_to(safe_ca(X[:, N]) >= self.safe_margin)
        opti.subject_to(psi_ca(X[:, N]) >= self.delta)
        opti.subject_to(opti.bounded(-U_MAX[0], U[0, :], U_MAX[0]))
        opti.subject_to(opti.bounded(-U_MAX[1], U[1, :], U_MAX[1]))
        J += self.omega * (self.M_psi - psi_ca(X[:, N]))
        opti.minimize(J)
        if U_warm is not None:
            Uw = np.asarray(U_warm)[:, :N]
            opti.set_initial(U, Uw); opti.set_initial(X, self.rollout(x0, Uw))
        else:
            tgt_full = np.array([*TARGET_POS, x0[2], 0.0])
            opti.set_initial(X, np.array([(1 - k / N) * x0 + k / N * tgt_full for k in range(N + 1)]).T)
            opti.set_initial(U, 0)
        opti.solver("ipopt", self.opts)
        try:
            sol = opti.solve()
            Us = np.atleast_2d(sol.value(U)).reshape(2, N)
            return True, Us, np.atleast_2d(sol.value(X)).reshape(4, N + 1)
        except Exception:
            return False, None, None

    def candidate(self, U_star, X_star):
        """{u*_1..u*_{N-1}, u^c(x*_N)} and whether its rollout touches X_T."""
        uN = uc_np(X_star[:, -1])
        Uc = np.hstack([U_star[:, 1:], uN[:, None]])
        Xc = self.rollout(X_star[:, 1], Uc)
        touch = np.any(tgt_np(*Xc) <= 0)
        # last (u^c) segment: check integrator sub-steps too
        xs = Xc[:, -2].copy()
        for _ in range(10):
            xs = rk4_np(xs, uN, self.dt / 10)
            touch = touch or float(tgt_np(*xs)) <= 0
        return Uc, Xc, bool(touch), bool(np.all(np.abs(uN) <= U_MAX))

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


def run_one(x0, dt=0.05, N=25, delta=0.1, omega=1000.0, T_max=10.0, plant="nominal"):
    """Picklable helper for joblib: closed loop + Corollary 1 check."""
    mpc = DubinsRAMPC(dt=dt, N=N, delta=delta, omega=omega)
    log = mpc.run(x0, T_max=T_max, plant=plant)
    log["check"] = corollary_check(log, dt)
    return log
