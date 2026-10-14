"""Algorithm 1 (RA-MPC) and the vanilla MPC for the reduced Panda with a certificate of cert_franka.CertFranka.
Same structure as ../arm_top/ra_mpc_top.py: RK4 sub-step model, state constraint psi(h(x)) >= margin at every sub-step
node, terminal set {Psi >= delta}, terminal cost omega (M_Psi - Psi), shifted candidate ending with u^c, keep-the-cheaper,
adaptive horizon, fallback initial guesses. The plant is MuJoCo (model/panda_task.xml)."""
import os, sys, time, numpy as np, casadi as ca, mujoco
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from franka_model import MJ, NQ, Y_f, J_f, ik
from cert_franka import CertFranka, FDYN, TARGET, RT, Z0, phi_xy, sample
Q_Y, R_U, QF_VAN = 500.0, 0.04, 8000.0
CL = os.path.join(os.path.dirname(sys.executable), "arm64-apple-darwin20.0.0-clang"); SDK = "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"   # conda clang, C code of the NLP
ELL_LOWER = Q_Y * RT**2
_TGT_Q = None


def target_q():
    global _TGT_Q
    if _TGT_Q is None:
        y = np.array([TARGET[0], TARGET[1], Z0, 0.0]); _TGT_Q, _ = ik(y, np.array([np.arctan2(-y[0], y[1]), 0.3, -2.0, 2.3]))
    return _TGT_Q


def yx(x): return np.array(Y_f(np.asarray(x)[:NQ])).ravel()
def tgt_np(x): y = yx(x); return float(phi_xy(y))
def stage_cost(x, u): y = yx(x); return Q_Y * float((y[0] - TARGET[0])**2 + (y[1] - TARGET[1])**2) + R_U * float(np.dot(u, u))


def rk4(x, u, h):
    f = lambda z: np.array(FDYN(z, u)).ravel()
    k1 = f(x); k2 = f(x + h / 2 * k1); k3 = f(x + h / 2 * k2); k4 = f(x + h * k3); return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


class Plant:
    """MuJoCo plant: the input is held over the period (motor torques), 0.5 ms RK4 steps."""
    def __init__(self): self.d = mujoco.MjData(MJ)
    def step(self, x, u, T):
        d = self.d; d.qpos[:] = x[:NQ]; d.qvel[:] = x[NQ:]; d.ctrl[:] = u
        for _ in range(int(round(T / MJ.opt.timestep))): mujoco.mj_step(MJ, d)
        return np.r_[d.qpos.copy(), d.qvel.copy()]


def uc_guess(cert, x0, N, u_max, step):
    x = np.asarray(x0, float).copy(); U = np.zeros((NQ, N))
    for k in range(N): U[:, k] = np.clip(cert.uc(x.reshape(-1, 1))[:, 0], -u_max, u_max); x = step(x, U[:, k])
    return U


def steer_guess(x0, N, u_max, step):
    x = np.asarray(x0, float).copy(); U = np.zeros((NQ, N)); qt = target_q(); Kp = u_max / 0.5; Kd = u_max / 1.0
    for k in range(N): U[:, k] = np.clip(Kp * (qt - x[:NQ]) - Kd * x[NQ:], -u_max, u_max); x = step(x, U[:, k])
    return U


class FrankaMPC:
    def __init__(self, cert, dt=0.01, N=60, dfrac=0.5, omega=None, safe_margin=0.005, u_max=(4.5, 6.5, 2.9, 0.3), nsub=2, terminal="cert", opts=None, hessian=None):
        self.c, self.dt, self.N0, self.terminal = cert, dt, N, terminal; self.u_max = np.asarray(u_max, float); self.NSUB = int(nsub); self.safe_margin = safe_margin
        self.M_psi = 1.0 - cert.shift; self.delta = dfrac * self.M_psi
        if terminal == "cert":
            self.eps, self.sup_ell = self.terminal_stats(); self.omega = omega if omega is not None else (1.2 * self.sup_ell / self.eps if self.eps > 0 else 1e5)
        else: self.eps, self.sup_ell, self.omega = np.nan, np.nan, 0.0
        xs, us = ca.MX.sym("x", 2 * NQ), ca.MX.sym("u", NQ); h = dt / self.NSUB            # MX: FDYN stays a call, its derivatives are generated once
        def rk(x, u):
            k1 = FDYN(x, u); k2 = FDYN(x + h / 2 * k1, u); k3 = FDYN(x + h / 2 * k2, u); k4 = FDYN(x + h * k3, u); return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        xn = xs; mids = []
        for i in range(self.NSUB):
            xn = rk(xn, us)
            if i < self.NSUB - 1: mids.append(xn)
        self.F = ca.Function("F", [xs, us], [xn]); self.Fmid = ca.Function("Fm", [xs, us], [ca.horzcat(*mids)]) if mids else None
        self.cscale = 1.0 / max(1.0, self.omega)
        self.opts = opts or {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 300, "ipopt.tol": 1e-5, "ipopt.constr_viol_tol": 1e-8,
                             "ipopt.acceptable_tol": 1e-4, "ipopt.acceptable_constr_viol_tol": 1e-7, "ipopt.acceptable_iter": 5}
        hessian = hessian or os.environ.get("MPC_HESSIAN", "exact")
        if hessian == "limited-memory": self.opts["ipopt.hessian_approximation"] = "limited-memory"

    def terminal_stats(self, n=6000):
        x, _ = sample(self.c, self.delta, n=n, seed=1); u = self.c.uc(x)
        from cert_franka import zoh
        xn = zoh(x, u, self.dt); keep = phi_xy(self.c.y(xn)) > 0; eps = float((self.c.Psi(xn) - self.c.Psi(x))[keep].min())
        Y = self.c.y(x); ell = Q_Y * ((Y[0] - TARGET[0])**2 + (Y[1] - TARGET[1])**2) + R_U * (u**2).sum(0); return eps, float(ell.max())

    def _step(self, x, u):
        for _ in range(self.NSUB): x = rk4(x, u, self.dt / self.NSUB)
        return x

    def rollout(self, x0, U):
        X = [np.asarray(x0, float)]
        for k in range(U.shape[1]): X.append(self._step(X[-1], U[:, k]))
        return np.array(X).T

    def _Psi(self, x): return float(self.c.Psi(np.asarray(x, float).reshape(-1, 1))[0])
    def safe_np(self, x): return float(self.c.psi(np.asarray(x, float).reshape(-1, 1))[0])

    def cost(self, X, U):
        J = sum(stage_cost(X[:, k], U[:, k]) for k in range(U.shape[1]))
        if self.terminal == "cert": return J + self.omega * (self.M_psi - self._Psi(X[:, -1]))
        y = yx(X[:, -1]); return J + QF_VAN * float((y[0] - TARGET[0])**2 + (y[1] - TARGET[1])**2)

    def feasible(self, X, U, tol=1e-6):
        if not np.all(np.abs(U) <= self.u_max[:, None] + tol): return False
        for k in range(1, X.shape[1]):
            if self.safe_np(X[:, k]) < self.safe_margin - tol: return False
        return self._Psi(X[:, -1]) >= self.delta - tol if self.terminal == "cert" else tgt_np(X[:, -1]) <= tol

    # ── one NLP of horizon N0, compiled once to C and cached on disk; the active horizon Na <= N0 is set at solve
    #    time through parameters (stage weights a_k, terminal indicator e_k) and constraint bounds (equivalent problem) ──
    def _key(self):
        import hashlib
        tag = f"{self.terminal}|{self.c.s:.10g}|{self.c.mu0:.10g}|{self.dt}|{self.N0}|{self.NSUB}|{self.omega:.10g}|{self.delta:.10g}|{self.safe_margin}|{tuple(self.u_max)}|{Q_Y}|{R_U}|{QF_VAN}|{self.cscale:.10g}"
        return hashlib.sha1(tag.encode()).hexdigest()[:12]

    def _build(self):
        """MX graph with one call per interval (dynamics, mid-state psi, node psi, outputs, Psi)."""
        N, n, m_ = self.N0, 2 * NQ, NQ; X = ca.MX.sym("X", n, N + 1); U = ca.MX.sym("U", m_, N); x0 = ca.MX.sym("x0", n); a = ca.MX.sym("a", N); e = ca.MX.sym("e", N + 1)
        psi = self.c.F["psi"]; Psi = self.c.F["Psi"]; yfun = self.c.F["y"]; tg = ca.DM(TARGET); g = [X[:, 0] - x0]; J = 0
        for k in range(N): g.append(X[:, k + 1] - self.F(X[:, k], U[:, k]))
        if self.Fmid is not None:
            for k in range(N):
                Xm = self.Fmid(X[:, k], U[:, k]); g += [psi(Xm[:, j]) for j in range(self.NSUB - 1)]
        g += [psi(X[:, k]) for k in range(1, N + 1)]
        for k in range(N): dy = yfun(X[:, k])[:2] - tg; J += a[k] * (Q_Y * ca.dot(dy, dy) + R_U * ca.dot(U[:, k], U[:, k]))
        if self.terminal == "cert":                                  # terminal state x_T = X e (e one-hot at the active horizon): Psi evaluated once
            xT = ca.MX.sym("xT", n); g += [xT - X @ e, Psi(xT)]; J += self.omega * (self.M_psi - Psi(xT))
            w = ca.vertcat(ca.vec(X), ca.vec(U), xT)
        else:
            yN = yfun(X[:, N]); g.append((yN[0] - tg[0])**2 + (yN[1] - tg[1])**2 - RT**2); dyN = yN[:2] - tg; J += QF_VAN * ca.dot(dyN, dyN); w = ca.vertcat(ca.vec(X), ca.vec(U))
        return dict(x=w, p=ca.vertcat(x0, a, e), f=self.cscale * J, g=ca.vertcat(*g))

    def _solver(self):
        """Compiled in process with the conda clang (CasADi JIT), once per controller object."""
        if not hasattr(self, "_S"):
            o = dict(self.opts); o.update({"jit": True, "compiler": "shell", "jit_options": {"compiler": CL, "linker": CL, "compiler_flags": ["-O1", "-fPIC", "-isysroot", SDK], "linker_flags": ["-isysroot", SDK, "-fPIC"], "verbose": False}})
            self._S = ca.nlpsol("S", "ipopt", self._build(), o)
        return self._S

    def solve(self, x0, N, U_warm=None):
        S = self._solver(); N0, n, m_ = self.N0, 2 * NQ, NQ; x0 = np.asarray(x0, float); Na = int(N)
        a = (np.arange(N0) < Na).astype(float); e = (np.arange(N0 + 1) == Na).astype(float)
        if U_warm is not None: Uw = np.zeros((m_, N0)); Uw[:, :min(Na, U_warm.shape[1])] = np.asarray(U_warm)[:, :Na]
        else: Uw = np.zeros((m_, N0))
        Xw = self.rollout(x0, Uw) if U_warm is not None else np.array([(1 - min(k, Na) / Na) * x0 + min(k, Na) / Na * np.r_[target_q(), np.zeros(NQ)] for k in range(N0 + 1)]).T
        INF = 1e20; lbg, ubg = [np.zeros(n * (N0 + 1))], [np.zeros(n * (N0 + 1))]
        if self.Fmid is not None:
            lm = np.array([[(min(self.safe_margin, self.safe_np(x0)) if k == 0 else self.safe_margin) if k < Na else -INF] * (self.NSUB - 1) for k in range(N0)]).ravel(); lbg.append(lm); ubg.append(INF * np.ones_like(lm))
        ln = np.where(np.arange(1, N0 + 1) <= Na, self.safe_margin, -INF); lbg.append(ln); ubg.append(INF * np.ones(N0))
        if self.terminal == "cert": lbg += [np.zeros(n), [self.delta]]; ubg += [np.zeros(n), [INF]]
        else: lbg.append([-INF]); ubg.append([0.0])
        ub = np.tile(self.u_max, N0); lbx = np.r_[-INF * np.ones(n * (N0 + 1)), -ub]; ubx = np.r_[INF * np.ones(n * (N0 + 1)), ub]; w0 = np.r_[Xw.T.ravel(), Uw.T.ravel()]
        if self.terminal == "cert": lbx = np.r_[lbx, -INF * np.ones(n)]; ubx = np.r_[ubx, INF * np.ones(n)]; w0 = np.r_[w0, Xw[:, Na]]
        try:
            r = S(x0=w0, p=np.r_[x0, a, e], lbg=np.concatenate(lbg), ubg=np.concatenate(ubg), lbx=lbx, ubx=ubx)
            ok = bool(S.stats()["success"]) or S.stats()["return_status"] in ("Solved_To_Acceptable_Level",)
        except Exception: ok = False
        if ok:
            w = np.array(r["x"]).ravel(); Us = w[n * (N0 + 1): n * (N0 + 1) + m_ * N0].reshape(N0, m_).T[:, :Na]; return True, Us, self.rollout(x0, Us)
        if U_warm is None and not getattr(self, "_retry", False):
            self._retry = True
            try:
                for Ug in self._guesses(x0, Na):
                    ok, Us, Xs = self.solve(x0, Na, U_warm=Ug)
                    if ok: return ok, Us, Xs
            finally: self._retry = False
        return False, None, None

    def _guesses(self, x0, N):
        yield uc_guess(self.c, x0, N, self.u_max, self._step); yield steer_guess(x0, N, self.u_max, self._step)
        if self.terminal == "cert":
            if not hasattr(self, "_van"): self._van = FrankaMPC(self.c, self.dt, N, self.delta / self.M_psi, 0.0, self.safe_margin, self.u_max, self.NSUB, "target"); self._van._retry = True
            ok, Uv, _ = self._van.solve(x0, N)
            if ok: yield Uv
            if not getattr(self, "_no_continuation", False):
                if not hasattr(self, "_cont"):
                    self._cont = {}
                    for f in (0.25, 0.5):
                        self._cont[f] = FrankaMPC(self.c.scaled(f), self.dt, N, self.delta / self.M_psi, self.omega, self.safe_margin, self.u_max, self.NSUB, "cert"); self._cont[f]._no_continuation = True
                Uw = None
                for f in (0.25, 0.5):
                    ok, Uw, _ = self._cont[f].solve(x0, N, U_warm=Uw)
                    if not ok: return
                yield Uw

    def candidate(self, U_star, X_star):
        xN = X_star[:, -1].reshape(-1, 1); uN = self.c.uc(xN)[:, 0]; Uc = np.hstack([U_star[:, 1:], uN[:, None]]); Xc = self.rollout(X_star[:, 1], Uc)
        touch = any(tgt_np(Xc[:, k]) <= 0 for k in range(Xc.shape[1])); return Uc, Xc, bool(touch), bool(np.all(np.abs(uN) <= self.u_max))

    def run(self, x0, T_max=6.0, keep_cheaper=True):
        plant = Plant(); x = np.asarray(x0, float); N = self.N0; U_warm = None
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
            else: U_warm = np.hstack([Us[:, 1:], Us[:, -1:]])
            x = plant.step(x, Us[:, 0], self.dt)
            if not np.all(np.isfinite(x)): log["status"] = "nonfinite"; break
            log["t"].append(log["t"][-1] + self.dt); log["x"].append(x.copy())
        else:
            if tgt_np(x) <= 0: log["status"] = "reached"
        for key in ("t", "J", "N", "ell"): log[key] = np.array(log[key])
        log["x"] = np.array(log["x"]); log["u"] = np.array(log["u"]).reshape(-1, NQ); J = log["J"]
        log["check"] = dict(n_steps=len(J), t_entry=log["t"][-1] if log["status"] == "reached" else np.nan, t_bound=(np.floor(J[0] / ELL_LOWER) + 1) * self.dt if len(J) else np.nan,
                            n_decrease_viol=int(((np.diff(J) + ELL_LOWER) > 1e-6).sum()) if len(J) > 1 else 0)
        log["check"]["bound_ok"] = bool(log["check"]["t_entry"] <= log["check"]["t_bound"] + 1e-9) if log["status"] == "reached" else False
        log["omega"], log["eps"] = self.omega, self.eps; return log
