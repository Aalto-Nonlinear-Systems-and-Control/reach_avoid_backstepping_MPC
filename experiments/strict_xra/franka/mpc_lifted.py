"""Lifted RK4 for the Panda MPC: the stage accelerations are decision variables with M(q_s) a_s + b(q_s, v_s) = u, so no
matrix inverse appears in the NLP (same RK4 scheme as the rollout). Subclass of mpc_franka.FrankaMPC; only _build and
the bounds/initial guess of the extra variables differ. Interpreted (no compilation)."""
import os, numpy as np, casadi as ca
import mpc_franka as mp
from mpc_franka import FrankaMPC, NQ, TARGET, RT, Q_Y, R_U, QF_VAN, target_q
from franka_model import M_sx, B_sx, q as Qs, qd as QDs

_x = ca.SX.sym("x", 2 * NQ); _a = ca.SX.sym("a", NQ); _u = ca.SX.sym("u", NQ)
RES = ca.Function("res", [_x, _a, _u], [ca.substitute(ca.substitute(M_sx, Qs, _x[:NQ]), QDs, _x[NQ:]) @ _a + ca.substitute(ca.substitute(B_sx, Qs, _x[:NQ]), QDs, _x[NQ:]) - _u])


class LiftedMPC(FrankaMPC):
    def _lift_fun(self):
        """Per interval: x_k, u_k, stage accelerations A (NQ x 4 NSUB) -> next state, residuals, mid states."""
        x = ca.SX.sym("x", 2 * NQ); u = ca.SX.sym("u", NQ); A = ca.SX.sym("A", NQ, 4 * self.NSUB); h = self.dt / self.NSUB; res = []; mids = []; z = x
        for i in range(self.NSUB):
            a1, a2, a3, a4 = (A[:, 4 * i + j] for j in range(4)); q0, v0 = z[:NQ], z[NQ:]
            x1 = z; x2 = ca.vertcat(q0 + h / 2 * v0, v0 + h / 2 * a1); v2 = v0 + h / 2 * a1
            x3 = ca.vertcat(q0 + h / 2 * v2, v0 + h / 2 * a2); v3 = v0 + h / 2 * a2
            x4 = ca.vertcat(q0 + h * v3, v0 + h * a3); v4 = v0 + h * a3
            res += [RES(x1, a1, u), RES(x2, a2, u), RES(x3, a3, u), RES(x4, a4, u)]
            z = ca.vertcat(q0 + h / 6 * (v0 + 2 * v2 + 2 * v3 + v4), v0 + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4))
            if i < self.NSUB - 1: mids.append(z)
        return ca.Function("lift", [x, u, A], [z, ca.vertcat(*res), ca.horzcat(*mids) if mids else ca.SX.zeros(2 * NQ, 0)])

    def _build(self):
        N, n, m_ = self.N0, 2 * NQ, NQ; L = self._lift_fun(); X = ca.MX.sym("X", n, N + 1); U = ca.MX.sym("U", m_, N); A = ca.MX.sym("A", m_, 4 * self.NSUB * N)
        x0 = ca.MX.sym("x0", n); a = ca.MX.sym("a", N); e = ca.MX.sym("e", N + 1); psi = self.c.F["psi"]; Psi = self.c.F["Psi"]; yfun = self.c.F["y"]; tg = ca.DM(TARGET)
        Ar = ca.horzsplit(A, [4 * self.NSUB * k for k in range(N + 1)]); Z, R, Mid = [], [], []
        Lm = L.map(N); Zs, Rs, Ms = Lm(X[:, :N], U, ca.reshape(A, m_ * 4 * self.NSUB, N) if False else ca.horzcat(*[ca.vec(Ar[k]) for k in range(N)]).reshape((m_ * 4 * self.NSUB, N)) if False else None) if False else (None, None, None)
        for k in range(N):
            zk, rk, mk = L(X[:, k], U[:, k], Ar[k]); Z.append(zk); R.append(rk); Mid.append(mk)
        g = [X[:, 0] - x0, ca.vec(X[:, 1:] - ca.horzcat(*Z)), ca.vertcat(*R)]
        if self.NSUB > 1: g.append(ca.vec(psi.map(N * (self.NSUB - 1))(ca.horzcat(*Mid))))
        g.append(ca.vec(psi.map(N)(X[:, 1:])))
        Y = yfun.map(N)(X[:, :N]); dY = Y[:2, :] - ca.repmat(tg, 1, N); ell = Q_Y * ca.sum1(dY**2) + R_U * ca.sum1(U**2); J = ca.dot(a, ell.T)
        if self.terminal == "cert": Pn = Psi.map(N + 1)(X); g.append(ca.vec(Pn)); J += self.omega * ca.dot(e, (self.M_psi - Pn).T)
        else: yN = yfun(X[:, N]); g.append((yN[0] - tg[0])**2 + (yN[1] - tg[1])**2 - RT**2); dyN = yN[:2] - tg; J += QF_VAN * ca.dot(dyN, dyN)
        self._nA = m_ * 4 * self.NSUB * N
        return dict(x=ca.vertcat(ca.vec(X), ca.vec(U), ca.vec(A)), p=ca.vertcat(x0, a, e), f=self.cscale * J, g=ca.vertcat(*g))

    def _solver(self):
        if not hasattr(self, "_S"): self._S = ca.nlpsol("S", "ipopt", self._build(), self.opts)
        return self._S

    def _accels(self, x0, U):
        """Stage accelerations of the explicit RK4 rollout (initial guess of A)."""
        h = self.dt / self.NSUB; f = lambda z, u: np.array(mp.FDYN(z, u)).ravel(); x = np.asarray(x0, float); cols = []
        for k in range(U.shape[1]):
            for _ in range(self.NSUB):
                k1 = f(x, U[:, k]); k2 = f(x + h / 2 * k1, U[:, k]); k3 = f(x + h / 2 * k2, U[:, k]); k4 = f(x + h * k3, U[:, k])
                cols += [k1[NQ:], k2[NQ:], k3[NQ:], k4[NQ:]]; x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return np.array(cols).T

    def solve(self, x0, N, U_warm=None):
        S = self._solver(); N0, n, m_ = self.N0, 2 * NQ, NQ; x0 = np.asarray(x0, float); Na = int(N)
        a = (np.arange(N0) < Na).astype(float); e = (np.arange(N0 + 1) == Na).astype(float)
        Uw = np.zeros((m_, N0))
        if U_warm is not None: Uw[:, :min(Na, U_warm.shape[1])] = np.asarray(U_warm)[:, :Na]
        else: Uw = np.zeros((m_, N0))
        Xw = self.rollout(x0, Uw); Aw = self._accels(x0, Uw)
        INF = 1e20; lbg = [np.zeros(n * (N0 + 1)), np.zeros(m_ * 4 * self.NSUB * N0)]; ubg = [np.zeros(n * (N0 + 1)), np.zeros(m_ * 4 * self.NSUB * N0)]
        if self.NSUB > 1:
            lm = np.array([[(min(self.safe_margin, self.safe_np(x0)) if k == 0 else self.safe_margin) if k < Na else -INF] * (self.NSUB - 1) for k in range(N0)]).ravel(); lbg.append(lm); ubg.append(INF * np.ones_like(lm))
        ln = np.where(np.arange(1, N0 + 1) <= Na, self.safe_margin, -INF); lbg.append(ln); ubg.append(INF * np.ones(N0))
        if self.terminal == "cert": lt = np.where(np.arange(N0 + 1) == Na, self.delta, -INF); lbg.append(lt); ubg.append(INF * np.ones(N0 + 1))
        else: lbg.append([-INF]); ubg.append([0.0])
        ub = np.tile(self.u_max, N0); nA = m_ * 4 * self.NSUB * N0
        lbx = np.r_[-INF * np.ones(n * (N0 + 1)), -ub, -INF * np.ones(nA)]; ubx = np.r_[INF * np.ones(n * (N0 + 1)), ub, INF * np.ones(nA)]
        try:
            r = S(x0=np.r_[Xw.T.ravel(), Uw.T.ravel(), Aw.T.ravel()], p=np.r_[x0, a, e], lbg=np.concatenate(lbg), ubg=np.concatenate(ubg), lbx=lbx, ubx=ubx)
            st = S.stats(); ok = bool(st["success"]) or st["return_status"] == "Solved_To_Acceptable_Level"
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
