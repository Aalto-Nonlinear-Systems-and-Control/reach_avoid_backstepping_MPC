"""Backstepping ECGBF certificate of the reduced Panda for one kitting task (default red_fwd_sq of ../arm_top).
Outputs y = (y1, y2, y3, y4) in metres and radians (franka_model). Safe-set function
    psi(y) = psi_xy(y1/ALPHA, y2/ALPHA) - CZ (y3 - Z0)^2 - CO y4^2,
psi_xy the normalised planar function of the layout and k1 = (ALPHA k1_xy(y/ALPHA), -KZ (y3 - Z0), -KO y4).
The first-layer condition holds because the planar part is the SOS certificate of the layout (spatial scaling is exact)
and the added terms satisfy (2 KZ + lambda) CZ (y3 - Z0)^2 >= 0, (2 KO + lambda) CO y4^2 >= 0.
Time scaling (k1, lambda, delta_s, mu) -> (s k1, s lambda, s delta_s, s^2 mu) is exact: no gravity, no damping."""
import os, sys, json, numpy as np, sympy as sp, casadi as ca
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from franka_model import Y_sx, J_sx, JDQD_sx, M_sx, B_sx, q as Qs, qd as QDs, NQ, MJ, Y_f, J_f
TASK = os.environ.get("TASK", "red_fwd_sq"); LAM_TAG = os.environ.get("LAM", "0.1")
ATOP = os.path.join(HERE, "..", "arm_top"); SPEC = json.load(open(os.path.join(ATOP, "tasks", f"{TASK}.json")))
K1FILE = os.path.join(ATOP, "k1", f"{TASK}_d8_l{LAM_TAG}.txt")
ALPHA, Z0 = 0.10, 0.30; CZ, CO, KZ, KO = 1 / 0.03**2, 1 / 0.10**2, 1.0, 1.0
TARGET = ALPHA * np.array(SPEC["target"]); RT = ALPHA * SPEC["RT"]; SOURCE = ALPHA * np.array(SPEC["source"]); BOX = ALPHA * np.array(SPEC["box"])
_u1, _u2 = sp.symbols("y1 y2")


def _p(txt): return sp.sympify(txt.replace("^", "**"), locals={"y1": _u1, "y2": _u2})


def s2c(e, mp):
    if isinstance(e, sp.Symbol): return mp[e]
    if e.is_Number: return float(e)
    if isinstance(e, sp.Add): return sum(s2c(a, mp) for a in e.args)
    if isinstance(e, sp.Mul):
        r = 1.0
        for a in e.args: r = r * s2c(a, mp)
        return r
    if isinstance(e, sp.Pow): return s2c(e.args[0], mp) ** float(e.args[1])
    raise ValueError(type(e))


_t = open(K1FILE).read().split("\n"); K1XY = [sp.expand(_p(_t[0])), sp.expand(_p(_t[1]))]; LAM0, DS0 = float(_t[2]), max(float(_t[3]), 0.0); PSIXY = sp.expand(_p(_t[4]) / float(_t[5]))


def phi_xy(Y): return (Y[0] - TARGET[0])**2 + (Y[1] - TARGET[1])**2 - RT**2            # target in the plane: phi < 0


class CertFranka:
    def __init__(self, mu=0.001, s=1.0):
        self.s, self.mu0 = s, mu; self.lam, self.ds, self.mu = s * LAM0, s * DS0, s * s * mu; self.shift = self.ds / self.lam
        ys = ca.SX.sym("y", 4); mp = {_u1: ys[0] / ALPHA, _u2: ys[1] / ALPHA}
        psi = s2c(PSIXY, mp) - CZ * (ys[2] - Z0)**2 - CO * ys[3]**2
        k1 = s * ca.vertcat(ALPHA * s2c(K1XY[0], mp), ALPHA * s2c(K1XY[1], mp), -KZ * (ys[2] - Z0), -KO * ys[3])
        self.psi_y = ca.Function("psi_y", [ys], [psi]); self.k1_y = ca.Function("k1_y", [ys], [k1]); grad_y = ca.Function("g", [ys], [ca.gradient(psi, ys)]); Jk_y = ca.Function("Jk", [ys], [ca.jacobian(k1, ys)])
        x = ca.SX.sym("x", 2 * NQ); qq, vv = x[:NQ], x[NQ:]; sub = lambda e: ca.substitute(ca.substitute(e, Qs, qq), QDs, vv)
        y = sub(Y_sx); Jm = sub(J_sx); jd = sub(JDQD_sx); Mm = sub(M_sx); bb = sub(B_sx); yd = Jm @ vv
        k = self.k1_y(y); e = yd - k; g = grad_y(y); Jk = Jk_y(y)
        Psi = self.psi_y(y) - self.shift - ca.dot(e, e) / (2 * self.mu)
        ydd_des = self.mu * g + Jk @ yd + 0.5 * self.lam * e; tau = Mm @ ca.solve(Jm, ydd_des - jd) + bb
        u = ca.SX.sym("u", NQ); ydd = Jm @ ca.solve(Mm, u - bb) + jd; Pdot = ca.dot(g, yd) - ca.dot(e, ydd - Jk @ yd) / self.mu
        self.F = dict(Psi=ca.Function("Psi", [x], [Psi]), uc=ca.Function("uc", [x], [tau]), Pdot=ca.Function("Pd", [x, u], [Pdot]), y=ca.Function("yx", [x], [y]),
                      psi=ca.Function("psix", [x], [self.psi_y(y)]), yd=ca.Function("ydx", [x], [yd]))
        self._maps = {}

    def scaled(self, s): return CertFranka(mu=self.mu0, s=self.s * s)

    def _b(self, name, *args):
        n = args[0].shape[1]; key = (name, n)
        if key not in self._maps: self._maps[key] = self.F[name].map(n, "thread", 8)
        return np.array(self._maps[key](*args))

    def Psi(self, X): return self._b("Psi", X).ravel()
    def uc(self, X): return self._b("uc", X)
    def Psidot(self, X, U): return self._b("Pdot", X, U).ravel()
    def y(self, X): return self._b("y", X)
    def psi(self, X): return self._b("psi", X).ravel()

    def psi_of_y(self, Y):
        n = Y.shape[1]; key = ("psiy", n)
        if key not in self._maps: self._maps[key] = self.psi_y.map(n, "thread", 8)
        return np.array(self._maps[key](Y)).ravel()

    def k1_of_y(self, Y):
        n = Y.shape[1]; key = ("k1y", n)
        if key not in self._maps: self._maps[key] = self.k1_y.map(n, "thread", 8)
        return np.array(self._maps[key](Y))


_YM, _JM = {}, {}


def ik_batch(Y, iters=40):
    """Batched Newton on h(q) = y, hand pointing along y4; seeds from the base yaw and a bent elbow."""
    n = Y.shape[1]
    if n not in _YM: _YM[n] = Y_f.map(n, "thread", 8); _JM[n] = J_f.map(n, "thread", 8)
    X, Yw = Y[1], -Y[0]; Q = np.vstack([np.arctan2(Yw, X), 0.3 * np.ones(n), -2.0 * np.ones(n), 2.3 * np.ones(n)])
    for _ in range(iters):
        E = np.array(_YM[n](Q)) - Y; Jb = np.array(_JM[n](Q)).reshape(4, n, 4).transpose(1, 0, 2)
        Q = Q - np.linalg.solve(Jb, E.T[:, :, None])[:, :, 0].T
    err = np.abs(np.array(_YM[n](Q)) - Y).max(0); return Q, err


def sample(cert, delta, n=20000, seed=0, boundary_frac=0.5, n_y=None):
    """States of {Psi >= delta} outside the target: planar positions uniform in the corridor, height and pitch uniform in
    the ellipse CZ dz^2 + CO dth^2 < rho allowed there, velocities in the ball of radius sqrt(2 mu rho') around k1."""
    r = np.random.default_rng(seed); m = 30 * n
    P = np.array([r.uniform(*BOX[0], m), r.uniform(*BOX[1], m), Z0 * np.ones(m), np.zeros(m)])
    rho = cert.psi_of_y(P) - cert.shift - delta; keep = (rho > 0) & (phi_xy(P) > 0); P, rho = P[:, keep][:, :n], rho[keep][:n]; n = P.shape[1]
    a = np.sqrt(r.random(n)) * 0.999; b = r.uniform(0, 2 * np.pi, n); Y = P.copy(); Y[2] += np.sqrt(rho / CZ) * a * np.cos(b); Y[3] += np.sqrt(rho / CO) * a * np.sin(b)
    rho2 = cert.psi_of_y(Y) - cert.shift - delta; rad = np.sqrt(2 * cert.mu * np.maximum(rho2, 0))
    d = r.normal(size=(4, n)); d /= np.linalg.norm(d, axis=0); sf = r.random(n)**0.25; sf[r.random(n) < boundary_frac] = 1.0
    Yd = cert.k1_of_y(Y) + rad * sf * d; Q, err = ik_batch(Y); ok = err < 1e-9; Q, Yd, Y = Q[:, ok], Yd[:, ok], Y[:, ok]
    Jb = np.array(J_f.map(Q.shape[1], "thread", 8)(Q)).reshape(4, Q.shape[1], 4).transpose(1, 0, 2); QD = np.linalg.solve(Jb, Yd.T[:, :, None])[:, :, 0].T
    return np.vstack([Q, QD]), Y


_xs, _us = ca.SX.sym("x", 2 * NQ), ca.SX.sym("u", NQ)
_qq, _vv = _xs[:NQ], _xs[NQ:]; _M = ca.substitute(ca.substitute(M_sx, Qs, _qq), QDs, _vv); _B = ca.substitute(ca.substitute(B_sx, Qs, _qq), QDs, _vv)
FDYN = ca.Function("f", [_xs, _us], [ca.vertcat(_vv, ca.solve(_M, _us - _B))]); _FM = {}


def rhs_batch(X, U):
    n = X.shape[1]
    if n not in _FM: _FM[n] = FDYN.map(n, "thread", 8)
    return np.array(_FM[n](X, U))


def zoh(X, U, T, n=None):
    """Hold U over [0, T] with RK4 sub-steps of at most 0.5 ms (the MuJoCo step)."""
    n = n or max(1, int(round(T / 0.0005))); h = T / n; X = X.copy()
    for _ in range(n):
        a1 = rhs_batch(X, U); a2 = rhs_batch(X + h / 2 * a1, U); a3 = rhs_batch(X + h / 2 * a2, U); a4 = rhs_batch(X + h * a3, U); X = X + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
    return X
