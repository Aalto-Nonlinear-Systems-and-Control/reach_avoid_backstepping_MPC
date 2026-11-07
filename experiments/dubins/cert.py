"""Backstepping certificate for the Dubins car built directly from k1(y).

  Psi(x) = psi(y) - delta_s/lambda - ||ydot - k1(y)||^2 / (2 mu),   ydot = v (cos th, sin th)
  yddot_des = mu grad psi + J_k1 ydot + (lambda/2)(ydot - k1)
  a = yddot_des . t,   omega = yddot_des . n / v,      t = (cos th, sin th), n = (-sin th, cos th)

The terminal set {Psi >= delta} is, in (y, ydot) coordinates, the velocity disc of radius
sqrt(2 mu (psi - delta_s/lambda - delta)) around k1(y): it is sampled directly (no state box).
"""
import numpy as np, sympy as sp

y1, y2 = sp.symbols("y1 y2")
PHI = y2**2 + ((y1 + sp.Rational(17, 10)) * 2)**2 - sp.Rational(2, 5)
H = -(y1**4 + y2**4 - 16) * (y1**4 + y2**4 - 4)


def _p(txt):
    return sp.sympify(txt.replace("^", "**"), locals={"y1": y1, "y2": y2})


class Cert:
    def __init__(self, k1_file=None, mu=0.1, psi=None, k1=None, lam=None, delta_s=None):
        if k1_file is not None:
            t = open(k1_file).read().split("\n")
            k1 = [_p(t[0]), _p(t[1])]; lam, delta_s = float(t[2]), float(t[3])
            if psi is None:
                psi = _p(t[4]) * H / float(t[5])
        self.mu, self.lam, self.ds = mu, lam, delta_s
        self.psi_sym, self.k1_sym = psi, k1
        Y = (y1, y2)
        self.psi = sp.lambdify(Y, psi, "numpy")
        self.gpsi = sp.lambdify(Y, [sp.diff(psi, a) for a in Y], "numpy")
        self.k1 = sp.lambdify(Y, k1, "numpy")
        self.J = sp.lambdify(Y, [[sp.diff(k, a) for a in Y] for k in k1], "numpy")
        self.phi = sp.lambdify(Y, PHI, "numpy")
        self.h = sp.lambdify(Y, H, "numpy")
        self.shift = delta_s / lam

    @staticmethod
    def _b(a, n):
        return np.broadcast_to(np.asarray(a, float), (n,))

    def parts(self, x):
        n = x.shape[1]
        y = x[:2]; c, s = np.cos(x[2]), np.sin(x[2]); yd = np.array([x[3] * c, x[3] * s])
        k = np.array([self._b(a, n) for a in self.k1(*y)])
        g = np.array([self._b(a, n) for a in self.gpsi(*y)])
        Jm = self.J(*y); J = np.array([[self._b(Jm[i][j], n) for j in range(2)] for i in range(2)])
        return y, c, s, yd, k, g, J

    def Psi(self, x):
        y, c, s, yd, k, g, J = self.parts(x)
        return self._b(self.psi(*y), x.shape[1]) - self.shift - ((yd - k)**2).sum(0) / (2 * self.mu)

    def uc(self, x):
        y, c, s, yd, k, g, J = self.parts(x)
        ydd = self.mu * g + np.einsum("ijn,jn->in", J, yd) + 0.5 * self.lam * (yd - k)
        return np.array([(-s * ydd[0] + c * ydd[1]) / x[3], c * ydd[0] + s * ydd[1]])  # (omega, a)

    def Psidot(self, x, u):
        """dPsi/dt along the plant with input u = (omega, a) (frozen or not)."""
        y, c, s, yd, k, g, J = self.parts(x)
        ydd = np.array([u[1] * c - x[3] * u[0] * s, u[1] * s + x[3] * u[0] * c])
        e = yd - k
        return (g * yd).sum(0) - (e * (ydd - np.einsum("ijn,jn->in", J, yd))).sum(0) / self.mu

    def scaled(self, s):
        """Time-rescaled certificate: (k1, lambda, delta_s) -> s (k1, lambda, delta_s), mu -> s^2 mu.
        The closed loop is y_s(t) = y(s t); the certified sampling period scales as 1/s."""
        return Cert(mu=self.mu * s * s, psi=self.psi_sym, k1=[s * k for k in self.k1_sym],
                    lam=self.lam * s, delta_s=self.ds * s)

    def psi_state_sym(self):
        x1, x2, th, v = sp.symbols("x1 x2 th v")
        sub = {y1: x1, y2: x2}
        k = [e.subs(sub, simultaneous=True) for e in self.k1_sym]
        return (self.psi_sym.subs(sub, simultaneous=True) - self.shift
                - ((v * sp.cos(th) - k[0])**2 + (v * sp.sin(th) - k[1])**2) / (2 * self.mu)), (x1, x2, th, v)

    # ── sampling of {Psi >= delta} \ X_T in (y, ydot) coordinates ────────────
    def sample(self, delta, n_y=200000, seed=0, boundary_frac=0.5):
        r = np.random.default_rng(seed)
        y = r.uniform(-2.05, 2.05, (2, n_y))
        rho = self._b(self.psi(*y), n_y) - self.shift - delta
        keep = (rho > 0) & (self._b(self.phi(*y), n_y) > 0) & (self._b(self.h(*y), n_y) > 0)
        y = y[:, keep]; rad = np.sqrt(2 * self.mu * rho[keep]); n = y.shape[1]
        k = np.array([self._b(a, n) for a in self.k1(*y)])
        sfrac = np.sqrt(r.random(n)); sfrac[r.random(n) < boundary_frac] = 1.0
        ang = r.uniform(0, 2 * np.pi, n)
        yd = k + rad * sfrac * np.array([np.cos(ang), np.sin(ang)])
        v = np.hypot(yd[0], yd[1]); th = np.arctan2(yd[1], yd[0])
        x = np.array([y[0], y[1], th, v])
        info = dict(n_y=n, w=rad**2, kn=np.hypot(k[0], k[1]), rad=rad, y=y)
        return x, info

    def zoh(self, x, u, T, n=40, track=False):
        """RK4 integration of the Dubins plant under the frozen input u; optionally track Psi and Psidot."""
        h = T / n; x = x.copy()
        f = lambda z: np.array([z[3] * np.cos(z[2]), z[3] * np.sin(z[2]), u[0], u[1]])
        if track:
            pd_prev = self.Psidot(x, u); M2 = np.zeros(x.shape[1]); pmin = self.Psi(x)
        for _ in range(n):
            a1 = f(x); a2 = f(x + h / 2 * a1); a3 = f(x + h / 2 * a2); a4 = f(x + h * a3)
            x = x + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            if track:
                pd = self.Psidot(x, u); M2 = np.maximum(M2, np.abs(pd - pd_prev) / h); pd_prev = pd
                pmin = np.minimum(pmin, self.Psi(x))
        return (x, M2, pmin) if track else x

    def margin_on_grid(self, n=1201):
        g = np.linspace(-2.05, 2.05, n); A, B = np.meshgrid(g, g)
        P = self.psi(A, B); F = self.phi(A, B); Hh = self.h(A, B); reg = (Hh >= 0) & (F >= 0)
        gp = self.gpsi(A, B); k = self.k1(A, B)
        m = gp[0] * k[0] + gp[1] * k[1] - self.lam * P + self.ds
        kmax = max(np.abs(np.broadcast_to(k[0], A.shape)[Hh >= 0]).max(), np.abs(np.broadcast_to(k[1], A.shape)[Hh >= 0]).max())
        out = reg & (F > 0)
        return dict(min_margin=float(m[reg].min()), kmax=float(kmax),
                    area_share=float(((P >= self.shift) & out).sum() / out.sum()))
