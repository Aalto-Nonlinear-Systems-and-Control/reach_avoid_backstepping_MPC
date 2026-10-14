"""Top-view two-link arm (horizontal plane, no gravity): dynamics, backstepping reach-avoid controller, ECGBF certificate,
time scaling, sampling of the certified sets and ZOH simulation. Copy of ../arm/cert_arm.py with GR = 0 and the task
(psi, target disc, sampling box) read from tasks/<TASK>.json, TASK from the environment (default red_fwd).
State x = (q1, q2, dq1, dq2), output y = end-effector position, input tau = (tau1, tau2)."""
import os, json, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
M1, M2, L1, L2, LC1, LC2, I1, I2, GR = 1.0, 1.0, 4.0, 4.0, 2.0, 2.0, 0.02, 0.02, 0.0
y1, y2 = sp.symbols("y1 y2")
TASK = os.environ.get("TASK", "red_fwd"); SPEC = json.load(open(os.path.join(HERE, "tasks", f"{TASK}.json")))
TARGET_POS = np.array(SPEC["target"], float); RT = float(SPEC["RT"]); BOX = tuple(tuple(b) for b in SPEC["box"])
PHI = (y1 - sp.nsimplify(TARGET_POS[0]))**2 + (y2 - sp.nsimplify(TARGET_POS[1]))**2 - sp.nsimplify(RT)**2
phi_np = sp.lambdify((y1, y2), PHI, "numpy")


def _p(txt):
    return sp.sympify(txt.replace("^", "**"), locals={"y1": y1, "y2": y2})


# ── kinematics and dynamics, vectorised over the columns of x (4, n) ──────────────────────────────────────────
def fk(q):
    q1, q2 = q[0], q[1]
    return np.array([L1 * np.cos(q1) + L2 * np.cos(q1 + q2), L1 * np.sin(q1) + L2 * np.sin(q1 + q2)])


def ik(y, elbow="down"):
    """Inverse kinematics, elbow-down branch (q2 >= 0) by default. y: (2, n)."""
    r2 = y[0]**2 + y[1]**2; c2 = np.clip((r2 - L1**2 - L2**2) / (2 * L1 * L2), -1, 1)
    q2 = np.arccos(c2) if elbow == "down" else -np.arccos(c2)
    q1 = np.arctan2(y[1], y[0]) - np.arctan2(L2 * np.sin(q2), L1 + L2 * np.cos(q2))
    return np.array([q1, q2])


def jac(q):
    q1, q2 = q[0], q[1]; s1, c1, s12, c12 = np.sin(q1), np.cos(q1), np.sin(q1 + q2), np.cos(q1 + q2)
    return np.array([[-L1 * s1 - L2 * s12, -L2 * s12], [L1 * c1 + L2 * c12, L2 * c12]])          # (2, 2, n)


def jacdot_qdot(q, dq):
    """(d/dt J) dq, i.e. the velocity-quadratic part of ydd."""
    q1, q2 = q[0], q[1]; d1, d2 = dq[0], dq[1]; s1, c1, s12, c12 = np.sin(q1), np.cos(q1), np.sin(q1 + q2), np.cos(q1 + q2)
    a = -L1 * c1 * d1**2 - L2 * c12 * (d1 + d2)**2
    b = -L1 * s1 * d1**2 - L2 * s12 * (d1 + d2)**2
    return np.array([a, b])


def mass(q):
    c2 = np.cos(q[1]); M11 = I1 + I2 + M1 * LC1**2 + M2 * (L1**2 + LC2**2 + 2 * L1 * LC2 * c2); M12 = M2 * (LC2**2 + L1 * LC2 * c2) + I2
    M22 = (M2 * LC2**2 + I2) * np.ones_like(c2)
    return np.array([[M11, M12], [M12, M22]])


def coriolis_qdot(q, dq):
    s2 = np.sin(q[1]); h = M2 * L1 * LC2 * s2
    return np.array([-h * dq[1] * dq[0] - h * (dq[0] + dq[1]) * dq[1], h * dq[0] * dq[0]])


def gravity(q):
    return np.array([(M1 * GR * LC1 + M2 * GR * L1) * np.cos(q[0]) + M2 * GR * LC2 * np.cos(q[0] + q[1]), M2 * GR * LC2 * np.cos(q[0] + q[1])])


def _solve2(A, b):
    """Solve A z = b for (2,2,n) A and (2,n) b."""
    det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    return np.array([(A[1, 1] * b[0] - A[0, 1] * b[1]) / det, (-A[1, 0] * b[0] + A[0, 0] * b[1]) / det])


def _mul2(A, v):
    return np.array([A[0, 0] * v[0] + A[0, 1] * v[1], A[1, 0] * v[0] + A[1, 1] * v[1]])


def f_np(x, tau):
    """Plant: x = (q, dq), tau (2, n) or (2,)."""
    q, dq = x[:2], x[2:]; tau = np.asarray(tau, float)
    if tau.ndim == 1: tau = tau[:, None] * np.ones((1,) + q.shape[1:]) if q.ndim > 1 else tau
    qdd = _solve2(mass(q), tau - coriolis_qdot(q, dq) - gravity(q))
    return np.concatenate([dq, qdd])


class CertArm:
    """Backstepping ECGBF certificate of the arm for a given k1 design (psi_raw/pmax, k1, lambda, delta_s) and mu."""
    def __init__(self, k1_file=None, mu=0.1, psi=None, k1=None, lam=None, delta_s=None):
        if k1_file is not None:
            t = open(k1_file).read().split("\n"); k1 = [_p(t[0]), _p(t[1])]; lam, delta_s = float(t[2]), float(t[3])
            if psi is None:
                psi = _p(t[4]) / float(t[5])
        self.mu, self.lam, self.ds, self.psi_sym, self.k1_sym = mu, lam, delta_s, psi, k1
        Y = (y1, y2)
        self.psi = sp.lambdify(Y, psi, "numpy"); self.gpsi = sp.lambdify(Y, [sp.diff(psi, a) for a in Y], "numpy")
        self.k1 = sp.lambdify(Y, k1, "numpy"); self.J1 = sp.lambdify(Y, [[sp.diff(k, a) for a in Y] for k in k1], "numpy")
        self.phi = phi_np; self.shift = delta_s / lam

    @staticmethod
    def _b(a, n):
        return np.broadcast_to(np.asarray(a, float), (n,))

    def parts(self, x):
        n = x.shape[1]; q, dq = x[:2], x[2:]; y = fk(q); J = jac(q); yd = _mul2(J, dq)
        k = np.array([self._b(a, n) for a in self.k1(*y)]); g = np.array([self._b(a, n) for a in self.gpsi(*y)])
        Jm = self.J1(*y); J1 = np.array([[self._b(Jm[i][j], n) for j in range(2)] for i in range(2)])
        return y, J, yd, k, g, J1

    def Psi(self, x):
        y, J, yd, k, g, J1 = self.parts(x)
        return self._b(self.psi(*y), x.shape[1]) - self.shift - ((yd - k)**2).sum(0) / (2 * self.mu)

    def ydd_des(self, x):
        y, J, yd, k, g, J1 = self.parts(x)
        return self.mu * g + _mul2(J1, yd) + 0.5 * self.lam * (yd - k)

    def uc(self, x, with_gravity=True):
        """tau = M J^{-1} (ydd_des - Jdot dq) + C dq + G; with_gravity=False returns the part that scales with s^2."""
        q, dq = x[:2], x[2:]; J = jac(q)
        qdd = _solve2(J, self.ydd_des(x) - jacdot_qdot(q, dq))
        tau = _mul2(mass(q), qdd) + coriolis_qdot(q, dq)
        return tau + gravity(q) if with_gravity else tau

    def Psidot(self, x, tau):
        y, J, yd, k, g, J1 = self.parts(x); q, dq = x[:2], x[2:]
        qdd = _solve2(mass(q), np.asarray(tau) - coriolis_qdot(q, dq) - gravity(q)); ydd = _mul2(J, qdd) + jacdot_qdot(q, dq)
        return (g * yd).sum(0) - ((yd - k) * (ydd - _mul2(J1, yd))).sum(0) / self.mu

    def scaled(self, s):
        return CertArm(mu=self.mu * s * s, psi=self.psi_sym, k1=[s * k for k in self.k1_sym], lam=self.lam * s, delta_s=self.ds * s)

    def psi_state_sym(self):
        """Psi as a sympy expression in the state (x1, x2, x3, x4) for the MPC."""
        x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4")
        h1 = L1 * sp.cos(x1) + L2 * sp.cos(x1 + x2); h2 = L1 * sp.sin(x1) + L2 * sp.sin(x1 + x2)
        yd1 = (-L1 * sp.sin(x1) - L2 * sp.sin(x1 + x2)) * x3 - L2 * sp.sin(x1 + x2) * x4
        yd2 = (L1 * sp.cos(x1) + L2 * sp.cos(x1 + x2)) * x3 + L2 * sp.cos(x1 + x2) * x4
        sub = {y1: h1, y2: h2}; k = [e.subs(sub, simultaneous=True) for e in self.k1_sym]
        return (self.psi_sym.subs(sub, simultaneous=True) - self.shift - ((yd1 - k[0])**2 + (yd2 - k[1])**2) / (2 * self.mu)), (x1, x2, x3, x4)

    # ── sampling of {Psi >= delta} \ X_T in (y, ydot), mapped to (q, dq) through the elbow-down branch ──────────
    def sample(self, delta, n_y=200000, seed=0, boundary_frac=0.5, box=None):
        r = np.random.default_rng(seed); box = BOX if box is None else box
        y = np.array([r.uniform(*box[0], n_y), r.uniform(*box[1], n_y)])
        rho = self._b(self.psi(*y), n_y) - self.shift - delta
        keep = (rho > 0) & (self._b(self.phi(*y), n_y) > 0)
        y = y[:, keep]; rad = np.sqrt(2 * self.mu * rho[keep]); n = y.shape[1]
        k = np.array([self._b(a, n) for a in self.k1(*y)])
        sfrac = np.sqrt(r.random(n)); sfrac[r.random(n) < boundary_frac] = 1.0
        ang = r.uniform(0, 2 * np.pi, n)
        yd = k + rad * sfrac * np.array([np.cos(ang), np.sin(ang)])
        q = ik(y); dq = _solve2(jac(q), yd)
        x = np.vstack([q, dq])
        return x, dict(n_y=n, w=rad**2, kn=np.hypot(k[0], k[1]), rad=rad, y=y)

    def zoh(self, x, tau, T, n=40, track=False, comp_gravity=False):
        """Hold tau over [0, T]. With comp_gravity the held quantity is v = tau - G(q_k) and the plant receives
        v + G(q(t)) (gravity compensated continuously by an inner loop)."""
        h = T / n; x = x.copy(); v = np.asarray(tau) - gravity(x[:2]) if comp_gravity else None
        inp = lambda z: v + gravity(z[:2]) if comp_gravity else tau
        if track:
            pd_prev = self.Psidot(x, inp(x)); M2 = np.zeros(x.shape[1]); pmin = self.Psi(x)
        for _ in range(n):
            a1 = f_np(x, inp(x)); a2 = f_np(x + h / 2 * a1, inp(x + h / 2 * a1)); a3 = f_np(x + h / 2 * a2, inp(x + h / 2 * a2)); a4 = f_np(x + h * a3, inp(x + h * a3))
            x = x + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            if track:
                pd = self.Psidot(x, inp(x)); M2 = np.maximum(M2, np.abs(pd - pd_prev) / h); pd_prev = pd; pmin = np.minimum(pmin, self.Psi(x))
        return (x, M2, pmin) if track else x

    def margin_on_grid(self, n=801):
        """ECGBF margin grad psi . k1 - lambda psi + delta_s on the safe set, the largest |k1_i|, and the set's extent."""
        g1 = np.linspace(*BOX[0], n); g2 = np.linspace(*BOX[1], n); A, B = np.meshgrid(g1, g2)
        P = self.psi(A, B); F = self.phi(A, B); safe = P >= 0
        gp = self.gpsi(A, B); k = self.k1(A, B)
        m = gp[0] * k[0] + gp[1] * k[1] - self.lam * P + self.ds
        kmax = max(np.abs(np.broadcast_to(k[0], A.shape)[safe]).max(), np.abs(np.broadcast_to(k[1], A.shape)[safe]).max())
        out = safe & (F > 0)
        return dict(min_margin=float(m[safe & (F > 0)].min()), kmax=float(kmax), area_share=float(((P >= self.shift) & out).sum() / out.sum()),
                    rmax=float(np.hypot(A, B)[safe].max()), rmin=float(np.hypot(A, B)[safe].min()))
