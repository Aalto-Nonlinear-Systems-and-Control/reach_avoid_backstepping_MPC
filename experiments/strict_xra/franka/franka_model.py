"""Simplified Panda-like arm (model/panda_simple.xml): Panda kinematics with q3 = q5 = 0, q7 = -pi/4, active joints
(q1, q2, q4, q6), point masses at the elbow (M_E), the wrist (M_W) and the tool centre point (M_H), yaw inertia I1,
armature 0.1, no gravity, no damping. Closed-form kinematics from the Panda DH parameters (modified DH), dynamics
M(q) = sum m_i J_i^T J_i + diag(I1, 0, 0, 0) + 0.1 I, b = Mdot qd - 1/2 d/dq (qd^T M qd). Outputs y = (-Y, X, Z, pitch) of
the tool centre point, pitch = -q2 + q4 + q6 (hand pointing down at 0)."""
import os, numpy as np, sympy as sp, casadi as ca, mujoco
HERE = os.path.dirname(os.path.abspath(__file__)); XML = os.path.join(HERE, "model", "panda_simple.xml")
MJ = mujoco.MjModel.from_xml_path(XML); NQ = MJ.nq
M_E, M_W, M_H, I1, ARM = 6.8, 2.9, 1.5, 0.05, 0.1
DH = [(0, 0, 0.333), (0, -sp.pi / 2, 0), (0, sp.pi / 2, 0.316), (0.0825, sp.pi / 2, 0), (-0.0825, -sp.pi / 2, 0.384), (0, sp.pi / 2, 0), (0.088, sp.pi / 2, 0)]
s1, s2, s4, s6 = sp.symbols("q1 q2 q4 q6", real=True); TH = [s1, s2, 0, s4, 0, s6, -sp.pi / 4]


def _T(a, al, d, th):
    ca_, sa, ct, st = sp.cos(al), sp.sin(al), sp.cos(th), sp.sin(th)
    return sp.Matrix([[ct, -st, 0, a], [st * ca_, ct * ca_, -sa, -sa * d], [st * sa, ct * sa, ca_, ca_ * d], [0, 0, 0, 1]])


_Ts = [sp.eye(4)]
for (a, al, d), th in zip(DH, TH): _Ts.append(sp.simplify(_Ts[-1] * _T(a, al, d, th)))
P_ELB = sp.simplify(_Ts[4][:3, 3]); P_WRI = sp.simplify(_Ts[6][:3, 3]); P_TCP = sp.simplify(_Ts[7][:3, 3] + _Ts[7][:3, 2] * (0.107 + 0.1034))


def s2c(e, mp):
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


q = ca.SX.sym("q", NQ); qd = ca.SX.sym("qd", NQ); _mp = {s1: q[0], s2: q[1], s4: q[2], s6: q[3]}
ELB = ca.vertcat(*[s2c(e, _mp) for e in P_ELB]); WRI = ca.vertcat(*[s2c(e, _mp) for e in P_WRI]); TCP_sx = ca.vertcat(*[s2c(e, _mp) for e in P_TCP])
M_sx = ARM * ca.SX.eye(NQ); M_sx[0, 0] += I1
for m, p in ((M_E, ELB), (M_W, WRI), (M_H, TCP_sx)): Jv = ca.jacobian(p, q); M_sx += m * Jv.T @ Jv
B_sx = ca.jacobian(M_sx @ qd, q) @ qd - 0.5 * ca.gradient(qd.T @ M_sx @ qd, q)
PITCH_C = [-1.0, 1.0, 1.0, 0.0]
Y_sx = ca.vertcat(-TCP_sx[1], TCP_sx[0], TCP_sx[2], -q[1] + q[2] + q[3]); J_sx = ca.jacobian(Y_sx, q); JDQD_sx = ca.jacobian(J_sx @ qd, q) @ qd
M_f = ca.Function("M", [q], [M_sx]); B_f = ca.Function("b", [q, qd], [B_sx]); TCP_f = ca.Function("tcp", [q], [TCP_sx]); Y_f = ca.Function("y", [q], [Y_sx]); J_f = ca.Function("J", [q], [J_sx]); JDQD_f = ca.Function("jdqd", [q, qd], [JDQD_sx])
_SH = ca.DM([0, 0, 0.333]); JOINTS_f = ca.Function("joints", [q], [ca.horzcat(ca.DM.zeros(3), _SH, ELB, WRI, TCP_sx)])


def fwd(x, tau):
    qq, vv = x[:NQ], x[NQ:]; return np.r_[vv, np.linalg.solve(np.array(M_f(qq)), np.asarray(tau) - np.array(B_f(qq, vv)).ravel())]


def ik(y, q0):
    qq = np.array(q0, float)
    for _ in range(60):
        e = np.array(Y_f(qq)).ravel() - y
        if np.linalg.norm(e) < 1e-12: break
        qq = qq - np.linalg.solve(np.array(J_f(qq)), e)
    return qq, float(np.linalg.norm(np.array(Y_f(qq)).ravel() - y))


if __name__ == "__main__":
    print("closed-form points: elbow", [sp.count_ops(e) for e in P_ELB], "wrist", [sp.count_ops(e) for e in P_WRI], "tcp", [sp.count_ops(e) for e in P_TCP], "operations")
    d = mujoco.MjData(MJ); r = np.random.default_rng(1); eM = eB = eT = eP = 0.0
    for _ in range(50):
        qq = np.array([r.uniform(-2, 2), r.uniform(-1.5, 1.5), r.uniform(-3.0, -0.2), r.uniform(0.1, 3.6)]); vv = r.normal(0, 1.5, 4)
        d.qpos[:] = qq; d.qvel[:] = vv; mujoco.mj_forward(MJ, d); Mm = np.zeros((NQ, NQ)); mujoco.mj_fullM(MJ, d, Mm)
        eM = max(eM, np.abs(Mm - np.array(M_f(qq))).max()); eB = max(eB, np.abs(d.qfrc_bias - np.array(B_f(qq, vv)).ravel()).max()); eT = max(eT, np.abs(d.site("tcp").xpos - np.array(TCP_f(qq)).ravel()).max())
        a = d.site("tcp").xmat.reshape(3, 3)[:, 2]; ang = np.arctan2(a[0] * np.cos(qq[0]) + a[1] * np.sin(qq[0]), -a[2]); eP = max(eP, abs((ang - (-qq[1] + qq[2] + qq[3]) + np.pi) % (2 * np.pi) - np.pi))
        for nm, P in (("link4", ELB), ("link6", WRI)): eT = max(eT, np.abs(d.body(nm).xpos - np.array(ca.Function("p", [q], [P])(qq)).ravel()).max())
    print(f"against MuJoCo (panda_simple.xml), 50 random states: |M| {eM:.1e}, |bias| {eB:.1e}, |positions| {eT:.1e}, |pitch| {eP:.1e}")
    qq = np.array([0.3, -0.4, -2.2, 1.9]); vv = np.array([0.5, -0.3, 0.8, -0.4]); tau = np.array([3.0, -2.0, 1.5, 0.3]); x = np.r_[qq, vv]
    d.qpos[:] = qq; d.qvel[:] = vv; d.ctrl[:] = tau
    for _ in range(20): mujoco.mj_step(MJ, d)
    h = 0.0005
    for _ in range(20):
        k1 = fwd(x, tau); k2 = fwd(x + h / 2 * k1, tau); k3 = fwd(x + h / 2 * k2, tau); k4 = fwd(x + h * k3, tau); x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    print(f"10 ms under a held torque: MuJoCo vs closed-form RK4, max |q| diff {np.abs(d.qpos - x[:4]).max():.1e}, |qd| diff {np.abs(d.qvel - x[4:]).max():.1e}")
