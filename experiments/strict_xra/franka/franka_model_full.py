"""Symbolic (CasADi) dynamics and outputs of the reduced Panda (joints 1, 2, 4, 6), built from the compiled MuJoCo model
model/panda_task.xml so that the controller model and the simulator agree. Outputs y = (y1, y2, y3, y4):
y1 = -Y_world (to the right), y2 = X_world (forward), y3 = Z_world of the tool centre point, y4 = pitch of the hand from
the downward vertical in the arm plane (linear in q2, q4, q6). Plant: M(q) qdd + b(q, qd) = tau (no gravity, no damping)."""
import os, numpy as np, casadi as ca, mujoco
HERE = os.path.dirname(os.path.abspath(__file__)); XML = os.path.join(HERE, "model", "panda_task.xml")
MJ = mujoco.MjModel.from_xml_path(XML); NQ = MJ.nq


def _qmat(qt):
    w, x, y, z = np.asarray(qt, float) / np.linalg.norm(qt)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def _rot(axis, a):
    k = ca.DM(axis); K = ca.vertcat(ca.horzcat(0, -k[2], k[1]), ca.horzcat(k[2], 0, -k[0]), ca.horzcat(-k[1], k[0], 0))
    return ca.DM.eye(3) + ca.sin(a) * K + (1 - ca.cos(a)) * (K @ K)


def build():
    q = ca.SX.sym("q", NQ); qd = ca.SX.sym("qd", NQ); R = {0: ca.SX.eye(3)}; P = {0: ca.SX.zeros(3)}; anc = {0: []}; zax = {}
    assert np.allclose(MJ.jnt_pos, 0), "joint anchors are assumed at the body origins"
    for b in range(1, MJ.nbody):
        p = MJ.body_parentid[b]; Rb = R[p] @ ca.DM(_qmat(MJ.body_quat[b])); Pb = P[p] + R[p] @ ca.DM(MJ.body_pos[b]); A = list(anc[p])
        for j in [j for j in range(MJ.njnt) if MJ.jnt_bodyid[j] == b]:
            i = MJ.jnt_dofadr[j]; zax[i] = Rb @ ca.DM(MJ.jnt_axis[j]); Rb = Rb @ _rot(MJ.jnt_axis[j], q[i]); A.append(i)
        R[b], P[b], anc[b] = Rb, Pb, A
    M = ca.SX.zeros(NQ, NQ)
    for b in range(1, MJ.nbody):
        if MJ.body_mass[b] <= 0: continue
        c = P[b] + R[b] @ ca.DM(MJ.body_ipos[b]); Jv = ca.jacobian(c, q)
        Jw = ca.horzcat(*[zax[i] if i in anc[b] else ca.SX.zeros(3) for i in range(NQ)])
        Ri = R[b] @ ca.DM(_qmat(MJ.body_iquat[b])); Iw = Ri @ ca.diag(ca.DM(MJ.body_inertia[b])) @ Ri.T
        M += MJ.body_mass[b] * Jv.T @ Jv + Jw.T @ Iw @ Jw
    M += ca.diag(ca.DM(MJ.dof_armature))
    bias = ca.jacobian(M @ qd, q) @ qd - 0.5 * ca.gradient(qd.T @ M @ qd, q)
    s = MJ.site("tcp").id; hb = MJ.site_bodyid[s]; tcp = P[hb] + R[hb] @ ca.DM(MJ.site_pos[s]); app = R[hb][:, 2]
    return q, qd, M, bias, tcp, app, R, P


q, qd, M_sx, B_sx, TCP_sx, APP_sx, _R, _P = build()
M_f = ca.Function("M", [q], [M_sx]); B_f = ca.Function("b", [q, qd], [B_sx]); TCP_f = ca.Function("tcp", [q], [TCP_sx]); APP_f = ca.Function("app", [q], [APP_sx])
JOINTS_f = ca.Function("joints", [q], [ca.horzcat(*[_P[MJ.body(n).id] for n in ("link1", "link2", "link4", "link6", "link7", "hand")])])

# pitch of the hand from the downward vertical in the arm plane: theta = c . (q2, q4, q6) + c0, exact for q3 = q5 = 0
def _theta(qq):
    a = np.array(APP_f(qq)).ravel(); return np.arctan2(a[0], -a[2])        # q1 = 0: the arm plane is x-z
_qr = np.array([0.0, -0.3, -2.0, 1.8]); _c = [round((_theta(_qr + 1e-6 * np.eye(4)[i]) - _theta(_qr - 1e-6 * np.eye(4)[i])) / 2e-6) for i in (1, 2, 3)]
_c0 = _theta(_qr) - (_c[0] * _qr[1] + _c[1] * _qr[2] + _c[2] * _qr[3]); _c0 = (_c0 + np.pi) % (2 * np.pi) - np.pi
PITCH_C = [float(_c[0]), float(_c[1]), float(_c[2]), float(_c0)]
_r = np.random.default_rng(0); _wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi
PITCH_FIT_ERR = max(abs(_wrap(_theta(qq) - (PITCH_C[0] * qq[1] + PITCH_C[1] * qq[2] + PITCH_C[2] * qq[3] + PITCH_C[3])))
                    for qq in np.c_[np.zeros(200), _r.uniform(-1.5, 1.5, 200), _r.uniform(-3.0, -0.2, 200), _r.uniform(0.0, 3.7, 200)])
Y_sx = ca.vertcat(-TCP_sx[1], TCP_sx[0], TCP_sx[2], PITCH_C[0] * q[1] + PITCH_C[1] * q[2] + PITCH_C[2] * q[3] + PITCH_C[3])
Y_f = ca.Function("y", [q], [Y_sx]); J_sx = ca.jacobian(Y_sx, q); J_f = ca.Function("J", [q], [J_sx])
JDQD_sx = ca.jacobian(J_sx @ qd, q) @ qd; JDQD_f = ca.Function("jdqd", [q, qd], [JDQD_sx])


def fwd(x, tau):
    """Plant right-hand side, x = (q, qd) (8,), tau (4,)."""
    qq, vv = x[:NQ], x[NQ:]; return np.r_[vv, np.linalg.solve(np.array(M_f(qq)), np.asarray(tau) - np.array(B_f(qq, vv)).ravel())]


def ik(y, q0):
    """Newton on h(q) = y from q0 (4,)."""
    qq = np.array(q0, float)
    for _ in range(60):
        e = np.array(Y_f(qq)).ravel() - y
        if np.linalg.norm(e) < 1e-12: break
        qq = qq - np.linalg.solve(np.array(J_f(qq)), e)
    return qq, float(np.linalg.norm(np.array(Y_f(qq)).ravel() - y))


if __name__ == "__main__":
    d = mujoco.MjData(MJ); r = np.random.default_rng(1); errM = errB = errT = 0.0
    for _ in range(50):
        qq = np.array([r.uniform(-2, 2), r.uniform(-1.5, 1.5), r.uniform(-3.0, -0.2), r.uniform(0.1, 3.6)]); vv = r.normal(0, 1.5, 4)
        d.qpos[:] = qq; d.qvel[:] = vv; mujoco.mj_forward(MJ, d); Mm = np.zeros((NQ, NQ)); mujoco.mj_fullM(MJ, d, Mm)
        errM = max(errM, np.abs(Mm - np.array(M_f(qq))).max()); errB = max(errB, np.abs(d.qfrc_bias - np.array(B_f(qq, vv)).ravel()).max())
        errT = max(errT, np.abs(d.site("tcp").xpos - np.array(TCP_f(qq)).ravel()).max())
    print(f"model check against MuJoCo on 50 random states: max |M| error {errM:.2e}, max |bias| error {errB:.2e}, max |tcp| error {errT:.2e}")
    print(f"pitch = {PITCH_C[0]:+g} q2 {PITCH_C[1]:+g} q4 {PITCH_C[2]:+g} q6 {PITCH_C[3]:+.6f}, largest deviation on 200 random poses {PITCH_FIT_ERR:.1e}")
    # one MuJoCo step against RK4 on the symbolic model
    qq = np.array([0.3, -0.4, -2.2, 1.9]); vv = np.array([0.5, -0.3, 0.8, -0.4]); tau = np.array([3.0, -2.0, 1.5, 0.3]); x = np.r_[qq, vv]
    d.qpos[:] = qq; d.qvel[:] = vv; d.ctrl[:] = tau
    for _ in range(20): mujoco.mj_step(MJ, d)
    h = 0.0005
    for _ in range(20):
        k1 = fwd(x, tau); k2 = fwd(x + h / 2 * k1, tau); k3 = fwd(x + h / 2 * k2, tau); k4 = fwd(x + h * k3, tau); x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    print(f"10 ms under a held torque: MuJoCo vs symbolic RK4, max |q| diff {np.abs(d.qpos - x[:4]).max():.2e}, |qd| diff {np.abs(d.qvel - x[4:]).max():.2e}")
