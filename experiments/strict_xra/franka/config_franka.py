"""Setting of the Panda kitting task (TASK, default red_fwd_sq): exact time scaling, torque held over T, MuJoCo plant."""
import os, sys, json, numpy as np, mujoco
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from franka_model import MJ, NQ, Y_f, J_f
from cert_franka import CertFranka, sample, phi_xy, ik_batch, TASK, BOX, Z0, TARGET, ALPHA, SPEC
import mpc_franka as mp
SC = json.load(open(os.path.join(HERE, "data", f"scale_{TASK}.json")))
MU0, DFRAC, T = SC["MU0"], SC["DFRAC"], SC["T"]; U_MAX = tuple(SC["U_MAX"]); SCALE = float(os.environ.get("SCALE", SC["s"])); URA_SCALE = float(os.environ.get("URA_SCALE", SC["URA_SCALE"]))
N = int(os.environ.get("N_HOR", 60)); T_MAX_MPC, T_MAX_URA = 6.0, 20.0; SAFE_MARGIN = 0.005; R0 = float(os.environ.get("R0", 0.3)); SEED_POS, N_INIT = 5, 150
DATA, FIGS = os.path.join(HERE, "data", TASK), os.path.join(HERE, "figs", TASK); os.makedirs(DATA, exist_ok=True); os.makedirs(FIGS, exist_ok=True)
C0 = CertFranka(mu=MU0); C = C0.scaled(SCALE); M_PSI = 1.0 - C.shift; DELTA = DFRAC * M_PSI


def positions(n=N_INIT, seed=SEED_POS):
    r = np.random.default_rng(seed); Y = np.array([r.uniform(*BOX[0], 200000), r.uniform(*BOX[1], 200000), Z0 * np.ones(200000), np.zeros(200000)])
    keep = (C.psi_of_y(Y) > 0) & (phi_xy(Y) > 0); return Y[:, keep][:, :n]


def initial_set(n=N_INIT, seed=SEED_POS, r0=R0):
    """Positions uniform in the corridor (height Z0, hand down), planar velocities uniform in a disc of radius
    URA_SCALE * ALPHA * r0 around the URA's field URA_SCALE * k1 (Dubins convention); no height or pitch rate."""
    Y = positions(n, seed); r = np.random.default_rng(seed + 100); m = Y.shape[1]; K = C0.k1_of_y(Y); rad = ALPHA * r0 * np.sqrt(r.random(m)); ang = r.uniform(0, 2 * np.pi, m)
    Yd = np.zeros_like(Y); Yd[:2] = URA_SCALE * (K[:2] + rad * np.array([np.cos(ang), np.sin(ang)])); Q, err = ik_batch(Y); assert err.max() < 1e-9
    Jb = np.array(J_f.map(m, "thread", 8)(Q)).reshape(4, m, 4).transpose(1, 0, 2); QD = np.linalg.solve(Jb, Yd.T[:, :, None])[:, :, 0].T
    return np.vstack([Q, QD])


def ura_run(x0, T_max=T_MAX_URA, scale=URA_SCALE, dense=False, every=1):
    """Continuous terminal controller at the URA scale: u^c recomputed at every MuJoCo step (0.5 ms)."""
    c = C0.scaled(scale); d = mujoco.MjData(MJ); d.qpos[:] = x0[:NQ]; d.qvel[:] = x0[NQ:]; dt = MJ.opt.timestep; st = "timeout"; cost = 0.0; up = np.zeros(NQ); tv = 0.0; tt, xs, us = [], [], []
    for k in range(int(round(T_max / dt))):
        x = np.r_[d.qpos, d.qvel]; y = np.array(Y_f(d.qpos)).ravel()
        if phi_xy(y) <= 0: st = "reached"; break
        if float(c.F["psi"](x)) < 0: st = "unsafe"; break
        u = np.array(c.F["uc"](x)).ravel(); d.ctrl[:] = u; up = np.maximum(up, np.abs(u)); tv += dt * float(np.any(np.abs(u) > np.asarray(U_MAX)))
        cost += dt * (mp.Q_Y * ((y[0] - TARGET[0])**2 + (y[1] - TARGET[1])**2) + mp.R_U * u @ u)
        if dense and k % every == 0: tt.append(k * dt); xs.append(x.copy()); us.append(u.copy())
        mujoco.mj_step(MJ, d)
    out = dict(status=st, t=k * dt, upeak=up, cost=cost, t_viol=tv, psi0=float(c.Psi(np.asarray(x0).reshape(-1, 1))[0]))
    if dense: out.update(tt=np.array(tt), x=np.array(xs), u=np.array(us), y=np.array([np.array(Y_f(z[:NQ])).ravel() for z in xs]))
    return out


def uc_zoh(X0, dt=T, T_max=T_MAX_URA, cert=None):
    """Sampled-data terminal controller alone (held torque) in MuJoCo; status 1 reached, 2 unsafe, 0 timeout."""
    c = C if cert is None else cert; st = np.zeros(X0.shape[1], int); tt = np.full(X0.shape[1], np.nan); dmin = np.full(X0.shape[1], np.inf); pl = mp.Plant()
    for i in range(X0.shape[1]):
        x = X0[:, i].copy(); P = float(c.Psi(x.reshape(-1, 1))[0])
        for k in range(int(round(T_max / dt))):
            u = np.array(c.F["uc"](x)).ravel(); x = pl.step(x, u, dt); y = np.array(Y_f(x[:NQ])).ravel()
            if phi_xy(y) <= 0: st[i], tt[i] = 1, (k + 1) * dt; break
            if float(c.F["psi"](x)) < 0: st[i] = 2; break
            Pn = float(c.F["Psi"](x)); dmin[i] = min(dmin[i], Pn - P); P = Pn
    return st, tt, dmin


_CTRL = {}


def mpc_run(x0, terminal="cert", dt=T, N=N, omega=None, T_max=T_MAX_MPC):
    """One controller object (with its compiled NLP) per process and setting, reused across runs."""
    key = (terminal, dt, N, omega)
    if key not in _CTRL: _CTRL[key] = mp.FrankaMPC(C, dt=dt, N=N, dfrac=DFRAC, omega=omega, safe_margin=SAFE_MARGIN, u_max=U_MAX, nsub=max(1, int(round(dt / 0.005))), terminal=terminal)
    return _CTRL[key].run(np.asarray(x0, float), T_max=T_max)


def safety_until_target(xs, us, dt, nsub=8):
    m = np.inf
    for k in range(len(us)):
        x = np.asarray(xs[k], float).copy()
        for _ in range(nsub):
            x = mp.rk4(x, us[k], dt / nsub); y = np.array(Y_f(x[:NQ])).ravel()
            if phi_xy(y) <= 0: return m, True
            m = min(m, float(C.F["psi"](x)))
    return m, False
