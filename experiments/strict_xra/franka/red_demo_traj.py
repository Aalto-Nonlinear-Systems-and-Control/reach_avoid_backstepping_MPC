"""Demo trajectory of the red task for the video, as a full pick and place: from rest above the red cube the hand
descends to the cube, grasps it and lifts it back to the transport height (pick), RA-MPC carries it from rest into the
target X_T (transport, the certified part), and the hand moves over the slot, lowers the cube onto the table, releases it
and retracts (place). Pick and place are not certified: they track minimum-jerk task-space references by computed torque,
clipped to the same torque bounds. Everything runs on the MuJoCo plant (model/panda_simple.xml, 0.5 ms RK4, inputs held
over T); the 0.5 ms states are saved with their phase. Writes data/red_fwd_sq/demo.npz."""
import time, numpy as np, mujoco
import config_franka as cf
from cert_franka import ik_batch
from franka_model import MJ, NQ, Y_f, J_f, M_f, B_f, JDQD_f
H_CUBE = 0.02; Z_LOW = H_CUBE; UMAX = np.array(cf.U_MAX); KP, KD = 400.0, 40.0
src = cf.ALPHA * np.array(cf.SPEC["source"]); slot = np.array([cf.TARGET[0], cf.TARGET[1]])
Y = np.array([[src[0]], [src[1]], [cf.Z0], [0.0]]); Q, err = ik_batch(Y); x0 = np.r_[Q[:, 0], np.zeros(NQ)]
print(f"pre-grasp pose q = {np.round(Q[:, 0], 4)} (IK error {err.max():.1e})", flush=True)
pl = mujoco.MjData(MJ); n_sub = int(round(cf.T / MJ.opt.timestep)); QS, PH, US = [x0[:NQ].copy()], [], []


def run(x, u, ph):
    pl.qpos[:] = x[:NQ]; pl.qvel[:] = x[NQ:]; pl.ctrl[:] = u
    for _ in range(n_sub): mujoco.mj_step(MJ, pl); QS.append(pl.qpos.copy()); PH.append(ph); US.append(np.array(u, float))
    return np.r_[pl.qpos.copy(), pl.qvel.copy()]


def task(x):
    q, v = x[:NQ], x[NQ:]; J = np.array(J_f(q)); return np.array(Y_f(q)).ravel(), J @ v, J


def move(x, y_goal, dur, ph, settle=0.6, tol=5e-4):
    """Minimum-jerk reference from the current output to y_goal over dur, tracked by clipped computed torque."""
    y_start = task(x)[0]; sat = 0; k = 0
    while True:
        t = k * cf.T; s = min(t / dur, 1.0); D = y_goal - y_start
        yr = y_start + D * (10 * s**3 - 15 * s**4 + 6 * s**5); vr = D * (30 * s**2 - 60 * s**3 + 30 * s**4) / dur * (t < dur); ar = D * (60 * s - 180 * s**2 + 120 * s**3) / dur**2 * (t < dur)
        y, yd, J = task(x); q, v = x[:NQ], x[NQ:]
        if t >= dur and (np.linalg.norm(y[:3] - y_goal[:3]) < tol and np.linalg.norm(yd[:3]) < 2 * tol or t > dur + settle): break
        a = ar + KP * (yr - y) + KD * (vr - yd); qdd = np.linalg.solve(J, a - np.array(JDQD_f(q, v)).ravel())
        u = np.array(M_f(q)) @ qdd + np.array(B_f(q, v)).ravel(); sat += int(np.any(np.abs(u) > UMAX)); x = run(x, np.clip(u, -UMAX, UMAX), ph); k += 1
    print(f"  phase {ph}: {k * cf.T:.2f} s, {sat} of {k} steps clipped, error {1e3 * np.linalg.norm(y[:3] - y_goal[:3]):.2f} mm", flush=True); return x


def hold(x, dur, ph):
    y = task(x)[0]; return move(x, y, dur, ph, settle=0.0)


t_mark = {}
x = move(x0, np.array([src[0], src[1], Z_LOW, 0.0]), 1.8, 0)                                   # pick: descend to the cube
t_mark["grasp"] = len(PH) * MJ.opt.timestep; x = hold(x, 0.3, 0)
x = move(x, np.array([src[0], src[1], cf.Z0, 0.0]), 1.2, 0, settle=1.0, tol=2e-4)                # lift back to the pre-grasp pose
t_mark["ra_start"] = len(PH) * MJ.opt.timestep
t0 = time.time(); L = cf.mpc_run(x, "cert"); print(f"RA-MPC from the lifted state: {L['status']} in {L['t'][-1]:.2f} s, {len(L['u'])} steps, solve median {np.median(L['solve_ms']):.0f} ms ({time.time()-t0:.0f}s)", flush=True)
assert L["status"] == "reached", "RA-MPC failed"
for u in L["u"]: x = run(x, u, 1)
print(f"  RA-MPC inputs replayed: state mismatch {np.abs(x - np.asarray(L['x'][-1], float)).max():.1e}", flush=True)
t_mark["entry"] = len(PH) * MJ.opt.timestep
x = move(x, np.array([slot[0], slot[1], cf.Z0, 0.0]), 0.6, 2)                                    # place: centre over the slot
x = move(x, np.array([slot[0], slot[1], Z_LOW, 0.0]), 1.8, 2)                                    # lower the cube onto the table
t_mark["release"] = len(PH) * MJ.opt.timestep; x = hold(x, 0.3, 2)
x = move(x, np.array([slot[0], slot[1], 0.18, 0.0]), 0.9, 2, settle=0.3)                          # retract
QS = np.array(QS); TS = np.arange(len(QS)) * MJ.opt.timestep; Ys = np.array([np.array(Y_f(q)).ravel() for q in QS[::4]]); ph4 = np.r_[PH[0], PH][::4]
tr = ph4 == 1; psi = [float(cf.C.F["psi"](np.r_[q, np.zeros(NQ)])) for q in QS[::4][tr]]
print(f"total {TS[-1]:.2f} s; marks {', '.join(f'{k} {v:.2f}' for k, v in t_mark.items())}; transport: min psi {min(psi):.4f}, max |z - z0| {1e3*np.abs(Ys[tr,2]-cf.Z0).max():.1f} mm, "
      f"max |pitch| {np.degrees(np.abs(Ys[tr,3]).max()):.2f} deg; final slot error {1e3*np.hypot(*(Ys[-1,:2]-slot)):.1f} mm; lowest tool point {1e3*Ys[:,2].min():.1f} mm (pick {1e3*Ys[ph4==0,2].min():.1f}, place {1e3*Ys[ph4==2,2].min():.1f})")
np.savez(f"{cf.DATA}/demo.npz", t=TS, q=QS, u=np.array(US), phase=np.array(PH), x0=x0, umax=UMAX, h_cube=H_CUBE, **{f"t_{k}": v for k, v in t_mark.items()})
