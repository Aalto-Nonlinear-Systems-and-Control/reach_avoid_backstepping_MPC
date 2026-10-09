"""Demo trajectory of the red task for the video: RA-MPC from rest at the pick pose of the red cube (hand down, transport
height), then the sampled-data terminal controller u^c until the tool centre point settles at the slot; the applied
inputs are re-simulated in MuJoCo (model/panda_simple.xml) at 0.5 ms and saved with their phase. Writes data/red_fwd_sq/demo.npz."""
import time, numpy as np, mujoco
import config_franka as cf
from cert_franka import ik_batch, phi_xy
from franka_model import MJ, NQ, Y_f, J_f
src = cf.ALPHA * np.array(cf.SPEC["source"]); Y = np.array([[src[0]], [src[1]], [cf.Z0], [0.0]]); Q, err = ik_batch(Y); x0 = np.r_[Q[:, 0], np.zeros(NQ)]
print(f"pick pose q = {np.round(Q[:, 0], 4)} (IK error {err.max():.1e}); Psi at rest {float(cf.C.Psi(x0.reshape(-1, 1))[0]):+.3f}", flush=True)
t0 = time.time(); L = cf.mpc_run(x0, "cert"); print(f"RA-MPC from rest: {L['status']} in {L['t'][-1]:.2f} s, {len(L['J'])} steps, solve median {np.median(L['solve_ms']):.0f} ms ({time.time()-t0:.0f}s)", flush=True)
assert L["status"] == "reached", "RA-MPC from rest failed"
from franka_model import M_f, B_f, JDQD_f
U = [u for u in L["u"]]; x = np.asarray(L["x"][-1], float); ph = [0] * len(U); pl = cf.mp.Plant(); yslot = np.array([cf.TARGET[0], cf.TARGET[1], cf.Z0, 0.0])
for k in range(100):                                                     # placement after the target is reached: task-space PD, clipped to the bounds
    q, v = x[:NQ], x[NQ:]; y = np.array(Y_f(q)).ravel(); J = np.array(J_f(q)); yd = J @ v
    if np.linalg.norm(y[:2] - yslot[:2]) < 0.001 and np.linalg.norm(yd[:2]) < 0.005: break
    a = 120.0 * (yslot - y) - 22.0 * yd; qdd = np.linalg.solve(J, a - np.array(JDQD_f(q, v)).ravel())
    u = np.clip(np.array(M_f(q)) @ qdd + np.array(B_f(q, v)).ravel(), -np.array(cf.U_MAX), np.array(cf.U_MAX)); U.append(u); ph.append(1); x = pl.step(x, u, cf.T)
print(f"settled after {k*cf.T:.2f} s more; final distance to the slot {1e3*np.hypot(y[0]-cf.TARGET[0], y[1]-cf.TARGET[1]):.1f} mm", flush=True)
d = mujoco.MjData(MJ); d.qpos[:] = x0[:NQ]; d.qvel[:] = x0[NQ:]; n = int(round(cf.T / MJ.opt.timestep)); qs, ts, uu, pp = [d.qpos.copy()], [0.0], [], []
for k, (u, p) in enumerate(zip(U, ph)):
    d.ctrl[:] = u
    for _ in range(n): mujoco.mj_step(MJ, d); qs.append(d.qpos.copy()); ts.append(d.time); uu.append(u.copy()); pp.append(p)
qs = np.array(qs); Ys = np.array([np.array(Y_f(q)).ravel() for q in qs[::4]])
t_entry = L["t"][-1]; print(f"re-simulated {len(qs)} states over {ts[-1]:.2f} s; entry into the target disc at {t_entry:.2f} s; min psi along the path {min(float(cf.C.F['psi'](np.r_[q, np.zeros(NQ)])) for q in qs[::4]):.4f}; max |z - z0| {1e3*np.abs(Ys[:,2]-cf.Z0).max():.1f} mm, max |pitch| {np.degrees(np.abs(Ys[:,3]).max()):.2f} deg")
np.savez(f"{cf.DATA}/demo.npz", t=np.array(ts), q=qs, u=np.array(uu), phase=np.array(pp), t_entry=t_entry, x0=x0, umax=np.array(cf.U_MAX))
