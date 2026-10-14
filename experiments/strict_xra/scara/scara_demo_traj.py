"""Demo of the red task on the SCARA, as a full pick and place. Pick: the arm rests above the red cube while the Z axis
lowers the gripper to the cube, the fingers close and the Z axis lifts it to the transport height (the planar joints do
not move, u = 0). Transport: RA-MPC of ../arm_top from rest into the target X_T (certified). Place: the tool point is
centred over the slot (minimum-jerk reference, clipped computed torque, not certified), the Z axis lowers the cube,
the fingers open and the Z axis retracts while the arm holds still. The planar part runs on the paper's plant (task
units, held torque over T = 10 ms, RK4 at 0.5 ms); the same inputs times C_TAU are re-simulated on model/scara_arm.xml
and compared. The Z axis, wrist and fingers are scripted (they do not enter the planar dynamics). Writes data/demo.npz."""
import time, numpy as np, mujoco
from scara_common import *
import config_top as cf
from cert_top import fk, ik, jac, jacdot_qdot, mass, coriolis_qdot, phi_np, SPEC, TARGET_POS
from ra_mpc_top import rk4_np
T = cf.T; UMAX = np.array(cf.U_MAX); NSUB = int(round(T / 0.0005)); KP, KD = 400.0, 40.0; D_DOWN, D_GRASP, D_UP = 1.0, 0.3, 1.0
src, slot = np.array(SPEC["source"], float), TARGET_POS; q0 = ik(src.reshape(2, 1))[:, 0]; x0 = np.r_[q0, 0.0, 0.0]
X, U, PH = [x0.copy()], [], []


def step(x, u, ph):
    for _ in range(NSUB): x = rk4_np(x, u, T / NSUB)
    X.append(x.copy()); U.append(np.array(u, float)); PH.append(ph); return x


def move(x, y_goal, dur, ph, settle=1.0, tol=1e-3):
    """Planar minimum-jerk move of the tool point (task units) by clipped computed torque; dur = 0 holds y_goal."""
    y0 = fk(x[:2]); k = 0; sat = 0
    while True:
        t = k * T; s = min(t / dur, 1.0) if dur > 0 else 1.0; D = y_goal - y0; on = dur > 0 and t < dur
        yr = y0 + D * minjerk(s); vr = D * (30 * s**2 - 60 * s**3 + 30 * s**4) / dur * on if dur > 0 else 0 * D; ar = D * (60 * s - 180 * s**2 + 120 * s**3) / dur**2 * on if dur > 0 else 0 * D
        q, dq = x[:2], x[2:]; y = fk(q); J = jac(q); yd = J @ dq
        if t >= dur and (np.linalg.norm(y - y_goal) < tol and np.linalg.norm(yd) < 2 * tol or t >= dur + settle): break
        qdd = np.linalg.solve(J, ar + KP * (yr - y) + KD * (vr - yd) - jacdot_qdot(q, dq)); u = mass(q) @ qdd + coriolis_qdot(q, dq)
        sat += int(np.any(np.abs(u) > UMAX)); x = step(x, np.clip(u, -UMAX, UMAX), ph); k += 1
    print(f"  planar move phase {ph}: {k * T:.2f} s, {sat} of {k} steps clipped, error {1e3 * ALPHA * np.linalg.norm(fk(x[:2]) - y_goal):.2f} mm", flush=True); return x


def hold(x, dur, ph):
    for _ in range(int(round(dur / T))): x = step(x, np.zeros(2), ph)
    return x


mk = {}
x = hold(x0, D_DOWN + D_GRASP + D_UP, 0); mk["grasp"] = D_DOWN; mk["ra_start"] = len(U) * T                    # pick: arm at rest, Z axis moves
t0 = time.time(); L = cf.mpc_run(x, "cert"); print(f"RA-MPC from rest above the cube: {L['status']} in {L['t'][-1]:.2f} s, {len(L['u'])} steps, solve median {np.median(L['solve_ms']):.0f} ms ({time.time() - t0:.0f}s)", flush=True)
assert L["status"] == "reached", "RA-MPC failed"
for u in L["u"]: x = step(x, u, 1)
print(f"  replayed RA-MPC inputs: state mismatch to the MPC log {np.abs(x - L['x'][-1]).max():.1e}", flush=True)
mk["entry"] = len(U) * T
x = move(x, slot, 0.6, 2); mk["lower"] = len(U) * T                                                               # place: centre over the slot
x = move(x, slot, D_DOWN + D_GRASP + D_UP, 2, settle=0.0); mk["release"] = mk["lower"] + D_DOWN; t_end = len(U) * T   # hold while the Z axis works
# MuJoCo re-simulation of the planar part with the scaled torques
ma = mujoco.MjModel.from_xml_path(f"{HERE}/model/scara_arm.xml"); da = mujoco.MjData(ma); da.qpos[:] = x0[:2]; Q, err = [da.qpos.copy()], 0.0
for k, u in enumerate(U):
    da.ctrl[:] = C_TAU * u
    for _ in range(NSUB): mujoco.mj_step(ma, da); Q.append(da.qpos.copy())
    err = max(err, np.abs(np.r_[da.qpos, da.qvel] - X[k + 1]).max())
Q = np.array(Q); TS = np.arange(len(Q)) * ma.opt.timestep
# scripted Z axis and fingers on the 0.5 ms grid
Z = np.full(len(TS), Z0); t_ = TS
for a, b, z_a, z_b in ((0.0, mk["grasp"], Z0, Z_LOW), (mk["grasp"], mk["grasp"] + D_GRASP, Z_LOW, Z_LOW), (mk["grasp"] + D_GRASP, mk["ra_start"], Z_LOW, Z0),
                       (mk["lower"], mk["release"], Z0, Z_LOW), (mk["release"], mk["release"] + D_GRASP, Z_LOW, Z_LOW), (mk["release"] + D_GRASP, t_end, Z_LOW, Z0)):
    sel = (t_ >= a) & (t_ <= b + 1e-12); Z[sel] = z_a + (z_b - z_a) * minjerk((t_[sel] - a) / max(b - a, 1e-9))
GAP = np.where(t_ < mk["grasp"], GAP_OPEN, np.where(t_ < mk["release"], GAP_OPEN + (GAP_CLOSED - GAP_OPEN) * minjerk((t_ - mk["grasp"]) / T_FINGER), GAP_CLOSED + (GAP_OPEN - GAP_CLOSED) * minjerk((t_ - mk["release"]) / T_FINGER)))
mk["attach"] = mk["grasp"] + T_FINGER
PHQ = np.r_[0, np.repeat(PH, NSUB)]
tr = PHQ == 1; Yt = fk(Q[tr].T); psi = np.broadcast_to(cf.C.psi(*Yt), (Yt.shape[1],))
print(f"MuJoCo re-simulation (torques x {C_TAU:g}): max state difference to the paper's plant {err:.1e} (rad, rad/s)")
print(f"total {t_end:.2f} s; marks {', '.join(f'{k} {v:.2f}' for k, v in mk.items())}; transport: min psi {psi.min():.4f}, entry {mk['entry'] - mk['ra_start']:.2f} s after the start; "
      f"final slot error {1e3 * ALPHA * np.linalg.norm(fk(Q[-1]) - slot):.2f} mm; peak |tau| world {np.round(C_TAU * np.abs(np.array(U)).max(0), 3)} N m within {np.round(C_TAU * UMAX, 3)}")
np.savez(f"{HERE}/data/demo.npz", t=TS, q=Q, z=Z, gap=GAP, phase=PHQ, u=C_TAU * np.array(U), umax=C_TAU * UMAX, x_paper=np.array(X), **{f"t_{k}": v for k, v in mk.items()}, t_end=t_end)
# clearance of the moving parts to the obstacle, the other cubes and the table
ms = mujoco.MjModel.from_xml_path(f"{HERE}/model/scene_scara.xml"); ds = mujoco.MjData(ms); ft = np.zeros(6)
mov = [ms.geom(n).id for n in ("shaft", "flange", "wrist", "finger_l", "finger_r", "red_cube")] + [i for i in range(ms.ngeom) if ms.body(ms.geom_bodyid[i]).name in ("link1", "link2")]
env = {n: ms.geom(n).id for n in ("obstacle", "obj1", "obj2", "obj3", "tabletop")}; best = {n: (9.0, "", 0.0) for n in env}
home = np.r_[to_world(src), H_CUBE]
for i in range(0, len(TS), 20):
    t = TS[i]; tip = set_moving(ms, ds, Q[i], Z[i], GAP[i], home)
    cube = home if t < mk["attach"] else tip if t < mk["release"] else np.r_[tip[:2], H_CUBE]; set_moving(ms, ds, Q[i], Z[i], GAP[i], cube)
    for n, e in env.items():
        for a in mov:
            if n == "tabletop" and ms.geom(a).name == "red_cube": continue
            dist = mujoco.mj_geomDistance(ms, ds, a, e, 0.5, ft)
            if dist < best[n][0]: best[n] = (dist, ms.geom(a).name or ms.body(ms.geom_bodyid[a]).name, t)
print("clearance: " + "; ".join(f"{n} {1e3 * v[0]:.1f} mm ({v[1]}, t {v[2]:.2f})" for n, v in best.items()))
