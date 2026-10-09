"""Local MuJoCo viewer: replays the red-task demo (data/red_fwd_sq/demo.npz) in model/scene_red.xml in a loop.
macOS:  mjpython view_red.py [speed]      (default speed 0.25; 1.0 = real time)
The red cube follows the gripper; the tool path is drawn as a trail (orange: RA-MPC, green: placement). Mouse: rotate,
zoom and pan as in the MuJoCo viewer; close the window to stop."""
import os, sys, time, numpy as np, mujoco, mujoco.viewer
HERE = os.path.dirname(os.path.abspath(__file__)); SPEED = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_red.xml")); d = mujoco.MjData(m); D = np.load(os.path.join(HERE, "data", "red_fwd_sq", "demo.npz"))
Q, TE = D["q"], float(D["t_entry"]); TEND = D["t"][-1]; tcp, hand, mocap = m.site("tcp").id, m.body("hand").id, m.body("red_cube").mocapid[0]


def set_state(ts):
    d.qpos[:] = Q[min(int(round(ts / m.opt.timestep)), len(Q) - 1)]; mujoco.mj_forward(m, d)
    d.mocap_pos[mocap] = d.site_xpos[tcp]; d.mocap_quat[mocap] = d.xquat[hand]; mujoco.mj_forward(m, d)


TRAIL = []
for ts in np.arange(0, TEND + 1e-9, 0.01): set_state(ts); TRAIL.append((ts, d.site_xpos[tcp].copy()))

with mujoco.viewer.launch_passive(m, d) as v:
    v.cam.lookat[:] = (0.40, 0.0, 0.26); v.cam.distance = 1.5; v.cam.azimuth = 12; v.cam.elevation = -50
    while v.is_running():
        t0 = time.time()
        while v.is_running():
            ts = (time.time() - t0) * SPEED - 0.3
            if ts > TEND + 0.6: break
            set_state(min(max(ts, 0.0), TEND))
            with v.lock():
                v.user_scn.ngeom = 0
                for t_, p in TRAIL:
                    if t_ > ts or v.user_scn.ngeom >= v.user_scn.maxgeom: break
                    mujoco.mjv_initGeom(v.user_scn.geoms[v.user_scn.ngeom], mujoco.mjtGeom.mjGEOM_SPHERE, np.array([0.0045, 0, 0]), p, np.eye(3).ravel(),
                                        np.array([0.95, 0.55, 0.05, 0.9]) if t_ <= TE else np.array([0.2, 0.6, 0.3, 0.9])); v.user_scn.ngeom += 1
            v.sync(); time.sleep(1 / 120)
