"""Local MuJoCo viewer: replays the red-task demo (data/red_fwd_sq/demo.npz) in model/scene_red.xml in a loop.
macOS:  mjpython view_red.py [speed]      (default speed 0.5; 1.0 = real time)
Pick, RA-MPC transport and place as in render_red.py; the red cube rests on the table, follows the gripper and rests at
the slot; the transport path is drawn as a trail. Shadows are off here, as the marks floating at the transport height
would cast them on the table. Mouse: rotate, zoom and pan as in the MuJoCo viewer; close the window to stop."""
import os, sys, time, numpy as np, mujoco, mujoco.viewer
HERE = os.path.dirname(os.path.abspath(__file__)); SPEED = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_red.xml")); m.light_active[m.light("top").id] = 0; d = mujoco.MjData(m); D = np.load(os.path.join(HERE, "data", "red_fwd_sq", "demo.npz"))   # the light of Menagerie's panda.xml off
Q, HC = D["q"], float(D["h_cube"]); TG, TS, TE, TR = (float(D[k]) for k in ("t_grasp", "t_ra_start", "t_entry", "t_release")); TEND = (len(Q) - 1) * m.opt.timestep
tcp, hand, mocap = m.site("tcp").id, m.body("hand").id, m.body("red_cube").mocapid[0]; m.light_castshadow[:] = 0


def pose(ts):
    d.qpos[:] = Q[min(int(round(ts / m.opt.timestep)), len(Q) - 1)]; mujoco.mj_forward(m, d); return d.site_xpos[tcp].copy(), d.xquat[hand].copy()
P_G, R_G = pose(TG); P_R, R_R = pose(TR); P_G[2] = P_R[2] = HC


def set_state(ts):
    p, r = pose(ts)
    if ts < TG: p, r = P_G, R_G
    elif ts >= TR: p, r = P_R, R_R
    d.mocap_pos[mocap] = p; d.mocap_quat[mocap] = r; mujoco.mj_forward(m, d)


TRAIL = [(ts, pose(ts)[0]) for ts in np.arange(TS, TE + 1e-9, 0.01)]

with mujoco.viewer.launch_passive(m, d) as v:
    v.cam.lookat[:] = (0.25, -0.10, 0.31); v.cam.distance = 1.95; v.cam.azimuth = 215; v.cam.elevation = -30
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
                    mujoco.mjv_initGeom(v.user_scn.geoms[v.user_scn.ngeom], mujoco.mjtGeom.mjGEOM_SPHERE, np.array([0.0045, 0, 0]), p, np.eye(3).ravel(), np.array([0.95, 0.55, 0.05, 1.0])); v.user_scn.ngeom += 1
            v.sync(); time.sleep(1 / 120)
