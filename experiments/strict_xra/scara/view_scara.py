"""Local MuJoCo viewer: replays the SCARA red-task demo (data/demo.npz) in model/scene_scara.xml in a loop.
macOS:  ./view_scara.sh [speed]      (default speed 0.5; 1.0 = real time)
Pick, RA-MPC transport and place as in render_scara.py; the transport path is drawn as a trail. Shadows are off here,
as the marks floating at the transport height would cast them on the table. Mouse: rotate, zoom and pan as in the
MuJoCo viewer; close the window to stop."""
import os, sys, time, numpy as np, mujoco, mujoco.viewer
from scara_common import HERE, H_CUBE, set_moving
SPEED = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_scara.xml")); d = mujoco.MjData(m); D = np.load(os.path.join(HERE, "data", "demo.npz")); m.light_castshadow[:] = 0
Q, Z, GAP = D["q"], D["z"], D["gap"]; TA, TS, TE, TR = (float(D[k]) for k in ("t_attach", "t_ra_start", "t_entry", "t_release")); TEND = (len(Q) - 1) * m.opt.timestep
idx = lambda ts: min(int(round(ts / m.opt.timestep)), len(Q) - 1)
HOME = set_moving(m, d, Q[0], H_CUBE, GAP[0], np.zeros(3)); REL = set_moving(m, d, Q[idx(TR)], H_CUBE, GAP[idx(TR)], np.zeros(3))


def set_state(ts):
    i = idx(ts); tip = set_moving(m, d, Q[i], Z[i], GAP[i], HOME)
    set_moving(m, d, Q[i], Z[i], GAP[i], HOME if ts < TA else tip if ts < TR else REL)


TRAIL = []
for ts in np.arange(TS, TE + 1e-9, 0.01): set_state(ts); TRAIL.append((ts, d.mocap_pos[m.body("quill").mocapid[0]].copy()))

with mujoco.viewer.launch_passive(m, d) as v:
    v.cam.lookat[:] = (0.25, -0.10, 0.40); v.cam.distance = 2.15; v.cam.azimuth = 215; v.cam.elevation = -30
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
