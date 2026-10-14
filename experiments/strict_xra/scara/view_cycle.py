"""Local MuJoCo viewer of the full kitting cycle on the SCARA (data/cycle.npz), in a loop. macOS: ./view_cycle.sh [speed]
(default 1.0 = real time). Corridor edge, target ring and trail of the current task as user geoms; shadows off."""
import os, sys, time, numpy as np, mujoco, mujoco.viewer
import cycle_common as cc
from scara_common import HERE
SPEED = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_scara_cycle.xml")); d = mujoco.MjData(m); m.light_castshadow[:] = 0
with mujoco.viewer.launch_passive(m, d) as v:
    v.cam.lookat[:] = (0.25, -0.10, 0.40); v.cam.distance = 2.15; v.cam.azimuth = 215; v.cam.elevation = -30
    while v.is_running():
        t0 = time.time()
        while v.is_running():
            ts = (time.time() - t0) * SPEED - 0.3
            if ts > cc.TEND + 0.8: break
            i, ph = cc.set_state(m, d, min(max(ts, 0.0), cc.TEND))
            with v.lock(): v.user_scn.ngeom = 0; cc.add_marks(v.user_scn, ts, i)
            v.sync(); time.sleep(1 / 120)
