"""Fig. 6 of the paper: the red task on the SCARA in MuJoCo during the RA-MPC transport, replayed from data/demo.npz in
model/scene_scara.xml, seen from an angled camera and cropped to the table and the arm. The corridor edge, the target ring and
the transport trail are drawn in a second pass without shadows and composited by depth. Writes ../figs/fig_scara_scene.png.
Usage: python render_still.py [t_after_transport_start]"""
import os, sys, numpy as np, mujoco, imageio
from scara_common import HERE, H_CUBE, set_moving
OUT = os.path.join(HERE, "..", "figs"); D = np.load(os.path.join(HERE, "data", "demo.npz"))
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_scara.xml")); d = mujoco.MjData(m); W, H = 1280, 720
Q, Z, GAP = D["q"], D["z"], D["gap"]; TA, TS, TE, TR = (float(D[k]) for k in ("t_attach", "t_ra_start", "t_entry", "t_release"))
idx = lambda ts: min(int(round(ts / m.opt.timestep)), len(Q) - 1)
HOME = set_moving(m, d, Q[0], H_CUBE, GAP[0], np.zeros(3)); REL = set_moving(m, d, Q[idx(TR)], H_CUBE, GAP[idx(TR)], np.zeros(3))


def set_state(ts):
    i = idx(ts); tip = set_moving(m, d, Q[i], Z[i], GAP[i], HOME)
    cube = HOME if ts < TA else tip if ts < TR else REL; set_moving(m, d, Q[i], Z[i], GAP[i], cube); return i


TRAIL = []
for ts in np.arange(TS, TE + 1e-9, 0.01): set_state(ts); TRAIL.append((ts, d.mocap_pos[m.body("quill").mocapid[0]].copy()))


def add_trail(scn, ts):
    for t_, p in TRAIL:
        if t_ > ts or scn.ngeom >= scn.maxgeom: break
        mujoco.mjv_initGeom(scn.geoms[scn.ngeom], mujoco.mjtGeom.mjGEOM_SPHERE, np.array([0.0045, 0, 0]), p, np.eye(3).ravel(), np.array([0.95, 0.55, 0.05, 1.0])); scn.ngeom += 1


cam = mujoco.MjvCamera(); cam.type = mujoco.mjtCamera.mjCAMERA_FREE; cam.lookat[:] = (0.25, 0.0, 0.30); cam.distance, cam.azimuth, cam.elevation = 2.0, 100, -32   # from the side of the corridor, so that the column does not hide it
dt = float(sys.argv[1]) if len(sys.argv) > 1 else 0.45; ts = TS + dt
ren = mujoco.Renderer(m, H, W); set_state(ts); out = []
for marks in (False, True):
    opt = mujoco.MjvOption(); opt.geomgroup[1] = int(marks); ren.update_scene(d, camera=cam, scene_option=opt)
    ren.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = int(not marks)
    if marks: add_trail(ren.scene, ts)
    rgb = ren.render().copy(); ren.enable_depth_rendering(); dep = ren.render().copy(); ren.disable_depth_rendering(); out.append((rgb, dep))
(a, da), (b, db) = out; mask = db < da - 3e-3; a[mask] = b[mask]
img = a[70:700, 230:1050]
imageio.imwrite(os.path.join(OUT, "fig_scara_scene.png"), img); print("wrote fig_scara_scene.png", img.shape, f"t = {ts:.2f} s")
