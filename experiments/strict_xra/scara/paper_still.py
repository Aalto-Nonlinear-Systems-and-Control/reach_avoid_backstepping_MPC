"""Clean still of the red task on the SCARA for the paper (no overlay): angled view during the RA-MPC transport, with
the corridor edge, the target ring and the transport trail, cropped to the table and the arm. Uses the scene and the
demo of render_scara.py. Writes figs/fig_scara_scene.png. Usage: python paper_still.py [t_after_transport_start]"""
import os, sys, numpy as np, mujoco, imageio
import render_scara as R
dt = float(sys.argv[1]) if len(sys.argv) > 1 else 0.45; ts = R.TS + dt
R.CAMS["angled"] = dict(lookat=(0.25, 0.0, 0.30), distance=2.0, azimuth=100, elevation=-32)   # from the side of the corridor, so that the column does not hide it
ren = mujoco.Renderer(R.m, R.H, R.W); R.set_state(ts); out = []
for marks in (False, True):
    opt = mujoco.MjvOption(); opt.geomgroup[1] = int(marks); ren.update_scene(R.d, camera=R.camera("angled"), scene_option=opt)
    ren.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = int(not marks)
    if marks: R.add_trail(ren.scene, ts)
    rgb = ren.render().copy(); ren.enable_depth_rendering(); dep = ren.render().copy(); ren.disable_depth_rendering(); out.append((rgb, dep))
(a, da), (b, db) = out; mask = db < da - 3e-3; a[mask] = b[mask]
img = a[70:700, 230:1050]
imageio.imwrite(os.path.join(R.OUT, "fig_scara_scene.png"), img); print("wrote fig_scara_scene.png", img.shape, f"t = {ts:.2f} s")
