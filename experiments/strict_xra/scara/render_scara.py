"""Video of the red task on the SCARA in MuJoCo: pick (Z axis down, grasp, up; the arm at rest), transport by RA-MPC
from rest into the target X_T (certified), place (centre over the slot, Z axis down, release, up; not certified).
Replayed from data/demo.npz in model/scene_scara.xml, first from an angled front camera, then from the fixed top camera
(the orientation of the paper figures). In the angled view the marks floating at the transport height and the transport
trail are drawn in a second pass without shadows and composited by depth; the top view shows the arm as a translucent
ghost without shadows, so that the corridor, the trail and the cubes below it stay visible (as in the paper figures). Pick and place play faster than the transport (SPEED).
Writes figs/scara_red_task.mp4 (and still frames with --stills)."""
import os, sys, numpy as np, mujoco, imageio
from PIL import Image, ImageDraw, ImageFont
from scara_common import HERE, H_CUBE, Z0, set_moving
OUT = os.path.join(HERE, "figs"); D = np.load(os.path.join(HERE, "data", "demo.npz"))
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_scara.xml")); d = mujoco.MjData(m); W, H = 1280, 720
Q, Z, GAP, PH = D["q"], D["z"], D["gap"], D["phase"]; TA, TS, TE, TR = (float(D[k]) for k in ("t_attach", "t_ra_start", "t_entry", "t_release")); TEND = (len(Q) - 1) * m.opt.timestep
FPS, HOLD0, HOLD1 = 60, 0.8, 1.4; SPEED = {"angled": (0.5, 0.25), "top": (1.0, 0.25)}       # (pick and place, transport) playback speeds
CAMS = {"angled": dict(lookat=(0.25, -0.10, 0.40), distance=2.15, azimuth=215, elevation=-30), "top": "top"}; GHOST = 0.30
PHASE = {0: "pick (Z axis; arm at rest)", 1: "transport by RA-MPC (certified reach-avoid)", 2: "place (not certified)"}
try: FONT = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 21); FONT_S = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 16)
except Exception: FONT = FONT_S = ImageFont.load_default()
idx = lambda ts: min(int(round(ts / m.opt.timestep)), len(Q) - 1)
ARM_MAT = [m.material(n).id for n in ("scara_white", "scara_dark", "scara_metal")]; ALPHA0 = m.mat_rgba[ARM_MAT, 3].copy()
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


def overlay(img, ts, cam):
    im = Image.fromarray(img).convert("RGBA"); band = Image.new("RGBA", im.size, (0, 0, 0, 0)); ImageDraw.Draw(band).rectangle([0, 0, W, 84], fill=(255, 255, 255, 190))
    im = Image.alpha_composite(im, band).convert("RGB"); dr = ImageDraw.Draw(im); ph = int(PH[idx(ts)]); sp = SPEED[cam][ph == 1]
    dr.text((20, 9), "RA-MPC on a SCARA arm in MuJoCo, red task", font=FONT, fill=(20, 20, 20))
    dr.text((20, 37), f"t = {min(ts, TEND):.2f} s   ({sp:g}x speed)   phase: {PHASE[ph]}", font=FONT_S, fill=(40, 40, 40))
    dr.text((20, 59), f"torque limits {', '.join(f'{u:.3g}' for u in D['umax'])} N m;  blue: edge of the safe corridor X_S at {Z0:.2f} m;  red ring: target X_T", font=FONT_S, fill=(70, 70, 70))
    dr.text((W - 130, H - 34), f"view: {cam}", font=FONT_S, fill=(230, 230, 230)); return np.asarray(im)


def camera(cam):
    c = mujoco.MjvCamera()
    if isinstance(CAMS[cam], str): c.type = mujoco.mjtCamera.mjCAMERA_FIXED; c.fixedcamid = m.camera(CAMS[cam]).id
    else:
        c.type = mujoco.mjtCamera.mjCAMERA_FREE
        for k, v in CAMS[cam].items(): setattr(c, k, v) if k != "lookat" else c.lookat.__setitem__(slice(None), v)
    return c


def frame(ren, cam, ts):
    """Angled view: two passes, the physical scene with shadows (marks hidden), then the marks and the trail without shadows,
    composited by depth. Top view: one pass without shadows, the arm translucent."""
    set_state(min(ts, TEND)); out = []
    if cam == "top":
        m.mat_rgba[ARM_MAT, 3] = GHOST; opt = mujoco.MjvOption(); ren.update_scene(d, camera=camera(cam), scene_option=opt)
        ren.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = 0; add_trail(ren.scene, ts); img = ren.render().copy(); m.mat_rgba[ARM_MAT, 3] = ALPHA0; return overlay(img, ts, cam)
    for marks in (False, True):
        opt = mujoco.MjvOption(); opt.geomgroup[1] = int(marks); ren.update_scene(d, camera=camera(cam), scene_option=opt)
        ren.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = int(not marks)
        if marks: add_trail(ren.scene, ts)
        rgb = ren.render().copy(); ren.enable_depth_rendering(); dep = ren.render().copy(); ren.disable_depth_rendering(); out.append((rgb, dep))
    (a, da), (b, db) = out; mask = db < da - 3e-3; a[mask] = b[mask]; return overlay(a, ts, cam)   # depth noise of the two passes is below 1 mm


def times(cam):
    out = [0.0] * int(HOLD0 * FPS); ts = 0.0
    while ts < TEND: out.append(ts); ts += SPEED[cam][int(PH[idx(ts)]) == 1] / FPS
    return out + [TEND] * int(HOLD1 * FPS)


if __name__ == "__main__":
    ren = mujoco.Renderer(m, H, W)
    if "--stills" in sys.argv:
        for cam in CAMS:
            for ts in (0.0, TA, 0.5 * (TS + TE), TR, TEND): imageio.imwrite(os.path.join(OUT, f"still_{cam}_{ts:.2f}.png"), frame(ren, cam, ts))
        print("stills written"); sys.exit()
    wr = imageio.get_writer(os.path.join(OUT, "scara_red_task.mp4"), fps=FPS, codec="libx264", quality=8, macro_block_size=8); nf = 0
    for cam in ("angled", "top"):
        for ts in times(cam): wr.append_data(frame(ren, cam, ts)); nf += 1
    wr.close(); print(f"wrote scara_red_task.mp4: {nf} frames, {nf / FPS:.1f} s at {FPS} fps, {W}x{H}")
