"""Video of the red task in MuJoCo: RA-MPC from rest at the pick pose (reach-avoid phase) and the placement phase,
replayed from the re-simulated states (data/red_fwd_sq/demo.npz) in model/scene_red.xml, at 0.25x speed, from an angled
camera and then from above. The red cube follows the gripper until it is released at the slot; the tool path is drawn
as a trail. Writes figs/red_fwd_sq/franka_red_task.mp4 (and still frames with --stills)."""
import os, sys, numpy as np, mujoco, imageio
from PIL import Image, ImageDraw, ImageFont
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figs", "red_fwd_sq"); D = np.load(os.path.join(HERE, "data", "red_fwd_sq", "demo.npz"))
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_red.xml")); d = mujoco.MjData(m); W, H = 1280, 720
T, Q, PH, TE = D["t"], D["q"], np.r_[D["phase"][0], D["phase"]], float(D["t_entry"]); TEND = T[-1]; SLOW, FPS, HOLD0, HOLD1 = 0.25, 60, 0.8, 1.4
CAMS = {"angled": dict(lookat=(0.40, 0.0, 0.26), distance=1.50, azimuth=12, elevation=-50), "top": "top"}      # "top": fixed camera of the scene (near-orthographic)
try: FONT = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 24); FONT_S = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 19)
except Exception: FONT = FONT_S = ImageFont.load_default()
tcp_id, hand_id, mocap = m.site("tcp").id, m.body("hand").id, m.body("red_cube").mocapid[0]


def set_state(ts):
    i = min(int(round(ts / m.opt.timestep)), len(Q) - 1); d.qpos[:] = Q[i]; mujoco.mj_forward(m, d)
    d.mocap_pos[mocap] = d.site_xpos[tcp_id]; d.mocap_quat[mocap] = d.xquat[hand_id]; mujoco.mj_forward(m, d); return i


TRAIL = []
for ts in np.arange(0, TEND + 1e-9, 0.01): set_state(ts); TRAIL.append((ts, d.site_xpos[tcp_id].copy()))


def add_trail(scn, ts):
    for t_, p in TRAIL:
        if t_ > ts or scn.ngeom >= scn.maxgeom: break
        g = scn.geoms[scn.ngeom]; rgba = np.array([0.95, 0.55, 0.05, 0.9]) if t_ <= TE else np.array([0.2, 0.6, 0.3, 0.9])
        mujoco.mjv_initGeom(g, mujoco.mjtGeom.mjGEOM_SPHERE, np.array([0.0045, 0, 0]), p, np.eye(3).ravel(), rgba); scn.ngeom += 1


def overlay(img, ts, cam):
    im = Image.fromarray(img).convert("RGBA"); band = Image.new("RGBA", im.size, (0, 0, 0, 0)); ImageDraw.Draw(band).rectangle([0, 0, W, 110], fill=(255, 255, 255, 200))
    im = Image.alpha_composite(im, band).convert("RGB"); dr = ImageDraw.Draw(im); phase = "reach-avoid (RA-MPC)" if ts < TE else "placement (target reached)"
    dr.text((24, 18), "RA-MPC on a Panda-like arm in MuJoCo, red task", font=FONT, fill=(20, 20, 20))
    dr.text((24, 52), f"t = {min(ts, TEND):.2f} s   ({SLOW:g}x speed)   phase: {phase}", font=FONT_S, fill=(40, 40, 40))
    dr.text((24, 78), f"task torque limits {', '.join(f'{u:g}' for u in D['umax'])} N m;  blue: safe corridor edge X_S at 0.30 m;  red disc: target X_T", font=FONT_S, fill=(70, 70, 70))
    dr.text((W - 150, H - 40), f"view: {cam}", font=FONT_S, fill=(90, 90, 90)); return np.asarray(im)


def frame(ren, cam, ts):
    set_state(min(ts, TEND)); c = mujoco.MjvCamera()
    if isinstance(CAMS[cam], str): c.type = mujoco.mjtCamera.mjCAMERA_FIXED; c.fixedcamid = m.camera(CAMS[cam]).id
    else:
        c.type = mujoco.mjtCamera.mjCAMERA_FREE
        for k, v in CAMS[cam].items(): setattr(c, k, v) if k != "lookat" else c.lookat.__setitem__(slice(None), v)
    opt = mujoco.MjvOption(); opt.geomgroup[:] = 1; ren.update_scene(d, camera=c, scene_option=opt); add_trail(ren.scene, ts); return overlay(ren.render(), ts, cam)


if __name__ == "__main__":
    ren = mujoco.Renderer(m, H, W)
    if "--stills" in sys.argv:
        for cam in CAMS:
            for ts in (0.0, 0.5, TEND): imageio.imwrite(os.path.join(OUT, f"still_{cam}_{ts:.2f}.png"), frame(ren, cam, ts))
        print("stills written"); sys.exit()
    wr = imageio.get_writer(os.path.join(OUT, "franka_red_task.mp4"), fps=FPS, codec="libx264", quality=8, macro_block_size=8); nf = 0
    for cam in ("angled", "top"):
        for tau in np.arange(0, HOLD0 + TEND / SLOW + HOLD1, 1 / FPS): wr.append_data(frame(ren, cam, max(0.0, tau - HOLD0) * SLOW)); nf += 1
    wr.close(); print(f"wrote franka_red_task.mp4: {nf} frames, {nf/FPS:.1f} s at {FPS} fps, {W}x{H}")
