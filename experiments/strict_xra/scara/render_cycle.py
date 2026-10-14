"""Video of the full kitting cycle on the SCARA (eight tasks, data/cycle.npz) in model/scene_scara_cycle.xml: the angled
front view with the top view (arm as a translucent ghost, orientation of the paper figures) as an inset at the top right. Corridor edge,
target ring and transport trail of the current task are drawn without shadows (second pass, composited by depth).
Pick, place and transit play at 1x, the certified transport at 0.5x. Writes figs/scara_cycle.mp4 (stills with --stills)."""
import os, sys, numpy as np, mujoco, imageio
from PIL import Image, ImageDraw, ImageFont
import cycle_common as cc
from scara_common import HERE, Z0
OUT = os.path.join(HERE, "figs"); m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_scara_cycle.xml")); d = mujoco.MjData(m)
W, H, WI, HI = 1280, 720, 448, 252; FPS, HOLD0, HOLD1 = 60, 0.8, 1.5; SPEED = (1.0, 0.5)                    # (pick, place, transit), transport
ANG = dict(lookat=(0.25, -0.10, 0.40), distance=2.15, azimuth=215, elevation=-30); GHOST = 0.30
ARM_MAT = [m.material(n).id for n in ("scara_white", "scara_dark", "scara_metal")]; ALPHA0 = m.mat_rgba[ARM_MAT, 3].copy()
try: FONT = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 21); FONT_S = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 16)
except Exception: FONT = FONT_S = ImageFont.load_default()


def cam_free():
    c = mujoco.MjvCamera(); c.type = mujoco.mjtCamera.mjCAMERA_FREE
    for k, v in ANG.items(): setattr(c, k, v) if k != "lookat" else c.lookat.__setitem__(slice(None), v)
    return c


def cam_top():
    c = mujoco.MjvCamera(); c.type = mujoco.mjtCamera.mjCAMERA_FIXED; c.fixedcamid = m.camera("top").id; return c


def frame(ren, ren_i, ts):
    ts = min(ts, cc.TEND); i, ph = cc.set_state(m, d, ts); out = []
    for marks in (False, True):
        ren.update_scene(d, camera=cam_free(), scene_option=mujoco.MjvOption()); ren.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = int(not marks)
        if marks: cc.add_marks(ren.scene, ts, i)
        rgb = ren.render().copy(); ren.enable_depth_rendering(); dep = ren.render().copy(); ren.disable_depth_rendering(); out.append((rgb, dep))
    (a, da), (b, db) = out; mask = db < da - 3e-3; a[mask] = b[mask]
    m.mat_rgba[ARM_MAT, 3] = GHOST; ren_i.update_scene(d, camera=cam_top(), scene_option=mujoco.MjvOption()); ren_i.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = 0
    cc.add_marks(ren_i.scene, ts, i); inset = ren_i.render().copy(); m.mat_rgba[ARM_MAT, 3] = ALPHA0
    im = Image.fromarray(a).convert("RGBA"); band = Image.new("RGBA", im.size, (0, 0, 0, 0)); dr = ImageDraw.Draw(band); dr.rectangle([0, 0, W, 84], fill=(255, 255, 255, 190))
    x0, y0 = W - WI - 18, 84 + 16; dr.rectangle([x0 - 3, y0 - 3, x0 + WI + 2, y0 + HI + 2], fill=(255, 255, 255, 230)); im = Image.alpha_composite(im, band).convert("RGB")
    im.paste(Image.fromarray(inset), (x0, y0)); dr = ImageDraw.Draw(im); sp_ = SPEED[ph == 1]
    dr.text((20, 9), "RA-MPC on a SCARA arm in MuJoCo: kitting cycle, eight certified transports", font=FONT, fill=(20, 20, 20))
    dr.text((20, 37), f"{cc.label(i)}   |   phase: {cc.PHASE[ph]}   |   t = {ts:.2f} s ({sp_:g}x speed)", font=FONT_S, fill=(40, 40, 40))
    dr.text((20, 59), f"torque limits {', '.join(f'{u:.3g}' for u in cc.UMAX_W)} N m for all tasks;  blue: edge of the current task's safe corridor X_S at {Z0:.2f} m;  red ring: its target X_T", font=FONT_S, fill=(70, 70, 70))
    dr.rectangle([x0, y0, x0 + 78, y0 + 24], fill=(255, 255, 255)); dr.text((x0 + 8, y0 + 4), "top view", font=FONT_S, fill=(60, 60, 60)); return np.asarray(im)


def times():
    out = [0.0] * int(HOLD0 * FPS); ts = 0.0
    while ts < cc.TEND: out.append(ts); ts += SPEED[int(cc.PH[cc.idx(ts)]) == 1] / FPS
    return out + [cc.TEND] * int(HOLD1 * FPS)


if __name__ == "__main__":
    ren, ren_i = mujoco.Renderer(m, H, W), mujoco.Renderer(m, HI, WI)
    if "--stills" in sys.argv:
        for i, inf in enumerate(cc.INFO):
            if i in (0, 1, 4, 7):
                for ts in (inf["t_grasp"], 0.5 * (inf["t_ra"] + inf["t_entry"]), inf["t_release"]): imageio.imwrite(os.path.join(OUT, f"cycle_still_{i}_{ts:.2f}.png"), frame(ren, ren_i, ts))
        print("stills written"); sys.exit()
    wr = imageio.get_writer(os.path.join(OUT, "scara_cycle.mp4"), fps=FPS, codec="libx264", quality=8, macro_block_size=8); nf = 0
    for ts in times(): wr.append_data(frame(ren, ren_i, ts)); nf += 1
    wr.close(); print(f"wrote scara_cycle.mp4: {nf} frames, {nf / FPS:.1f} s at {FPS} fps, {W}x{H}")
