"""Video of the red task in MuJoCo as a full pick and place: pick (descend, grasp, lift; not certified), transport by
RA-MPC from rest into the target X_T (certified), place (centre over the slot, lower, release, retract; not certified).
Replayed from the 0.5 ms plant states (data/red_fwd_sq/demo.npz) in model/scene_red.xml, first from an angled front
camera, then from the fixed top camera (the orientation of the paper figures). The red cube rests on the table until it
is grasped, follows the gripper and rests at the slot after release; the transport path is drawn as a trail. Pick and
place play faster than the transport (SPEED). Writes figs/red_fwd_sq/franka_red_task.mp4 (and still frames with --stills)."""
import os, sys, numpy as np, mujoco, imageio
from PIL import Image, ImageDraw, ImageFont
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figs", "red_fwd_sq"); D = np.load(os.path.join(HERE, "data", "red_fwd_sq", "demo.npz"))
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_red.xml")); m.light_active[m.light("top").id] = 0; d = mujoco.MjData(m); W, H = 1280, 720   # the light of Menagerie's panda.xml off
Q, PH = D["q"], np.r_[D["phase"][0], D["phase"]]; HC = float(D["h_cube"]); TG, TS, TE, TR = (float(D[k]) for k in ("t_grasp", "t_ra_start", "t_entry", "t_release")); TEND = (len(Q) - 1) * m.opt.timestep
FPS, HOLD0, HOLD1 = 60, 0.8, 1.4; SPEED = {"angled": (0.5, 0.25), "top": (1.0, 0.25)}       # (pick and place, transport) playback speeds
CAMS = {"angled": dict(lookat=(0.25, -0.10, 0.31), distance=1.95, azimuth=215, elevation=-30), "top": "top"}
PHASE = {0: "pick (not certified)", 1: "transport by RA-MPC (certified reach-avoid)", 2: "place (not certified)"}
try: FONT = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 21); FONT_S = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 16)
except Exception: FONT = FONT_S = ImageFont.load_default()
tcp_id, hand_id, mocap = m.site("tcp").id, m.body("hand").id, m.body("red_cube").mocapid[0]


def pose(ts):
    d.qpos[:] = Q[min(int(round(ts / m.opt.timestep)), len(Q) - 1)]; mujoco.mj_forward(m, d); return d.site_xpos[tcp_id].copy(), d.xquat[hand_id].copy()
P_G, R_G = pose(TG); P_R, R_R = pose(TR); P_G[2] = P_R[2] = HC                       # the cube rests on the table before the grasp and after the release


def set_state(ts):
    i = min(int(round(ts / m.opt.timestep)), len(Q) - 1); p, r = pose(ts)
    if ts < TG: p, r = P_G, R_G
    elif ts >= TR: p, r = P_R, R_R
    d.mocap_pos[mocap] = p; d.mocap_quat[mocap] = r; mujoco.mj_forward(m, d); return i


TRAIL = []
for ts in np.arange(TS, TE + 1e-9, 0.01): TRAIL.append((ts, pose(ts)[0]))          # the certified transport only


def add_trail(scn, ts):
    for t_, p in TRAIL:
        if t_ > ts or scn.ngeom >= scn.maxgeom: break
        mujoco.mjv_initGeom(scn.geoms[scn.ngeom], mujoco.mjtGeom.mjGEOM_SPHERE, np.array([0.0045, 0, 0]), p, np.eye(3).ravel(), np.array([0.95, 0.55, 0.05, 1.0])); scn.ngeom += 1


def overlay(img, ts, cam):
    im = Image.fromarray(img).convert("RGBA"); band = Image.new("RGBA", im.size, (0, 0, 0, 0)); ImageDraw.Draw(band).rectangle([0, 0, W, 84], fill=(255, 255, 255, 190))
    im = Image.alpha_composite(im, band).convert("RGB"); dr = ImageDraw.Draw(im); ph = int(PH[min(int(round(ts / m.opt.timestep)), len(PH) - 1)]); sp = SPEED[cam][ph == 1]
    dr.text((20, 9), "RA-MPC on a Panda-like arm in MuJoCo, red task", font=FONT, fill=(20, 20, 20))
    dr.text((20, 37), f"t = {min(ts, TEND):.2f} s   ({sp:g}x speed)   phase: {PHASE[ph]}", font=FONT_S, fill=(40, 40, 40))
    dr.text((20, 59), f"torque limits {', '.join(f'{u:g}' for u in D['umax'])} N m;  blue: edge of the safe corridor X_S at 0.30 m;  red ring: target X_T", font=FONT_S, fill=(70, 70, 70))
    dr.text((W - 130, H - 34), f"view: {cam}", font=FONT_S, fill=(230, 230, 230)); return np.asarray(im)


def camera(cam):
    c = mujoco.MjvCamera()
    if isinstance(CAMS[cam], str): c.type = mujoco.mjtCamera.mjCAMERA_FIXED; c.fixedcamid = m.camera(CAMS[cam]).id
    else:
        c.type = mujoco.mjtCamera.mjCAMERA_FREE
        for k, v in CAMS[cam].items(): setattr(c, k, v) if k != "lookat" else c.lookat.__setitem__(slice(None), v)
    return c


def frame(ren, cam, ts):
    """Two passes: the physical scene with shadows (marks hidden), then the marks and the trail without shadows, composited by depth.
    Default geom groups 0-2: the arm's visual meshes are group 2; its grey collision meshes (group 3) stay hidden."""
    set_state(min(ts, TEND)); out = []
    for marks in (False, True):
        opt = mujoco.MjvOption(); opt.geomgroup[1] = int(marks); ren.update_scene(d, camera=camera(cam), scene_option=opt)
        ren.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = int(not marks)
        if marks: add_trail(ren.scene, ts)
        rgb = ren.render().copy(); ren.enable_depth_rendering(); dep = ren.render().copy(); ren.disable_depth_rendering(); out.append((rgb, dep))
    (a, da), (b, db) = out; mask = db < da - 3e-3; a[mask] = b[mask]; return overlay(a, ts, cam)   # depth noise of the two passes is below 1 mm


def times(cam):
    """Simulation time of each video frame: hold, play at the phase's speed, hold."""
    out = [0.0] * int(HOLD0 * FPS); ts = 0.0
    while ts < TEND:
        out.append(ts); ph = int(PH[min(int(round(ts / m.opt.timestep)), len(PH) - 1)]); ts += SPEED[cam][ph == 1] / FPS
    return out + [TEND] * int(HOLD1 * FPS)


if __name__ == "__main__":
    ren = mujoco.Renderer(m, H, W)
    if "--stills" in sys.argv:
        for cam in CAMS:
            for ts in (0.0, TG, 0.5 * (TS + TE), TE + 0.3, TEND): imageio.imwrite(os.path.join(OUT, f"still_{cam}_{ts:.2f}.png"), frame(ren, cam, ts))
        print("stills written"); sys.exit()
    wr = imageio.get_writer(os.path.join(OUT, "franka_red_task.mp4"), fps=FPS, codec="libx264", quality=8, macro_block_size=8); nf = 0
    for cam in ("angled", "top"):
        for ts in times(cam): wr.append_data(frame(ren, cam, ts)); nf += 1
    wr.close(); print(f"wrote franka_red_task.mp4: {nf} frames, {nf/FPS:.1f} s at {FPS} fps, {W}x{H}")
