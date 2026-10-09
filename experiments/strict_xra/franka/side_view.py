"""Side view of the Panda-like arm in the vertical plane of the arm (radial distance from the base axis against height),
for poses along the red task's RA-MPC trajectory from the farthest common start; joints and link offsets from MuJoCo's
forward kinematics of model/panda_simple.xml. Also a top view of the same poses for reference."""
import os, pickle, numpy as np, mujoco, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config_franka as cf
from franka_model import MJ, NQ

M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); logs = M["logs"]; X0 = M["X0"]; Y0 = cf.C.y(X0)
rr = np.array([L["status"] == "reached" for L in logs]); d2 = ((Y0[:2] - cf.TARGET[:, None])**2).sum(0); i0 = int(np.where(rr)[0][np.argmax(d2[rr])])
L = logs[i0]; Xs = np.asarray(L["x"]); idx = np.unique(np.linspace(0, len(Xs) - 1, 5).round().astype(int))
CHAIN = ["link0", "link1", "link2", "link3", "link4", "link5", "link6", "link7", "hand"]; NAMES = {"link1": "shoulder", "link4": "elbow", "link6": "wrist"}
d = mujoco.MjData(MJ)


def frames(x):
    d.qpos[:] = x[:NQ]; mujoco.mj_forward(MJ, d)
    P = np.array([d.body(b).xpos for b in CHAIN] + [d.site("tcp").xpos]); return P


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 6.4), gridspec_kw=dict(width_ratios=[1.25, 1]), constrained_layout=True)
cols = plt.cm.Purples(np.linspace(0.35, 0.95, len(idx)))
for c, k in zip(cols, idx):
    P = frames(Xs[k]); az = np.arctan2(P[-1, 1], P[-1, 0]); r = P[:, 0] * np.cos(az) + P[:, 1] * np.sin(az); z = P[:, 2]
    ax.plot(r, z, "-", color=c, lw=4, solid_capstyle="round", alpha=0.9, zorder=3, label=f"t = {k*cf.T:.2f} s")
    ax.plot(r[[1, 4, 6]], z[[1, 4, 6]], "o", color=c, mec="k", mew=0.6, ms=9, zorder=4); ax.plot(r[-1], z[-1], "s", color=c, mec="k", mew=0.6, ms=7, zorder=4)
    ax2.plot([0, -P[-1, 1]], [0, P[-1, 0]], "-", color=c, lw=4, alpha=0.6, solid_capstyle="round"); ax2.plot(-P[6, 1], P[6, 0], "o", color=c, mec="k", mew=0.5, ms=7)
P = frames(Xs[idx[0]]); az = np.arctan2(P[-1, 1], P[-1, 0]); r = P[:, 0] * np.cos(az) + P[:, 1] * np.sin(az)
for j, nm in ((1, "shoulder (joint 2)"), (4, "elbow (joint 4)"), (6, "wrist (joint 6)"), (9, "tool centre point")):
    ax.annotate(nm, (r[j], P[j, 2]), xytext=(r[j] - 0.21 if j != 9 else r[j] + 0.03, P[j, 2] + 0.06), fontsize=10, arrowprops=dict(arrowstyle="-", color="0.4", lw=0.8))
ax.axhspan(cf.Z0 - 0.03, cf.Z0 + 0.03, color="#e6f2ff", zorder=1); ax.axhline(cf.Z0, color="#2166ac", lw=1.2, ls="--", zorder=2); ax.text(0.62, cf.Z0 - 0.065, "transport height 0.30 m (tolerance 3 cm)", fontsize=9, color="#2166ac", ha="right")
ax.add_patch(plt.Rectangle((-0.08, -0.02), 0.16, 0.02, color="0.2")); ax.axhline(0, color="0.5", lw=1); ax.text(-0.25, 0.01, "bench", fontsize=9, color="0.4")
ax.legend(loc="upper right", fontsize=9, title="poses along the RA-MPC run", title_fontsize=9); ax.set_aspect("equal"); ax.set_xlim(-0.3, 0.65); ax.set_ylim(-0.05, 0.85); ax.set_xlabel("horizontal distance from the base axis, along the arm [m]"); ax.set_ylabel("height [m]"); ax.grid(alpha=0.25)
ax.set_title("side view in the arm plane: joints 2, 4, 6 bend in this plane, joint 1 turns the plane", fontsize=10)
ax2.plot(*cf.C.y(Xs.T)[:2], "k--", lw=1.5, alpha=0.6); ax2.plot(0, 0, "ks", ms=8); ax2.set_aspect("equal"); ax2.set_xlim(-0.36, 0.36); ax2.set_ylim(-0.06, 0.66); ax2.grid(alpha=0.25)
ax2.set_xlabel("$y_1$ [m]"); ax2.set_ylabel("$y_2$ [m]"); ax2.set_title("top view of the same poses: all links on the ray to the gripper", fontsize=10)
fig.savefig(os.path.join(cf.FIGS, "franka_side_view.png"), dpi=150); print("saved franka_side_view.png from run", i0)

# MuJoCo rendering from the side (if an OpenGL context is available)
try:
    mid = Xs[idx[len(idx) // 2]]; d.qpos[:] = mid[:NQ]; mujoco.mj_forward(MJ, d); ren = mujoco.Renderer(MJ, 480, 640); cam = mujoco.MjvCamera()
    P = frames(mid); az = np.arctan2(P[-1, 1], P[-1, 0]); cam.lookat[:] = [0.25 * np.cos(az), 0.25 * np.sin(az), 0.35]; cam.distance = 1.6; cam.azimuth = np.degrees(az) - 90; cam.elevation = -5
    ren.update_scene(d, cam); img = ren.render(); plt.imsave(os.path.join(cf.FIGS, "franka_side_render.png"), img); print("saved franka_side_render.png")
except Exception as e: print("MuJoCo rendering not available here:", type(e).__name__, str(e)[:120])
