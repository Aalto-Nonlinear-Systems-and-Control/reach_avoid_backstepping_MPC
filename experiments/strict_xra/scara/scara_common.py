"""Shared constants and helpers of the SCARA video. The horizontal two-link arm of ../arm_top (task red_fwd_sq) is the
planar part of a SCARA: both hinges turn about vertical axes, so gravity does not enter its dynamics. Lengths are scaled
to metres by ALPHA = 0.1 (links 0.4 m) and the link masses by MSCALE = 10, so M'(q) = C_TAU M(q) with C_TAU = MSCALE ALPHA^2
and the torques scale by C_TAU; the joint trajectories are exactly those of the paper's model. Task coordinates (y1 right,
y2 forward) map to the world as X = ALPHA y2, Y = -ALPHA y1 (a rotation, so joint angles and torques keep their sign)."""
import os, sys, numpy as np
os.environ.setdefault("TASK", "red_fwd_sq"); os.environ.setdefault("LAM", "0.1")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "arm_top"))
ALPHA, MSCALE = 0.1, 10.0; C_TAU = MSCALE * ALPHA**2
Z0, H_CUBE = 0.30, 0.02                     # transport height of the tool point (cube centre), half size of the cubes
Z_LOW = H_CUBE                              # tool point when the cube rests on the table
SHOULDER_Z, LINK2_DZ = 0.50, 0.07           # height of the shoulder hinge, rise of the elbow hinge above it
GAP_OPEN, GAP_CLOSED, T_FINGER = 0.032, 0.020, 0.15   # inner half gap of the fingers, closing time


def to_world(y):
    """Task coordinates (2, ...) -> world (X, Y)."""
    y = np.asarray(y, float); return np.array([ALPHA * y[1], -ALPHA * y[0]])


def minjerk(s):
    s = np.clip(s, 0.0, 1.0); return 10 * s**3 - 15 * s**4 + 6 * s**5


def set_moving(m, d, q, z, gap, cube):
    """Arm joints q, tool point height z, inner half gap of the fingers, cube centre (world); the gripper keeps a fixed
    world orientation (the wrist axis compensates the yaw of link 2), fingers close along the world Y axis."""
    import mujoco
    d.qpos[:2] = q; mujoco.mj_forward(m, d); x, y = d.site_xpos[m.site("tip").id][:2]
    d.mocap_pos[m.body("quill").mocapid[0]] = (x, y, z)
    for name, sgn in (("finger_l", 1), ("finger_r", -1)): d.mocap_pos[m.body(name).mocapid[0]] = (x, y + sgn * (gap + 0.005), z - 0.012 + 0.045)
    d.mocap_pos[m.body("red_cube").mocapid[0]] = cube; mujoco.mj_forward(m, d); return np.array([x, y, z])
