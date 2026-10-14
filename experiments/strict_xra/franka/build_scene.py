"""3D scene of the kitting workcell for the video and the local viewer, in the style of MuJoCo Menagerie's scene.xml
(blue checker floor, gradient skybox, moderate headlight, directional light). One table: the robot base stands on the
tabletop (z = 0), the cubes and the obstacle stand on the table, the slots, the corridor edge and the target are marks on
the tabletop. Visual only (contype = conaffinity = 0, the red cube is a mocap body): the arm's dynamics are those of
model/panda_simple.xml. Task coordinates (y1 right, y2 forward) map to the world as X = y2, Y = -y1. Writes model/scene_red.xml."""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cert_franka import CertFranka, ALPHA, Z0, SPEC, TARGET, RT
HERE = os.path.dirname(os.path.abspath(__file__))
X0, YS, SIG, FX = ALPHA * SPEC["X0"], [ALPHA * v for v in SPEC["YS"]], SPEC["SIGMA"], SPEC["fixture"]
COL = ["0.80 0.12 0.12", "0.10 0.35 0.75", "0.12 0.58 0.18", "0.95 0.50 0.05"]; H = 0.02; NOC = 'contype="0" conaffinity="0"'; TH, LEG = 0.04, 0.75; OBST_H = float(os.environ.get("OBST_H", 0.30))             # top at the transport height: blocks the carried cube, the hand passes above (16.5 mm)


def box(name, y1r, y2r, zr, rgba=None, mat=None):
    (a, b), (c, d), (e, f) = sorted(y1r), sorted(y2r), zr; look = f'material="{mat}"' if mat else f'rgba="{rgba}"'
    return f'<geom name="{name}" type="box" pos="{(c + d) / 2:.4f} {-(a + b) / 2:.4f} {(e + f) / 2:.4f}" size="{(d - c) / 2:.4f} {(b - a) / 2:.4f} {(f - e) / 2:.4f}" {look} {NOC}/>'


def dashed_square(name, cx, cy, half, rgba, n=3, w=0.0025, z=0.0006):
    """Dashed outline on the tabletop (task coordinates of the centre)."""
    out = []; L = 2 * half; seg = L / (2 * n - 1)
    for side in range(4):
        for i in range(n):
            s0 = -half + 2 * i * seg; s1 = s0 + seg
            if side == 0: r = ((cx + s0, cx + s1), (cy - half, cy - half + w))
            elif side == 1: r = ((cx + s0, cx + s1), (cy + half - w, cy + half))
            elif side == 2: r = ((cx - half, cx - half + w), (cy + s0, cy + s1))
            else: r = ((cx + half - w, cx + half), (cy + s0, cy + s1))
            out.append(box(f"{name}_{side}_{i}", r[0], r[1], (0, z), rgba))
    return out


g = [box("tabletop", (-0.62, 0.62), (-0.32, 0.82), (-TH, 0), mat="table")]
for i, (u, v) in enumerate(((-0.58, -0.28), (-0.58, 0.78), (0.58, -0.28), (0.58, 0.78))): g.append(box(f"leg{i}", (u - 0.025, u + 0.025), (v - 0.025, v + 0.025), (-LEG, -TH), "0.22 0.22 0.24 1"))
for k in (1, 2, 3): g.append(box(f"obj{k}", (-X0 - H, -X0 + H), (YS[k] - H, YS[k] + H), (0, 2 * H), COL[k] + " 1"))
for j in range(4):
    if j != SIG[0]: g += dashed_square(f"slot{j}", X0, YS[j], H + 0.006, "0.25 0.25 0.25 1")
g += dashed_square("slot_red", TARGET[0], TARGET[1], H + 0.006, COL[0] + " 1")
for k in range(4): g += dashed_square(f"home{k}", -X0, YS[k], H + 0.006, "0.55 0.55 0.55 1")
g.append(box("obstacle", (-ALPHA * FX["w"] / 2, ALPHA * FX["w"] / 2), (ALPHA * FX["y0"], ALPHA * FX["top"]), (0, OBST_H), "0.30 0.34 0.40 1"))
c = CertFranka(); g1 = np.linspace(-0.36, 0.36, 600); g2 = np.linspace(0.1, 0.6, 420); A, B = np.meshgrid(g1, g2)
P = c.psi_of_y(np.vstack([A.ravel(), B.ravel(), Z0 * np.ones(A.size), np.zeros(A.size)])).reshape(A.shape); cs = plt.contour(A, B, P, levels=[0]); segs = cs.allsegs[0]; plt.close(); n = 0
for seg in segs:
    seg = seg[::3]
    for p, q in zip(seg[:-1], seg[1:]):
        g.append(f'<geom name="edge{n}" type="capsule" fromto="{p[1]:.4f} {-p[0]:.4f} {Z0} {q[1]:.4f} {-q[0]:.4f} {Z0}" size="0.0022" rgba="0.20 0.55 0.95 1" {NOC} group="1"/>'); n += 1
for i, a in enumerate(np.linspace(0, 2 * np.pi, 41)[:-1]):
    p, q = TARGET + RT * np.array([np.cos(a), np.sin(a)]), TARGET + RT * np.array([np.cos(a + np.pi / 20), np.sin(a + np.pi / 20)])
    g.append(f'<geom name="target{i}" type="capsule" fromto="{p[1]:.4f} {-p[0]:.4f} {Z0} {q[1]:.4f} {-q[0]:.4f} {Z0}" size="0.0022" rgba="0.90 0.25 0.20 1" {NOC} group="1"/>')
Xs, Ys = ALPHA * SPEC["source"][1], -ALPHA * SPEC["source"][0]
xml = f"""<mujoco model="kitting red task">
  <include file="panda_simple.xml"/>
  <visual>
    <global offwidth="1280" offheight="720" azimuth="120" elevation="-20"/>
    <headlight diffuse="0.40 0.40 0.40" ambient="0.30 0.30 0.30" specular="0 0 0"/>
    <rgba haze="0.15 0.25 0.35 1"/>
    <quality shadowsize="4096"/>
  </visual>
  <asset>
    <texture type="skybox" builtin="gradient" rgb1="0.3 0.5 0.7" rgb2="0 0 0" width="512" height="3072"/>
    <texture type="2d" name="groundplane" builtin="checker" mark="edge" rgb1="0.2 0.3 0.4" rgb2="0.1 0.2 0.3" markrgb="0.8 0.8 0.8" width="300" height="300"/>
    <material name="groundplane" texture="groundplane" texuniform="true" texrepeat="5 5" reflectance="0.2"/>
    <material name="table" rgba="0.58 0.50 0.40 1" reflectance="0.03" specular="0" shininess="0.2"/>
  </asset>
  <worldbody>
    <light pos="0.3 0 1.8" dir="0 0 -1" directional="true" diffuse="0.65 0.65 0.65" specular="0.10 0.10 0.10" castshadow="true"/>
    <geom name="floor" size="0 0 0.05" type="plane" pos="0 0 {-LEG}" material="groundplane" {NOC}/>
    <camera name="top" pos="0.30 0 3.0" xyaxes="0 -1 0 1 0 0" fovy="14"/>
    {chr(10).join('    ' + s for s in g)}
    <body name="red_cube" mocap="true" pos="{Xs:.4f} {Ys:.4f} {H}">
      <geom type="box" size="{H} {H} {H}" rgba="{COL[0]} 1" {NOC}/>
    </body>
  </worldbody>
</mujoco>
"""
open(os.path.join(HERE, "model", "scene_red.xml"), "w").write(xml)
import mujoco
m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_red.xml")); mr = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "panda_simple.xml"))
d1, d2 = mujoco.MjData(m), mujoco.MjData(mr); q = np.array([0.3, -0.4, -2.2, 1.9]); d1.qpos[:] = q; d2.qpos[:] = q; mujoco.mj_forward(m, d1); mujoco.mj_forward(mr, d2)
M1, M2 = np.zeros((4, 4)), np.zeros((4, 4)); mujoco.mj_fullM(m, d1, M1); mujoco.mj_fullM(mr, d2, M2)
print(f"wrote model/scene_red.xml: {m.ngeom} geoms, arm mass matrix identical to panda_simple: {np.allclose(M1, M2)}, contacts {d1.ncon}")
