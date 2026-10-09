"""3D scene of the kitting workcell for the video and the local viewer (visual only: every added geom has contype =
conaffinity = 0, the red cube is a mocap body, so the arm's dynamics are those of model/panda_simple.xml).
Task coordinates (y1 to the right, y2 forward, metres) map to the world as X = y2, Y = -y1. Writes model/scene_red.xml."""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cert_franka import CertFranka, ALPHA, Z0, SPEC, TARGET, RT
HERE = os.path.dirname(os.path.abspath(__file__))
X0, YS, SIG = ALPHA * SPEC["X0"], [ALPHA * v for v in SPEC["YS"]], SPEC["SIGMA"]; RD = ALPHA * SPEC.get("rack_d", 1.2); FX = SPEC["fixture"]
COL = ["0.84 0.15 0.16", "0.12 0.47 0.71", "0.17 0.63 0.17", "1.00 0.50 0.05"]; H = 0.02; PLAT = Z0 - H; NOC = 'contype="0" conaffinity="0"'


def w(y1, y2): return y2, -y1                                     # task (y1, y2) -> world (X, Y)


def box(name, y1r, y2r, zr, rgba, extra=""):
    (a, b), (c, d), (e, f) = sorted(y1r), sorted(y2r), zr; X = (c + d) / 2; Y = -(a + b) / 2
    return f'<geom name="{name}" type="box" pos="{X:.4f} {Y:.4f} {(e + f) / 2:.4f}" size="{(d - c) / 2:.4f} {(b - a) / 2:.4f} {(f - e) / 2:.4f}" rgba="{rgba}" {NOC} {extra}/>'


g = []
for s, nm in ((-1, "left"), (1, "right")):                         # trays at the transport height: platform, legs, dividers, back wall
    y1a, y1b = s * (X0 - 0.03), s * (X0 + RD); y2a, y2b = YS[0] - 0.05, YS[-1] + 0.05
    g.append(box(f"{nm}_platform", (y1a, y1b), (y2a, y2b), (PLAT - 0.012, PLAT), "0.80 0.74 0.64 1"))
    for i, (u, v) in enumerate(((y1a, y2a), (y1a, y2b), (y1b, y2a), (y1b, y2b))):
        g.append(box(f"{nm}_leg{i}", (u - 0.006, u + 0.006), (v - 0.006, v + 0.006), (0, PLAT - 0.012), "0.55 0.50 0.42 1"))
    for i, yy in enumerate([y2a] + [(YS[j] + YS[j + 1]) / 2 for j in range(3)] + [y2b]):
        g.append(box(f"{nm}_div{i}", (y1a, y1b), (yy - 0.002, yy + 0.002), (PLAT, PLAT + 0.035), "0.62 0.56 0.46 1"))
    g.append(box(f"{nm}_back", (y1b - 0.003 * s, y1b), (y2a, y2b), (PLAT, PLAT + 0.05), "0.62 0.56 0.46 1"))
for k in (1, 2, 3):                                               # objects in the left tray (red is the mocap body below)
    g.append(box(f"obj{k}", (-X0 - H, -X0 + H), (YS[k] - H, YS[k] + H), (PLAT, PLAT + 2 * H), COL[k] + " 1"))
j_red = SIG[0]
for j in range(4):                                                # slot markers in the right tray (dashed look: four short bars)
    rgba = (COL[0] + " 1") if j == j_red else "0.35 0.35 0.35 1"
    for t, (u0, u1, v0, v1) in enumerate(((X0 - H, X0 + H, YS[j] - H, YS[j] - H + 0.003), (X0 - H, X0 + H, YS[j] + H - 0.003, YS[j] + H), (X0 - H, X0 - H + 0.003, YS[j] - H, YS[j] + H), (X0 + H - 0.003, X0 + H, YS[j] - H, YS[j] + H))):
        g.append(box(f"slot{j}_{t}", (u0, u1), (v0, v1), (PLAT, PLAT + 0.002), rgba))
g.append(box("obstacle", (-ALPHA * FX["w"] / 2, ALPHA * FX["w"] / 2), (ALPHA * FX["y0"], ALPHA * FX["top"]), (0, 0.42), "0.72 0.72 0.74 1"))
# corridor edge (planar section of X_S at the transport height) and the target disc, translucent
c = CertFranka(); g1 = np.linspace(-0.36, 0.36, 600); g2 = np.linspace(0.1, 0.6, 420); A, B = np.meshgrid(g1, g2)
P = c.psi_of_y(np.vstack([A.ravel(), B.ravel(), Z0 * np.ones(A.size), np.zeros(A.size)])).reshape(A.shape)
cs = plt.contour(A, B, P, levels=[0]); segs = cs.allsegs[0]; plt.close()
n = 0
for seg in segs:
    seg = seg[::3]
    for p, q in zip(seg[:-1], seg[1:]):
        (Xa, Ya), (Xb, Yb) = w(*p), w(*q); g.append(f'<geom name="edge{n}" type="capsule" fromto="{Xa:.4f} {Ya:.4f} {Z0:.4f} {Xb:.4f} {Yb:.4f} {Z0:.4f}" size="0.0022" rgba="0.13 0.40 0.67 0.55" {NOC} group="1"/>'); n += 1
Xt, Yt = w(*TARGET); g.append(f'<geom name="target" type="cylinder" pos="{Xt:.4f} {Yt:.4f} {Z0:.4f}" size="{RT:.4f} 0.0015" rgba="0.84 0.38 0.30 0.35" {NOC} group="1"/>')
Xs, Ys = w(*(ALPHA * np.array(SPEC["source"])))
xml = f"""<mujoco model="kitting red task">
  <include file="panda_simple.xml"/>
  <visual>
    <global offwidth="1280" offheight="720"/>
    <quality shadowsize="4096"/>
    <headlight ambient="0.38 0.38 0.38" diffuse="0.45 0.45 0.45" specular="0.1 0.1 0.1"/>
  </visual>
  <asset>
    <texture name="skybox" type="skybox" builtin="gradient" rgb1="1 1 1" rgb2="0.86 0.89 0.93" width="512" height="3072"/>
    <texture name="bench" type="2d" builtin="checker" rgb1="0.94 0.92 0.88" rgb2="0.91 0.89 0.85" width="512" height="512" mark="none"/>
    <material name="bench" texture="bench" texrepeat="10 10" reflectance="0.04"/>
  </asset>
  <worldbody>
    <light name="key" pos="0.6 0.8 1.8" dir="-0.3 -0.4 -1" directional="false" castshadow="true" diffuse="0.55 0.55 0.55"/>
    <geom name="bench" type="plane" size="1.2 1.2 0.01" pos="0.25 0 0" material="bench" {NOC}/>
    <camera name="top" pos="0.30 0 3.0" xyaxes="0 -1 0 1 0 0" fovy="14"/>
    {chr(10).join('    ' + s for s in g)}
    <body name="red_cube" mocap="true" pos="{Xs:.4f} {Ys:.4f} {Z0:.4f}">
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
print(f"wrote model/scene_red.xml: {m.ngeom} geoms ({n} corridor-edge segments), nq {m.nq}, mocap {m.nmocap}; arm mass matrix identical to panda_simple: {np.allclose(M1, M2)}, contacts at a test pose {d1.ncon}")
