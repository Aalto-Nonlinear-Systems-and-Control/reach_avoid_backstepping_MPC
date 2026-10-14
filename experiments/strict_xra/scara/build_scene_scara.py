"""MuJoCo models of the SCARA video. model/scara_arm.xml: the arm alone (static base, link 1 and link 2 on vertical hinges,
torque motors, point-like link inertias from ../arm_top scaled by scara_common, RK4 at 0.5 ms as the paper's plant).
model/scene_scara.xml: the arm in the workcell of the Franka scene (Menagerie floor and sky, one table holding the arm,
the cubes and a 0.30 m obstacle, slots and cube homes as dashed marks, edge of X_S and ring of X_T at the transport height
in geom group 1) plus mocap bodies for the moving parts that do not enter the planar dynamics: the Z-axis shaft with the
gripper flange, the two fingers and the red cube. Checks the mass matrix against C_TAU M(q) and the tip against the model."""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scara_common import *
import config_top as cf
from cert_top import fk, mass, M1, M2, L1, L2, LC1, LC2, I1, I2, SPEC, TARGET_POS, RT
NOC = 'contype="0" conaffinity="0"'; TH, LEG = 0.04, 0.75; OBST_H = 0.30; COL = ["0.80 0.12 0.12", "0.10 0.35 0.75", "0.12 0.58 0.18", "0.95 0.50 0.05"]
l1, l2, c1, c2 = ALPHA * L1, ALPHA * L2, ALPHA * LC1, ALPHA * LC2; m1, m2 = MSCALE * M1, MSCALE * M2; i1, i2 = C_TAU * I1, C_TAU * I2
H = H_CUBE

arm = f"""<mujoco model="scara arm (two-link example of the paper, scaled)">
  <compiler angle="radian" inertiafromgeom="false"/>
  <option timestep="0.0005" integrator="RK4"/>
  <default><geom {NOC}/><joint damping="0" armature="0" frictionloss="0"/></default>
  <asset>
    <material name="scara_white" rgba="0.92 0.92 0.91 1" specular="0.3" shininess="0.5"/>
    <material name="scara_dark" rgba="0.20 0.21 0.23 1" specular="0.2" shininess="0.3"/>
    <material name="scara_metal" rgba="0.62 0.64 0.67 1" specular="0.6" shininess="0.8"/>
  </asset>
  <worldbody>
    <body name="base" quat="0.7071067811865476 0 0 -0.7071067811865476">
      <geom type="cylinder" size="0.10 {(SHOULDER_Z - 0.06) / 2:.4f}" pos="0 0 {(SHOULDER_Z - 0.06) / 2:.4f}" material="scara_dark"/>
      <geom type="box" size="0.14 0.14 0.012" pos="0 0 0.012" material="scara_dark"/>
      <geom type="cylinder" size="0.115 0.03" pos="0 0 {SHOULDER_Z - 0.03:.4f}" material="scara_white"/>
      <body name="link1" pos="0 0 {SHOULDER_Z}">
        <inertial pos="{c1} 0 0" mass="{m1}" diaginertia="{i1} {i1} {i1}"/>
        <joint name="j1" type="hinge" axis="0 0 1"/>
        <geom type="box" size="{l1 / 2} 0.065 0.035" pos="{l1 / 2} 0 0.035" material="scara_white"/>
        <geom type="cylinder" size="0.065 0.035" pos="0 0 0.035" material="scara_white"/>
        <geom type="cylinder" size="0.065 0.035" pos="{l1} 0 0.035" material="scara_white"/>
        <geom type="cylinder" size="0.045 0.004" pos="0 0 0.074" material="scara_dark"/>
        <body name="link2" pos="{l1} 0 {LINK2_DZ}">
          <inertial pos="{c2} 0 0" mass="{m2}" diaginertia="{i2} {i2} {i2}"/>
          <joint name="j2" type="hinge" axis="0 0 1"/>
          <geom type="cylinder" size="0.055 0.025" pos="0 0 0.025" material="scara_dark"/>
          <geom type="box" size="{l2 / 2} 0.055 0.03" pos="{l2 / 2} 0 0.08" material="scara_white"/>
          <geom type="cylinder" size="0.055 0.03" pos="0 0 0.08" material="scara_white"/>
          <geom type="cylinder" size="0.055 0.03" pos="{l2} 0 0.08" material="scara_white"/>
          <geom type="cylinder" size="0.045 0.11" pos="{l2} 0 0.13" material="scara_white"/>
          <geom type="cylinder" size="0.035 0.006" pos="{l2} 0 0.246" material="scara_dark"/>
          <site name="tip" pos="{l2} 0 0" size="0.005"/>
        </body>
      </body>
    </body>
  </worldbody>
  <actuator><motor name="m1" joint="j1" ctrllimited="false"/><motor name="m2" joint="j2" ctrllimited="false"/></actuator>
</mujoco>
"""
open(os.path.join(HERE, "model", "scara_arm.xml"), "w").write(arm)


def box(name, y1r, y2r, zr, rgba=None, mat=None):
    (a, b), (c, d), (e, f) = sorted(y1r), sorted(y2r), zr; look = f'material="{mat}"' if mat else f'rgba="{rgba}"'
    return f'<geom name="{name}" type="box" pos="{(c + d) / 2:.4f} {-(a + b) / 2:.4f} {(e + f) / 2:.4f}" size="{(d - c) / 2:.4f} {(b - a) / 2:.4f} {(f - e) / 2:.4f}" {look} {NOC}/>'


def dashed_square(name, cx, cy, half, rgba, n=3, w=0.0025, z=0.0006):
    out = []; seg = 2 * half / (2 * n - 1)
    for side in range(4):
        for i in range(n):
            s0 = -half + 2 * i * seg; s1 = s0 + seg
            r = [((cx + s0, cx + s1), (cy - half, cy - half + w)), ((cx + s0, cx + s1), (cy + half - w, cy + half)), ((cx - half, cx - half + w), (cy + s0, cy + s1)), ((cx + half - w, cx + half), (cy + s0, cy + s1))][side]
            out.append(box(f"{name}_{side}_{i}", r[0], r[1], (0, z), rgba))
    return out


X0, YS, SIG, FX = ALPHA * SPEC["X0"], [ALPHA * v for v in SPEC["YS"]], SPEC["SIGMA"], SPEC["fixture"]; TGT = ALPHA * TARGET_POS; R_T = ALPHA * RT
g = [box("tabletop", (-0.62, 0.62), (-0.32, 0.82), (-TH, 0), mat="table")]
for i, (u, v) in enumerate(((-0.58, -0.28), (-0.58, 0.78), (0.58, -0.28), (0.58, 0.78))): g.append(box(f"leg{i}", (u - 0.025, u + 0.025), (v - 0.025, v + 0.025), (-LEG, -TH), "0.22 0.22 0.24 1"))
for k in (1, 2, 3): g.append(box(f"obj{k}", (-X0 - H, -X0 + H), (YS[k] - H, YS[k] + H), (0, 2 * H), COL[k] + " 1"))
for j in range(4):
    if j != SIG[0]: g += dashed_square(f"slot{j}", X0, YS[j], H + 0.006, "0.25 0.25 0.25 1")
g += dashed_square("slot_red", TGT[0], TGT[1], H + 0.006, COL[0] + " 1")
for k in range(4): g += dashed_square(f"home{k}", -X0, YS[k], H + 0.006, "0.55 0.55 0.55 1")
g.append(box("obstacle", (-ALPHA * FX["w"] / 2, ALPHA * FX["w"] / 2), (ALPHA * FX["y0"], ALPHA * FX["top"]), (0, OBST_H), "0.30 0.34 0.40 1"))
G1, G2 = np.meshgrid(np.linspace(-3.6, 3.6, 600), np.linspace(1.0, 6.0, 420)); P = np.broadcast_to(cf.C0.psi(G1, G2), G1.shape)
cs = plt.contour(ALPHA * G1, ALPHA * G2, P, levels=[0]); segs = cs.allsegs[0]; plt.close(); n = 0
for seg in segs:
    seg = seg[::3]
    for p, q in zip(seg[:-1], seg[1:]):
        g.append(f'<geom name="edge{n}" type="capsule" fromto="{p[1]:.4f} {-p[0]:.4f} {Z0} {q[1]:.4f} {-q[0]:.4f} {Z0}" size="0.0022" rgba="0.20 0.55 0.95 1" {NOC} group="1"/>'); n += 1
for i, a in enumerate(np.linspace(0, 2 * np.pi, 41)[:-1]):
    p, q = TGT + R_T * np.array([np.cos(a), np.sin(a)]), TGT + R_T * np.array([np.cos(a + np.pi / 20), np.sin(a + np.pi / 20)])
    g.append(f'<geom name="target{i}" type="capsule" fromto="{p[1]:.4f} {-p[0]:.4f} {Z0} {q[1]:.4f} {-q[0]:.4f} {Z0}" size="0.0022" rgba="0.90 0.25 0.20 1" {NOC} group="1"/>')
SRC = to_world(SPEC["source"]); FH = 0.045                      # finger half length (z)
scene = f"""<mujoco model="kitting red task, SCARA">
  <include file="scara_arm.xml"/>
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
    <body name="quill" mocap="true" pos="{SRC[0]:.4f} {SRC[1]:.4f} {Z0}">
      <geom name="shaft" type="cylinder" size="0.012 0.375" pos="0 0 0.475" material="scara_metal" {NOC}/>
      <geom name="flange" type="box" size="0.022 0.045 0.008" pos="0 0 0.092" material="scara_dark" {NOC}/>
      <geom name="wrist" type="cylinder" size="0.022 0.012" pos="0 0 0.112" material="scara_dark" {NOC}/>
    </body>
    <body name="finger_l" mocap="true" pos="{SRC[0]:.4f} {SRC[1] + GAP_OPEN + 0.005:.4f} {Z0 - 0.012 + FH:.4f}">
      <geom name="finger_l" type="box" size="0.010 0.005 {FH}" material="scara_dark" {NOC}/>
    </body>
    <body name="finger_r" mocap="true" pos="{SRC[0]:.4f} {SRC[1] - GAP_OPEN - 0.005:.4f} {Z0 - 0.012 + FH:.4f}">
      <geom name="finger_r" type="box" size="0.010 0.005 {FH}" material="scara_dark" {NOC}/>
    </body>
    <body name="red_cube" mocap="true" pos="{SRC[0]:.4f} {SRC[1]:.4f} {H}">
      <geom name="red_cube" type="box" size="{H} {H} {H}" rgba="{COL[0]} 1" {NOC}/>
    </body>
  </worldbody>
</mujoco>
"""
open(os.path.join(HERE, "model", "scene_scara.xml"), "w").write(scene)

# kitting cycle scene: all four cubes are mocap bodies at their homes, slots dashed in the colour of the cube assigned to
# them, homes dashed in a lighter shade; the corridor of the current task is drawn by the renderer (user geoms)
pastel = lambda c: " ".join(f"{0.45 * float(v) + 0.55:.3f}" for v in c.split())
gc = g[:5] + [s for s in g if 'name="obstacle"' in s]
for j in range(4): gc += dashed_square(f"slot{j}", X0, YS[j], H + 0.006, COL[SIG.index(j)] + " 1")
for k in range(4): gc += dashed_square(f"home{k}", -X0, YS[k], H + 0.006, pastel(COL[k]) + " 1")
cubes = "\n".join(f'''    <body name="cube{k}" mocap="true" pos="{ALPHA * SPEC["YS"][k]:.4f} {X0:.4f} {H}">
      <geom name="cube{k}" type="box" size="{H} {H} {H}" rgba="{COL[k]} 1" {NOC}/>
    </body>''' for k in range(4))
i0, i1 = scene.index("    <body name=\"red_cube\""), scene.index("  </worldbody>")
i_g0 = scene.index("    <camera name=\"top\"") + len(scene[scene.index("    <camera name=\"top\""):].split("\n")[0]) + 1; i_q = scene.index("    <body name=\"quill\"")
cyc = scene[:i_g0] + "\n".join("    " + s for s in gc) + "\n" + scene[i_q:i0] + cubes + "\n" + scene[i1:]
open(os.path.join(HERE, "model", "scene_scara_cycle.xml"), "w").write(cyc.replace('model="kitting red task, SCARA"', 'model="kitting cycle, SCARA"'))
import mujoco
ma = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scara_arm.xml")); ms = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_scara.xml"))
da, ds = mujoco.MjData(ma), mujoco.MjData(ms); r = np.random.default_rng(0); em, et = 0.0, 0.0
for _ in range(200):
    q = r.uniform(-np.pi, np.pi, 2); da.qpos[:] = q; mujoco.mj_forward(ma, da); Mm = np.zeros((2, 2)); mujoco.mj_fullM(ma, da, Mm)
    em = max(em, np.abs(Mm - C_TAU * mass(q.reshape(2, 1))[:, :, 0]).max()); et = max(et, np.abs(da.site_xpos[0][:2] - to_world(fk(q))).max())
mc = mujoco.MjModel.from_xml_path(os.path.join(HERE, "model", "scene_scara_cycle.xml")); print(f"cycle scene: {mc.ngeom} geoms, nq {mc.nq}, nmocap {mc.nmocap}")
print(f"wrote model/scara_arm.xml and model/scene_scara.xml ({ms.ngeom} geoms, nq {ms.nq}, nmocap {ms.nmocap}); max |M_mujoco - C_TAU M| {em:.1e}, max tip error {et:.1e} m; C_TAU {C_TAU:g}")
