"""Panda-like arm with simplified (lumped-mass) dynamics for the planar pick-and-place task. Kinematics and meshes of
MuJoCo Menagerie's panda.xml (not modified); joints 3, 5, 7 and the fingers fixed (q3 = q5 = 0, q7 = -pi/4, fingers
0.02 m); joints 1, 2, 4, 6 with armature 0.1, no damping, torque motors; gravity off (compensated, as in Franka's
torque interface). Inertia replaced by point masses: elbow (joint-4 axis) M_E, wrist (joint-6 axis) M_W, tool centre
point M_H, and a yaw inertia I1 of link 1; all other bodies massless. Writes model/panda_simple.xml."""
import os, xml.etree.ElementTree as ET, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, "third_party", "mujoco_menagerie", "franka_emika_panda")
FIXED = {"joint3": 0.0, "joint5": 0.0, "joint7": -np.pi / 4}; FINGER = 0.02; ACTIVE = ["joint1", "joint2", "joint4", "joint6"]


def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2, w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def qrot(q, v):
    q = np.asarray(q, float) / np.linalg.norm(q); w, x, y, z = q
    R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])
    return R @ np.asarray(v, float)
M_E, M_W, M_H, I1, EPS = 6.8, 2.9, 1.5, 0.05, 1e-9
# link 1 turns about z only, so its x, y inertia does not enter the dynamics
LUMP = {"link1": ((0, 0, 0), 1e-3, (I1 / 2, I1 / 2, I1)), "link4": ((0, 0, 0), M_E, (EPS, EPS, EPS)), "link6": ((0, 0, 0), M_W, (EPS, EPS, EPS)), "hand": ((0, 0, 0.1034), M_H, (EPS, EPS, EPS))}

tree = ET.parse(os.path.join(SRC, "panda.xml")); root = tree.getroot()
root.find("compiler").set("meshdir", os.path.join(SRC, "assets"))
opt = root.find("option"); opt.set("integrator", "RK4"); opt.set("timestep", "0.0005"); opt.set("gravity", "0 0 0")
for d in root.iter("default"):
    if d.get("class") == "panda": d.find("joint").set("damping", "0")
for body in list(root.iter("body")):
    for j in list(body.findall("joint")):
        name = j.get("name")
        if name in FIXED:
            axis = np.array([float(a) for a in (j.get("axis") or "0 0 1").split()]); ang = FIXED[name]
            q0 = np.array([float(a) for a in (body.get("quat") or "1 0 0 0").split()]); q0 /= np.linalg.norm(q0)
            body.set("quat", " ".join(f"{v:.12g}" for v in qmul(q0, np.r_[np.cos(ang / 2), np.sin(ang / 2) * axis]))); body.remove(j)
        elif name.startswith("finger_joint"):
            q0 = np.array([float(a) for a in (body.get("quat") or "1 0 0 0").split()]); p0 = np.array([float(a) for a in (body.get("pos") or "0 0 0").split()])
            body.set("pos", " ".join(f"{v:.12g}" for v in p0 + qrot(q0, np.array([0, 1, 0]) * FINGER))); body.remove(j)
    for e in body.findall("inertial"): body.remove(e)
    name = body.get("name")
    if name in LUMP:
        pos, mass, I = LUMP[name]; body.insert(0, ET.Element("inertial", pos=" ".join(map(str, pos)), mass=str(mass), diaginertia=" ".join(map(str, I))))
    elif name != "link0":
        body.insert(0, ET.Element("inertial", pos="0 0 0", mass=str(EPS), diaginertia=f"{EPS} {EPS} {EPS}"))
    else:
        body.insert(0, ET.Element("inertial", pos="0 0 0", mass="1", diaginertia="0.01 0.01 0.01"))
    if name == "hand": ET.SubElement(body, "site", name="tcp", pos="0 0 0.1034", size="0.01", rgba="1 0 0 1")
for tag in ("tendon", "equality", "actuator", "keyframe"):
    for e in root.findall(tag): root.remove(e)
act = ET.SubElement(root, "actuator")
for j in ACTIVE: ET.SubElement(act, "motor", name=f"m_{j}", joint=j, ctrllimited="false")
out = os.path.join(HERE, "model", "panda_simple.xml"); tree.write(out)
import mujoco
m = mujoco.MjModel.from_xml_path(out)
print(f"wrote {out}: nq {m.nq}, masses {dict((m.body(i).name, round(float(m.body_mass[i]), 9)) for i in range(m.nbody) if m.body_mass[i] > 1e-6)}, armature {m.dof_armature}, damping {m.dof_damping}, gravity {m.opt.gravity}")
