"""Build the reduced Panda model for the planar pick-and-place task from MuJoCo Menagerie's panda.xml (not modified):
joints 3, 5, 7 and the fingers are removed with their angles baked into the body frames (q3 = q5 = 0, q7 = -pi/4,
fingers at 0.02 m), joints 1, 2, 4, 6 keep armature 0.1 with zero damping and get torque motors without limits,
gravity is off (Franka's torque interface compensates gravity), integrator RK4 at 0.5 ms; a site 'tcp' 0.1034 m below
the flange. Writes model/panda_task.xml."""
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


tree = ET.parse(os.path.join(SRC, "panda.xml")); root = tree.getroot()
root.find("compiler").set("meshdir", os.path.join(SRC, "assets"))
opt = root.find("option"); opt.set("integrator", "RK4"); opt.set("timestep", "0.0005"); opt.set("gravity", "0 0 0")
for d in root.iter("default"):
    if d.get("class") == "panda": d.find("joint").set("damping", "0")
parent = {c: p for p in root.iter() for c in p}
for body in list(root.iter("body")):
    for j in list(body.findall("joint")):
        name = j.get("name")
        if name in FIXED:
            axis = np.array([float(a) for a in (j.get("axis") or "0 0 1").split()]); ang = FIXED[name]
            q0 = np.array([float(a) for a in (body.get("quat") or "1 0 0 0").split()]); q0 /= np.linalg.norm(q0)
            qa = np.r_[np.cos(ang / 2), np.sin(ang / 2) * axis]; body.set("quat", " ".join(f"{v:.12g}" for v in qmul(q0, qa))); body.remove(j)
        elif name.startswith("finger_joint"):
            q0 = np.array([float(a) for a in (body.get("quat") or "1 0 0 0").split()]); p0 = np.array([float(a) for a in (body.get("pos") or "0 0 0").split()])
            body.set("pos", " ".join(f"{v:.12g}" for v in p0 + qrot(q0, np.array([0, 1, 0]) * FINGER))); body.remove(j)
    if body.get("name") == "hand": ET.SubElement(body, "site", name="tcp", pos="0 0 0.1034", size="0.01", rgba="1 0 0 1")
for tag in ("tendon", "equality", "actuator", "keyframe"):
    for e in root.findall(tag): root.remove(e)
act = ET.SubElement(root, "actuator")
for j in ACTIVE: ET.SubElement(act, "motor", name=f"m_{j}", joint=j, ctrllimited="false")
os.makedirs(os.path.join(HERE, "model"), exist_ok=True); out = os.path.join(HERE, "model", "panda_task.xml"); tree.write(out)
import mujoco
m = mujoco.MjModel.from_xml_path(out); d = mujoco.MjData(m)
print(f"wrote {out}: nq {m.nq}, nv {m.nv}, nu {m.nu}, joints {[m.joint(i).name for i in range(m.njnt)]}, armature {m.dof_armature}, damping {m.dof_damping}, gravity {m.opt.gravity}, integrator {m.opt.integrator}, timestep {m.opt.timestep}")
print("joint ranges:", {m.joint(i).name: np.round(m.jnt_range[i], 4).tolist() for i in range(m.njnt)})
