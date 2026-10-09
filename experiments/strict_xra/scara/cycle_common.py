"""State of the kitting-cycle scene (model/scene_scara_cycle.xml) at time t from data/cycle.npz: arm joints, Z axis,
fingers, the four cubes (at home, carried, at the slot, carried back, at home), the active task and its phase; corridor
edge (zero level of the task's psi, drawn at the transport height) and target ring of each task, added as user geoms."""
import os, json, numpy as np, sympy as sp, matplotlib, mujoco
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scara_common import *
from cert_top import fk
D = np.load(os.path.join(HERE, "data", "cycle.npz")); Q, Z, GAP, PH, TK, EV = D["q"], D["z"], D["gap"], D["phase"], D["task"], D["events"]
TASKS = [str(t) for t in D["tasks"]]; INFO = json.loads(str(D["info"])); DT = 0.0005; TEND = (len(Q) - 1) * DT; UMAX_W = D["umax"]
AT = os.path.join(HERE, "..", "arm_top"); SPECS = {t: json.load(open(os.path.join(AT, "tasks", f"{t}.json"))) for t in TASKS}
NAMES = ("red", "blue", "green", "orange"); S0 = SPECS[TASKS[0]]; HOME = [np.r_[to_world([-S0["X0"], S0["YS"][k]]), H_CUBE] for k in range(4)]
PHASE = {0: "pick (Z axis; arm holds)", 1: "transport by RA-MPC (certified reach-avoid)", 2: "place (not certified)", 3: "transit to the next cube (not certified)"}
idx = lambda ts: min(max(int(round(ts / DT)), 0), len(Q) - 1)
tip = lambda ts: to_world(fk(Q[idx(ts)]))


def label(i):
    n, d = TASKS[i].split("_")[:2]; return f"task {i + 1}/8: {n} cube {'to its slot' if d == 'fwd' else 'back to the rack'}"


def cube_pos(k, ts):
    rest = HOME[k]
    for kk, t_att, t_rel, _ in EV:
        if int(kk) != k: continue
        if ts < t_att: return rest
        if ts < t_rel: return np.r_[tip(ts), Z[idx(ts)]]
        rest = np.r_[tip(t_rel), H_CUBE]
    return rest


def set_state(m, d, ts):
    i = idx(ts); d.qpos[:2] = Q[i]; mujoco.mj_forward(m, d); x, y = d.site_xpos[m.site("tip").id][:2]; z, g = Z[i], GAP[i]
    d.mocap_pos[m.body("quill").mocapid[0]] = (x, y, z)
    for name, sgn in (("finger_l", 1), ("finger_r", -1)): d.mocap_pos[m.body(name).mocapid[0]] = (x, y + sgn * (g + 0.005), z - 0.012 + 0.045)
    for k in range(4): d.mocap_pos[m.body(f"cube{k}").mocapid[0]] = cube_pos(k, ts)
    mujoco.mj_forward(m, d); return int(TK[i]), int(PH[i])


def _edges(task):
    S = SPECS[task]; f = sp.lambdify(sp.symbols("y1 y2"), sp.sympify(S["psi"]), "numpy")
    G1, G2 = np.meshgrid(np.linspace(-3.6, 3.6, 600), np.linspace(1.0, 6.0, 420)); cs = plt.contour(ALPHA * G1, ALPHA * G2, np.broadcast_to(f(G1, G2), G1.shape), levels=[0]); segs = cs.allsegs[0]; plt.close()
    out = []
    for seg in segs:
        seg = seg[::3]; out += [((p[1], -p[0], Z0), (q[1], -q[0], Z0)) for p, q in zip(seg[:-1], seg[1:])]
    c = to_world(S["target"]); r = ALPHA * S["RT"]; a = np.linspace(0, 2 * np.pi, 41)
    ring = [((c[0] + r * np.cos(u), c[1] + r * np.sin(u), Z0), (c[0] + r * np.cos(v), c[1] + r * np.sin(v), Z0)) for u, v in zip(a[:-1], a[1:])]
    return out, ring


EDGES = {t: _edges(t) for t in TASKS}
TRAILS = {i: [(ts, np.r_[tip(ts), Z0]) for ts in np.arange(inf["t_ra"], inf["t_entry"] + 1e-9, 0.01)] for i, inf in enumerate(INFO)}


def add_marks(scn, ts, i):
    """Corridor edge and target ring of task i, and the trail of its transport up to ts."""
    def cap(a, b, rgba, w=0.0022):
        if scn.ngeom >= scn.maxgeom: return
        g = scn.geoms[scn.ngeom]; mujoco.mjv_initGeom(g, mujoco.mjtGeom.mjGEOM_CAPSULE, np.zeros(3), np.zeros(3), np.eye(3).ravel(), np.array(rgba, float))
        mujoco.mjv_connector(g, mujoco.mjtGeom.mjGEOM_CAPSULE, w, np.array(a, float), np.array(b, float)); scn.ngeom += 1
    edge, ring = EDGES[TASKS[i]]
    for a, b in edge: cap(a, b, (0.20, 0.55, 0.95, 1))
    for a, b in ring: cap(a, b, (0.90, 0.25, 0.20, 1))
    for t_, p in TRAILS[i]:
        if t_ > ts or scn.ngeom >= scn.maxgeom: break
        mujoco.mjv_initGeom(scn.geoms[scn.ngeom], mujoco.mjtGeom.mjGEOM_SPHERE, np.array([0.0045, 0, 0]), p, np.eye(3).ravel(), np.array([0.95, 0.55, 0.05, 1.0])); scn.ngeom += 1
