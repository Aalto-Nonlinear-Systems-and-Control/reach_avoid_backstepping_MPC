"""Full kitting cycle on the SCARA: the four cubes go from the left rack to their slots (forward tasks, crossing assignment
SIGMA) and back (reverse tasks), eight tasks in all. Each task: pick (Z axis from the transit height down to the cube,
grasp, up to the transport height; the arm holds its position), transport by RA-MPC with the task's own certificate and
scale under the common torque bounds (certified; ra_run_task.py, from the actual state), place (centre over the slot,
Z axis down, release, up to the transit height), transit to the next cube (point-to-point joint-space minimum-jerk
move at 0.40 m, above the obstacle, shortest duration whose feedforward torque stays within 90 % of the bounds; not
certified). Pick, place and transit use clipped computed torque with the same bounds. The planar part
runs on the paper's plant; each held-input period is re-simulated on model/scara_arm.xml (torques times C_TAU) from the
plant state at its start, which checks the model step by step and gives the 0.5 ms states of the video (an open-loop
replay of the whole 70 s would amplify round-off by about 40 per task, as the recorded inputs carry no feedback).
Writes data/cycle.npz and data/cycle_ra/<task>.pkl."""
import os, sys, json, time, pickle, subprocess, numpy as np, sympy as sp, mujoco
from scara_common import *
from cert_top import fk, ik, jac, jacdot_qdot, mass, coriolis_qdot
from ra_mpc_top import rk4_np
AT = os.path.join(HERE, "..", "arm_top"); PY = sys.executable
TASKS = [(n, d) for d in ("fwd", "rev") for n in ("red", "blue", "green", "orange")]
SPECS = {f"{n}_{d}_sq": json.load(open(os.path.join(AT, "tasks", f"{n}_{d}_sq.json"))) for n, d in TASKS}
T = 0.01; UMAX = np.array(json.load(open(os.path.join(AT, "data", "scale_red_fwd_sq.json")))["U_MAX"]); NSUB = int(round(T / 0.0005)); KP, KD = 400.0, 40.0
D_DOWN, D_GRASP, D_UP, Z_TRANSIT = 1.0, 0.3, 1.0, 0.40; NAMES = ("red", "blue", "green", "orange"); os.makedirs(os.path.join(HERE, "data", "cycle_ra"), exist_ok=True)
X, U, PH, TK = [], [], [], []; ZSEG = []; EV = []; INFO = []


def now(): return len(U) * T


def step(x, u, ph, tk):
    for _ in range(NSUB): x = rk4_np(x, u, T / NSUB)
    X.append(x.copy()); U.append(np.array(u, float)); PH.append(ph); TK.append(tk); return x


def ct(x, yr, vr, ar):
    q, dq = x[:2], x[2:]; J = jac(q); yd = J @ dq; qdd = np.linalg.solve(J, ar + KP * (yr - fk(q)) + KD * (vr - yd) - jacdot_qdot(q, dq))
    u = mass(q) @ qdd + coriolis_qdot(q, dq); return np.clip(u, -UMAX, UMAX), bool(np.any(np.abs(u) > UMAX))


def move(x, y_goal, dur, ph, tk, settle=1.0, tol=1e-3):
    """Planar minimum-jerk move of the tool point (task units), then regulation until settled (or settle seconds)."""
    y0 = fk(x[:2]); D = y_goal - y0; k = 0; sat = 0
    while True:
        t = k * T; s = min(t / dur, 1.0); on = t < dur
        yr = y0 + D * minjerk(s); vr = D * (30 * s**2 - 60 * s**3 + 30 * s**4) / dur * on; ar = D * (60 * s - 180 * s**2 + 120 * s**3) / dur**2 * on
        if t >= dur and ((np.linalg.norm(fk(x[:2]) - y_goal) < tol and np.linalg.norm(jac(x[:2]) @ x[2:]) < 2 * tol) or t >= dur + settle): break
        u, c = ct(x, yr, vr, ar); sat += c; x = step(x, u, ph, tk); k += 1
    return x, k * T, sat


def move_joint(x, q_goal, ph, tk, settle=1.0, tol=1e-4, KPJ=400.0, KDJ=40.0):
    """Point-to-point joint-space minimum-jerk move (SCARA transit), duration from the feedforward torque."""
    q0 = x[:2].copy(); D = q_goal - q0; S = np.linspace(0, 1, 201)
    ref = lambda s, dur: (q0 + D * minjerk(s), D * (30 * s**2 - 60 * s**3 + 30 * s**4) / dur, D * (60 * s - 180 * s**2 + 120 * s**3) / dur**2)
    for dur in np.arange(1.0, 8.01, 0.25):
        if all(np.all(np.abs(mass(qr) @ ar + coriolis_qdot(qr, vr)) <= 0.9 * UMAX) for qr, vr, ar in (ref(s, dur) for s in S)): break
    k = 0; sat = 0
    while True:
        t = k * T; s = min(t / dur, 1.0); qr, vr, ar = ref(s, dur)
        if t >= dur: vr, ar = 0 * vr, 0 * ar
        q, dq = x[:2], x[2:]
        if t >= dur and ((np.abs(q - q_goal).max() < tol and np.abs(dq).max() < 10 * tol) or t >= dur + settle): break
        u = mass(q) @ (ar + KPJ * (qr - q) + KDJ * (vr - dq)) + coriolis_qdot(q, dq); sat += bool(np.any(np.abs(u) > UMAX)); x = step(x, np.clip(u, -UMAX, UMAX), ph, tk); k += 1
    return x, k * T, sat


def hold(x, y_goal, dur, ph, tk):
    for _ in range(int(round(dur / T))): x = step(x, ct(x, y_goal, 0 * y_goal, 0 * y_goal)[0], ph, tk)
    return x


src0 = np.array(SPECS["red_fwd_sq"]["source"]); x = np.r_[ik(src0.reshape(2, 1))[:, 0], 0.0, 0.0]; X.append(x.copy()); t_start = time.time()
for i, (n, d) in enumerate(TASKS):
    task = f"{n}_{d}_sq"; S = SPECS[task]; src, dst = np.array(S["source"]), np.array(S["target"]); k_obj = NAMES.index(n)
    t0 = now(); ZSEG += [(t0, t0 + D_DOWN, Z_TRANSIT, Z_LOW), (t0 + D_DOWN, t0 + D_DOWN + D_GRASP, Z_LOW, Z_LOW), (t0 + D_DOWN + D_GRASP, t0 + D_DOWN + D_GRASP + D_UP, Z_LOW, Z0)]
    x = hold(x, src, D_DOWN + D_GRASP + D_UP, 0, i); t_grasp = t0 + D_DOWN; t_ra = now()
    np.save(os.path.join(HERE, "data", "cycle_ra", f"{task}_x0.npy"), x); out = os.path.join(HERE, "data", "cycle_ra", f"{task}.pkl")
    r = subprocess.run([PY, os.path.join(HERE, "ra_run_task.py"), task, os.path.join(HERE, "data", "cycle_ra", f"{task}_x0.npy"), out], capture_output=True, text=True)
    print("  " + (r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-800:]), flush=True)
    L = pickle.load(open(out, "rb")); assert L["status"] == "reached", f"RA-MPC failed for {task}"
    for u in L["u"]: x = step(x, u, 1, i)
    mism = float(np.abs(x - L["x"][-1]).max()); t_entry = now(); ZSEG.append((t_ra, t_entry, Z0, Z0))
    x, dc, sat_c = move(x, dst, 0.6, 2, i); t_low = now(); ZSEG.append((t_entry, t_low, Z0, Z0))
    ZSEG += [(t_low, t_low + D_DOWN, Z0, Z_LOW), (t_low + D_DOWN, t_low + D_DOWN + D_GRASP, Z_LOW, Z_LOW), (t_low + D_DOWN + D_GRASP, t_low + D_DOWN + D_GRASP + D_UP, Z_LOW, Z_TRANSIT)]
    x = hold(x, dst, D_DOWN + D_GRASP + D_UP, 2, i); t_rel = t_low + D_DOWN
    EV.append((k_obj, t_grasp + T_FINGER, t_rel, t_grasp)); sat_t = 0; dt_ = 0.0
    if i + 1 < len(TASKS):
        nxt = np.array(SPECS[f"{TASKS[i + 1][0]}_{TASKS[i + 1][1]}_sq"]["source"]); t_tr = now()
        x, dt_, sat_t = move_joint(x, ik(nxt.reshape(2, 1))[:, 0], 3, i); ZSEG.append((t_tr, now(), Z_TRANSIT, Z_TRANSIT))
        assert np.linalg.norm(fk(x[:2]) - nxt) < 1e-3, f"transit to task {i + 2} did not arrive"
    INFO.append(dict(task=task, obj=k_obj, t_pick=t0, t_grasp=t_grasp, t_ra=t_ra, t_entry=t_entry, t_lower=t_low, t_release=t_rel, ra_time=t_entry - t_ra, ra_steps=len(L["u"]),
                     ra_solve_ms=float(np.median(L["solve_ms"])), scale=float(L["scale"]), replay_mismatch=mism, centre_s=dc, centre_clipped=int(sat_c), transit_s=dt_, transit_clipped=int(sat_t)))
    print(f"{i + 1}/8 {task}: RA-MPC {t_entry - t_ra:.2f} s (scale {L['scale']:.3f}, replay mismatch {mism:.0e}), centring {dc:.2f} s, transit {dt_:.2f} s ({sat_t} clipped) | t = {now():.2f} s ({time.time() - t_start:.0f}s)", flush=True)
t_end = now()
ma = mujoco.MjModel.from_xml_path(f"{HERE}/model/scara_arm.xml"); da = mujoco.MjData(ma); Q, err = [X[0][:2].copy()], 0.0
for k, u in enumerate(U):
    da.qpos[:], da.qvel[:] = X[k][:2], X[k][2:]; da.ctrl[:] = C_TAU * u
    for _ in range(NSUB): mujoco.mj_step(ma, da); Q.append(da.qpos.copy())
    err = max(err, np.abs(np.r_[da.qpos, da.qvel] - X[k + 1]).max())
Q = np.array(Q); TS = np.arange(len(Q)) * ma.opt.timestep; Z = np.full(len(TS), Z_TRANSIT)
for a, b, za, zb in ZSEG:
    sel = (TS >= a - 1e-12) & (TS <= b + 1e-12); Z[sel] = za + (zb - za) * minjerk((TS[sel] - a) / max(b - a, 1e-9))
GAP = np.full(len(TS), GAP_OPEN)
for k_obj, t_att, t_rel, t_gr in EV:
    sel = (TS >= t_gr) & (TS < t_rel); GAP[sel] = GAP_OPEN + (GAP_CLOSED - GAP_OPEN) * minjerk((TS[sel] - t_gr) / T_FINGER)
    sel = (TS >= t_rel) & (TS < t_rel + T_FINGER); GAP[sel] = GAP_CLOSED + (GAP_OPEN - GAP_CLOSED) * minjerk((TS[sel] - t_rel) / T_FINGER)
PHQ = np.r_[PH[0], np.repeat(PH, NSUB)]; TKQ = np.r_[TK[0], np.repeat(TK, NSUB)]
y1s, y2s = sp.symbols("y1 y2"); PSI = {t: sp.lambdify((y1s, y2s), sp.sympify(S["psi"]), "numpy") for t, S in SPECS.items()}
for i, inf in enumerate(INFO):
    sel = (PHQ == 1) & (TKQ == i); Yt = fk(Q[sel].T); inf["min_psi_raw"] = float(np.min(PSI[inf["task"]](*Yt)) / SPECS[inf["task"]]["pmax"])
print(f"MuJoCo, period by period (torques x {C_TAU:g}): max one-period state difference to the paper's plant {err:.1e}; total {t_end:.2f} s; peak |tau| world {np.round(C_TAU * np.abs(np.array(U)).max(0), 3)} within {C_TAU * UMAX}")
print("min psi/pmax along each transport: " + ", ".join(f"{inf['task']} {inf['min_psi_raw']:.4f}" for inf in INFO))
np.savez(f"{HERE}/data/cycle.npz", t=TS, q=Q, z=Z, gap=GAP, phase=PHQ, task=TKQ, u=C_TAU * np.array(U), umax=C_TAU * UMAX, x_paper=np.array(X), events=np.array(EV), t_end=t_end,
         tasks=np.array([inf["task"] for inf in INFO]), info=json.dumps(INFO))
json.dump(dict(tasks=INFO, mujoco_max_one_period_state_diff=err, total_s=t_end, umax_world=(C_TAU * UMAX).tolist()), open(f"{HERE}/data/cycle_summary.json", "w"), indent=1)
