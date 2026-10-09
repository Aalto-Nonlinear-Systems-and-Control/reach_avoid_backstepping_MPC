"""Top-view figures of the Panda kitting task in the Dubins style (white background, sets and points first, workcell and
robot light grey): fig_franka_vanilla/uncon/rampc (+ _far), fig_franka_compare (torques against the bounds), and a
four-panel outcome check."""
import os, pickle, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config_franka as cf
from config_franka import C, SPEC, ALPHA, Z0, TARGET
from cert_franka import phi_xy, RT
from franka_model import JOINTS_f, NQ

M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0, logs, van = M["X0"], M["logs"], M["van"]; n = X0.shape[1]; Y0 = C.y(X0)
ura = np.array([u["status"] == "reached" for u in M["ura"]]); fv = np.array([len(v["J"]) > 0 for v in van]); vr = np.array([v["status"] == "reached" for v in van]); rr = np.array([L["status"] == "reached" for L in logs])
ura_viol = np.array([(np.asarray(u["upeak"]) > np.asarray(cf.U_MAX)).any() for u in M["ura"]])
COL = dict(van="#e08214", ura="#2166ac", ra="#8856a7"); C_ARM = (146 / 255, 149 / 255, 145 / 255); C_START = (1.0, 0.549, 0.0); C_END = (0.102, 0.596, 0.314); C_GRAY = (178 / 255, 187 / 255, 205 / 255, 0.45)
save = lambda fig, name, pdf=False: [fig.savefig(os.path.join(cf.FIGS, f"{name}.{e}"), dpi=200) for e in (("png", "pdf") if pdf else ("png",))] and print("saved", name)
VX, VY = (-3.6 * ALPHA, 3.6 * ALPHA), (-0.6 * ALPHA, 6.6 * ALPHA)
g1 = np.linspace(*VX, 500); g2 = np.linspace(*VY, 500); G1, G2 = np.meshgrid(g1, g2); GY = np.vstack([G1.ravel(), G2.ravel(), Z0 * np.ones(G1.size), np.zeros(G1.size)])
PS = C.psi_of_y(GY).reshape(G1.shape); PH = phi_xy(GY).reshape(G1.shape)
fv1 = np.linspace(*VX, 30); fv2 = np.linspace(*VY, 30); F1, F2 = np.meshgrid(fv1, fv2); FK = C.k1_of_y(np.vstack([F1.ravel(), F2.ravel(), Z0 * np.ones(F1.size), np.zeros(F1.size)])); FU, FW = FK[0].reshape(F1.shape), FK[1].reshape(F1.shape)


def workcell(ax):
    X0w, YS, sig, cols, fx, s = ALPHA * SPEC["X0"], [ALPHA * v for v in SPEC["YS"]], SPEC["SIGMA"], ("#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"), SPEC["fixture"], 0.04
    ax.add_patch(plt.Rectangle((-ALPHA * fx["w"] / 2, ALPHA * fx["y0"]), ALPHA * fx["w"], ALPHA * (fx["top"] - fx["y0"]), facecolor="#e6e6e6", edgecolor="#a6a6a6", lw=1.0, zorder=1.5))
    ax.text(0, ALPHA * (fx["y0"] + fx["top"]) / 2, "obstacle", rotation=90, ha="center", va="center", fontsize=15, color="0.35", zorder=7)
    j_src = int(np.argmin([abs(v - SPEC["source"][1]) for v in SPEC["YS"]])); j_dst = int(np.argmin([abs(v - SPEC["target"][1]) for v in SPEC["YS"]])); fwd = SPEC["source"][0] < 0
    objs = [(-X0w, YS[k], cols[k]) for k in range(4) if k != j_src] if fwd else [(X0w, YS[sig[k]], cols[k]) for k in range(4) if sig[k] != j_src]
    for xx, yy, cc in objs: ax.add_patch(plt.Rectangle((xx - s / 2, yy - s / 2), s, s, facecolor=cc, edgecolor="k", lw=0.8, zorder=7))
    for j in range(4):
        if j != j_dst: xx = X0w if fwd else -X0w; ax.add_patch(plt.Rectangle((xx - s / 2, YS[j] - s / 2), s, s, facecolor="none", edgecolor="0.35", lw=1.4, ls=(0, (3, 2)), zorder=7))
    k_task = j_src if fwd else int(np.where(np.array(sig) == j_src)[0][0]); ct = cols[k_task]; sx, sy = ALPHA * np.array(SPEC["source"]); tx, ty = TARGET
    ax.add_patch(plt.Rectangle((sx - s / 2, sy - s / 2), s, s, facecolor=ct, edgecolor="k", lw=0.8, zorder=7.5)); ax.add_patch(plt.Rectangle((tx - s / 2, ty - s / 2), s, s, facecolor="none", edgecolor=ct, lw=1.8, ls=(0, (3, 2)), zorder=7.5))


def robot(ax, x):
    """Top-view footprint: with q3 = q5 = 0 the arm lies in one vertical plane, so all links project onto the ray from the
    base to the gripper; drawn as a light bar from the base to the tool centre point with a wrist marker."""
    P = np.array(JOINTS_f(np.asarray(x)[:NQ])); tcp = np.array([-P[1, 4], P[0, 4]]); wr = np.array([-P[1, 3], P[0, 3]])
    ax.plot([0, tcp[0]], [0, tcp[1]], "-", alpha=0.4, lw=6, color=C_ARM, solid_capstyle="round", zorder=6); ax.plot(*wr, "o", ms=8, color=C_ARM, alpha=0.7, zorder=6.2)
    ax.plot(0, 0, "ks", ms=8, zorder=6.5)


def panel(name, mask, yt, field, xr):
    px = 1 / plt.rcParams["figure.dpi"]; FS = 30; fig, ax = plt.subplots(figsize=(650 * px, 600 * px), layout="constrained"); fig.set_dpi(200)
    if field: ax.streamplot(fv1, fv2, FU, FW, color=C_GRAY, linewidth=0.9, density=0.9, arrowsize=0.9, zorder=1)
    ax.contourf(G1, G2, PS, levels=[0, np.inf], colors=["#e6f2ff"], alpha=0.3, zorder=2); ax.contour(G1, G2, PS, levels=[0], colors=["#2166ac"], linewidths=2, zorder=3)
    ax.contourf(G1, G2, PH, levels=[-np.inf, 0], colors=["#fee0d2"], zorder=4); ax.contour(G1, G2, PH, levels=[0], colors=["#d6604d"], linewidths=2, zorder=5)
    workcell(ax); robot(ax, xr)
    ax.scatter(*Y0[:2, ~mask], s=50, color="#d73027", alpha=0.85, marker="x", linewidths=3, zorder=8); ax.scatter(*Y0[:2, mask], s=50, color="#2166ac", alpha=0.35, marker="o", zorder=9)
    ax.plot(yt[:, 0], yt[:, 1], color="black", lw=3, alpha=0.65, ls="--", zorder=10); ax.scatter([yt[0, 0]], [yt[0, 1]], s=70, c=[C_START], zorder=20); ax.scatter([yt[-1, 0]], [yt[-1, 1]], s=70, c=[C_END], marker="x", linewidths=3, zorder=20)
    ax.text(-0.155, 0.175, r"$\mathcal{X}_S$", fontsize=FS, zorder=100); ax.text(TARGET[0] + 0.033, TARGET[1] + 0.038, r"$\mathcal{X}_T$", fontsize=FS, zorder=100)
    ax.set_xlabel("$y_1$ [m]", fontsize=FS - 6); ax.set_ylabel("$y_2$ [m]", fontsize=FS - 6); ax.tick_params(labelsize=18); ax.set_xlim(*VX); ax.set_ylim(*VY); ax.set_aspect("equal"); save(fig, name, pdf=True); plt.close(fig)


def ytraj(L): return C.y(np.asarray(L["x"]).T).T


d2 = ((Y0[:2] - TARGET[:, None])**2).sum(0); both = ura & rr
pref = both & vr if (both & vr).any() else both; pref = pref & ura_viol if (pref & ura_viol).any() else pref
for suffix, pick in (("", lambda m: int(np.where(m)[0][np.argmin(np.abs(np.sqrt(d2[m]) - np.median(np.sqrt(d2[m]))))])), ("_far", lambda m: int(np.where(m)[0][np.argmax(d2[m])]))):
    i0 = pick(pref) if suffix == "" else pick(both); iv = i0 if vr[i0] else pick(vr)
    U = cf.ura_run(X0[:, i0], dense=True, every=10); R = logs[i0]; V = van[iv]
    panel("fig_franka_vanilla" + suffix, vr, ytraj(V), False, V["x"][0]); panel("fig_franka_uncon" + suffix, ura, U["y"], True, X0[:, i0]); panel("fig_franka_rampc" + suffix, rr, ytraj(R), True, X0[:, i0])
    print(f"{suffix or 'shared'}: start {i0} y0 {np.round(Y0[:2, i0], 3)}; URA {U['t']:.2f} s peaks {np.round(U['upeak'], 1)}; RA-MPC {R['t'][-1]:.2f} s; vanilla from {iv} {V['t'][-1]:.2f} s")
    px = 1 / plt.rcParams["figure.dpi"]; fig, axs = plt.subplots(2, 2, figsize=(1000 * px, 650 * px), layout="constrained", sharex=True); fig.set_dpi(200)
    for j, (a, ub, lbl) in enumerate(zip(axs.ravel(), cf.U_MAX, (r"$\tau_1$", r"$\tau_2$", r"$\tau_4$", r"$\tau_6$"))):
        a.axhspan(-ub, ub, color="#f0f0f0", alpha=0.6, zorder=0); a.axhline(ub, color="k", lw=1.5, ls="--"); a.axhline(-ub, color="k", lw=1.5, ls="--"); a.axhline(0, color="gray", lw=0.8, ls=":")
        a.plot(U["tt"], U["u"][:, j], color=COL["ura"], lw=2.5, label="URA"); ts = np.append(R["t"][:-1], R["t"][-1]); a.step(ts, np.append(R["u"][:, j], R["u"][-1, j]), where="post", color=COL["ra"], lw=2.5, label="RA-MPC")
        ym = max(ub, np.abs(U["u"][:, j]).max()); a.set_ylim(-1.25 * ym, 1.25 * ym); a.set_ylabel(lbl + " [N m]", fontsize=16); a.grid(alpha=0.25)
    axs[1, 0].set_xlabel("Time [s]", fontsize=16); axs[1, 1].set_xlabel("Time [s]", fontsize=16); axs[0, 0].legend(fontsize=11); save(fig, "fig_franka_compare" + suffix, pdf=True); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(18, 5.6), constrained_layout=True)
for a, ok, ttl, col in zip(ax, (vr, ura, rr), ("vanilla MPC", f"URA (s = {cf.URA_SCALE:.2f})", "RA-MPC"), (COL["van"], COL["ura"], COL["ra"])):
    a.contourf(G1, G2, PS, levels=[0, np.inf], colors=["#e6eef7"]); a.contour(G1, G2, PS, levels=[0], colors="k", linewidths=1.0); a.contourf(G1, G2, -PH, levels=[0, 1e9], colors=["#bfe5bf"])
    a.scatter(*Y0[:2, ok], s=14, c=col); a.scatter(*Y0[:2, ~ok], s=18, c="0.45", marker="x"); a.set_aspect("equal"); a.set_xlim(-0.35, 0.36); a.set_ylim(0.15, 0.52); a.set_title(f"{ttl}: {ok.sum()}/{n} reach the slot", fontsize=11)
    if ttl == "vanilla MPC": later = fv & ~vr; a.scatter(*Y0[:2, later], s=34, facecolors="none", edgecolors=COL["van"], linewidths=1.2, zorder=5); a.set_title(f"vanilla MPC: {ok.sum()}/{n} reach the slot (feasible at k = 0: {fv.sum()}, open circles lose feasibility)", fontsize=10)
fig.suptitle(f"Panda (simplified dynamics, MuJoCo), task {cf.TASK}: s = {cf.SCALE:.2f}, T = {1e3*cf.T:g} ms, N = {cf.N}, bounds {cf.U_MAX}", fontsize=11); save(fig, "verify_outcomes"); plt.close(fig)
