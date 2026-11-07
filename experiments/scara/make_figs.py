"""Figures of Example 2 in the paper, drawn in the workcell: fig_manip_vanilla, fig_manip_uncon, fig_manip_rampc (Fig. 7) and fig_manip_compare (Fig. 8). Output: figs/."""
import os, pickle, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as cf
from config import C, fk, ik, phi_np, ra, TARGET_POS, CertArm, SPEC

M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0, P0, logs, van = M["X0"], M["Psi0"], M["logs"], M["van"]; n = X0.shape[1]; Y0 = fk(X0[:2])
SW = pickle.load(open(f"{cf.DATA}/sweep.pkl", "rb")) if os.path.exists(f"{cf.DATA}/sweep.pkl") else None
ura = np.array([u["status"] == "reached" for u in M["ura"]]); vr = np.array([v["status"] == "reached" for v in van]); rr = np.array([L["status"] == "reached" for L in logs])
ura_viol = np.array([(np.asarray(u["upeak"]) > np.asarray(cf.U_MAX)).any() for u in M["ura"]])
COL = dict(van="#e08214", ura="#2166ac", ra="#8856a7"); save = lambda fig, name: fig.savefig(os.path.join(cf.FIGS, f"{name}.png"), dpi=200) or print("saved", name)
_V = SPEC.get("view", [[-5.3, 5.3], [-0.6, 5.9]]); XL = tuple(_V[0]); YL = (max(_V[1][0], 0.6), _V[1][1] - (0 if SPEC.get('layout') == 'square' else 0)); PVIEW = (tuple(_V[0]), (-0.6, 5.9)) if SPEC.get('layout') != 'square' else ((-3.6, 3.6), (-0.6, 6.6)); g1 = np.linspace(*XL, 640); g2 = np.linspace(*YL, 330); G1, G2 = np.meshgrid(g1, g2); PS = np.broadcast_to(C.psi(G1, G2), G1.shape); PH = phi_np(G1, G2)

# shared start: farthest (in the output space) among the states where all three methods succeed
all3 = ura & vr & rr; d2 = ((Y0 - TARGET_POS[:, None])**2).sum(0)
pref = all3 if all3.any() else (ura & rr)                      # without a state where all three succeed, URA and RA-MPC share the start
if (pref & ura_viol).any(): pref = pref & ura_viol
i0 = int(np.where(pref)[0][np.argmax(d2[pref])]); x0s = X0[:, i0]
U = cf.ura_run(x0s, dense=True); Rr = logs[i0]; V = van[i0] if vr[i0] else van[int(np.where(vr)[0][np.argmax(d2[vr])])]
print(f"shared x0 = state {i0}, y0 {np.round(Y0[:, i0], 3)}, q0 {np.round(x0s[:2], 3)}; all-three successes {all3.sum()}; vanilla panel start {'shared' if vr[i0] else 'its own farthest reached state'}; times: vanilla {V['t'][-1]:.2f} s, URA {U['t']:.2f} s (peaks {np.round(U['upeak'], 1)}), RA-MPC {Rr['t'][-1]:.2f} s")

# ── paper figures ──
import sympy as _sp; _c = _sp.lambdify(_sp.Symbol('y1'), _sp.sympify(SPEC['path'])); XS_LABEL = (-0.62 * SPEC['X0'], float(_c(-0.55 * SPEC['X0'])) - 1.25); px = 1 / plt.rcParams["figure.dpi"]; FS_P = 30; C_ARM = (146 / 255, 149 / 255, 145 / 255); C_START = (1.0, 0.549, 0.0); C_END = (0.102, 0.596, 0.314); C_FIELD, FIELD_LW = "0.55", 1.2


def workcell(ax):
    """Dubins-style context: white background, the objects that are not being carried drawn as plain squares, the slots
    as dashed squares, the fixture in light grey; no bench colour, no racks."""
    X0w, YS, sig, cols = SPEC["X0"], SPEC["YS"], SPEC["SIGMA"], ("#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"); fx = SPEC["fixture"]; s = 0.4
    ax.add_patch(plt.Rectangle((-fx["w"] / 2, fx["y0"]), fx["w"], fx.get("top", 2.28) - fx["y0"], facecolor="#e6e6e6", edgecolor="#a6a6a6", lw=1.0, zorder=1.5))
    fwd = SPEC["source"][0] < 0; j_src = int(np.argmin([abs(v - SPEC["source"][1]) for v in YS])); j_dst = int(np.argmin([abs(v - SPEC["target"][1]) for v in YS]))
    objs = [(-X0w, YS[k], cols[k]) for k in range(len(YS)) if k != j_src] if fwd else [(X0w, YS[sig[k]], cols[k]) for k in range(len(YS)) if sig[k] != j_src]
    slots = [(X0w if fwd else -X0w, YS[j]) for j in range(len(YS)) if j != j_dst]
    for xx, yy, cc in objs: ax.add_patch(plt.Rectangle((xx - s / 2, yy - s / 2), s, s, facecolor=cc, edgecolor="k", lw=0.8, zorder=7))
    for xx, yy in slots: ax.add_patch(plt.Rectangle((xx - s / 2, yy - s / 2), s, s, facecolor="none", edgecolor="0.35", lw=1.4, ls=(0, (3, 2)), zorder=7))
    # the task's own object at its pick position (inside X_S) and its slot (inside X_T)
    k_task = j_src if fwd else int(np.where(np.array(sig) == j_src)[0][0]); ct = cols[k_task]; sx, sy = SPEC["source"]; tx, ty = SPEC["target"]
    ax.add_patch(plt.Rectangle((sx - s / 2, sy - s / 2), s, s, facecolor=ct, edgecolor="k", lw=0.8, zorder=7.5))
    ax.add_patch(plt.Rectangle((tx - s / 2, ty - s / 2), s, s, facecolor="none", edgecolor=ct, lw=1.8, ls=(0, (3, 2)), zorder=7.5))


def draw_robot(ax, x):
    q1, q2 = float(x[0]), float(x[1]); jx, jy = 4 * np.cos(q1), 4 * np.sin(q1); ex, ey = jx + 4 * np.cos(q1 + q2), jy + 4 * np.sin(q1 + q2)
    ax.plot([0, jx], [0, jy], "o-", alpha=0.4, lw=5, ms=7, color=C_ARM, zorder=6); ax.plot([jx, ex], [jy, ey], "o-", alpha=0.4, lw=5, ms=7, color=C_ARM, zorder=6); ax.plot(0, 0, "ks", markersize=8, zorder=6.5)


FV1 = np.linspace(*PVIEW[0], 30); FV2 = np.linspace(*PVIEW[1], 30); _A, _B = np.meshgrid(FV1, FV2); _k = C.k1(_A, _B)
FU = np.broadcast_to(np.asarray(_k[0], float), _A.shape).astype(float); FW = np.broadcast_to(np.asarray(_k[1], float), _A.shape).astype(float)


def paper_figs(V, U, Rr, x0s):
    PAN = {"fig_manip_vanilla": (vr, fk(V["x"].T[:2]).T, False, V["x"][0]), "fig_manip_uncon": (ura, U["y"], True, x0s), "fig_manip_rampc": (rr, fk(Rr["x"].T[:2]).T, True, x0s)}
    for name, (mask, yt, field, xr) in PAN.items():
        fig, ax = plt.subplots(figsize=(650 * px, 600 * px), layout="constrained"); fig.set_dpi(200)
        # field k1: opaque mid grey, drawn above the safe-set fill and below its boundary so that it survives printing
        if field: ax.streamplot(FV1, FV2, FU, FW, color=C_FIELD, linewidth=FIELD_LW, density=0.75, arrowsize=1.5, zorder=2.5)
        ax.contourf(G1, G2, PS, levels=[0, np.inf], colors=["#e6f2ff"], alpha=0.3, zorder=2); ax.contour(G1, G2, PS, levels=[0], colors=["#2166ac"], linewidths=2, zorder=3)
        ax.contourf(G1, G2, PH, levels=[-np.inf, 0], colors=["#fee0d2"], alpha=1.0, zorder=4); ax.contour(G1, G2, PH, levels=[0], colors=["#d6604d"], linewidths=2, zorder=5)
        workcell(ax); draw_robot(ax, xr)
        ax.scatter(*Y0[:, ~mask], s=50, color="#d73027", alpha=0.85, marker="x", linewidths=3, zorder=8); ax.scatter(*Y0[:, mask], s=50, color="#2166ac", alpha=0.35, marker="o", zorder=9)
        ax.plot(yt[:, 0], yt[:, 1], color="black", lw=3, alpha=0.65, linestyle="--", zorder=10); ax.scatter([yt[0, 0]], [yt[0, 1]], s=70, c=[C_START], marker="o", linewidths=0.4, zorder=20); ax.scatter([yt[-1, 0]], [yt[-1, 1]], s=70, c=[C_END], marker="x", linewidths=3, zorder=20)
        ax.text(XS_LABEL[0], XS_LABEL[1], r"$\mathcal{X}_S$", zorder=100, fontsize=FS_P); ax.text(TARGET_POS[0] + 0.33, TARGET_POS[1] + 0.38, r"$\mathcal{X}_T$", zorder=100, fontsize=FS_P)
        ax.set_xlabel("$y_1$", fontsize=FS_P - 4); ax.set_ylabel("$y_2$", fontsize=FS_P - 4); ax.xaxis.set_tick_params(labelsize=22); ax.yaxis.set_tick_params(labelsize=22); ax.set_xlim(*PVIEW[0]); ax.set_ylim(*PVIEW[1]); ax.set_aspect("equal")
        save(fig, name); plt.close(fig)
    FS_C = 24; fig, (axR, axW, axA) = plt.subplots(3, 1, figsize=(650 * px, 600 * px), layout="constrained", gridspec_kw={"height_ratios": [0.7, 1.0, 1.0]}); fig.set_dpi(200)
    if SW is not None:
        rows = SW["rows"]; Tsw = np.array([r["T"] for r in rows]); rho = np.array([100 * r["ra"].mean() for r in rows]); axR.plot(Tsw, rho, color=COL["ra"], marker="s", lw=2.5, ms=9, alpha=0.9)
        axR.set_ylim(5 * np.floor((rho.min() - 3) / 5), 5 * np.ceil((rho.max() + 3) / 5)); axR.set_xlim(Tsw.min() * 0.8, Tsw.max() * 1.25)
    axR.axvline(cf.T, color="0.2", lw=1.6, ls="--", alpha=0.9, zorder=2)
    if SW is not None:  # mark the success rate at the sampling period of the main run
        iT = int(np.argmin(np.abs(Tsw - cf.T))); ylo, yhi = axR.get_ylim()
        axR.plot(Tsw[iT], rho[iT], marker="s", ms=14, mfc="none", mec=COL["ra"], mew=2.5, zorder=5)
        axR.annotate(rf"$\rho = {rho[iT]:.1f}\,\%$", xy=(Tsw[iT], rho[iT]), xytext=(cf.T * 1.7, ylo + 0.2 * (yhi - ylo)), fontsize=15, color=COL["ra"], ha="left", va="center",
                     arrowprops=dict(arrowstyle="-", color="0.3", lw=1.0, shrinkA=0, shrinkB=7), bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9), zorder=6)
    axR.set_xscale("log"); axR.set_xlabel("$T$ [s]", fontsize=FS_C - 6); axR.set_ylabel(r"$\rho$ [%]", fontsize=FS_C - 2); axR.xaxis.set_tick_params(labelsize=14); axR.yaxis.set_tick_params(labelsize=14); axR.grid(True, alpha=0.3, which="both")
    tu, uu = U["tt"], U["u"]; tr_, ur = Rr["t"], Rr["u"]
    for ax, ch, ub, lbl in zip((axW, axA), (0, 1), cf.U_MAX, (r"$\tau_1$", r"$\tau_2$")):
        ax.axhspan(-ub, ub, color="#f0f0f0", alpha=0.6, zorder=0); ax.axhline(ub, color="black", lw=1.5, ls="--", alpha=0.85, zorder=5); ax.axhline(-ub, color="black", lw=1.5, ls="--", alpha=0.85, zorder=5); ax.axhline(0, color="gray", lw=0.8, ls=":", zorder=4)
        LOGT = tu[-1] > 6 * tr_[-1]; ax.plot(np.maximum(tu, 1e-3) if LOGT else tu, uu[:, ch], color=COL["ura"], lw=3, alpha=0.9, zorder=10)
        ts_ = np.append(tr_[:-1], tr_[-1]); ax.step(np.maximum(ts_, 1e-3) if LOGT else ts_, np.append(ur[:, ch], ur[-1, ch]), where="post", color=COL["ra"], lw=3, alpha=0.9, zorder=10)
        if LOGT: ax.set_xscale("log"); ax.set_xlim(cf.T / 2, tu[-1] * 1.15)
        else: ax.set_xlim(0, tu[-1] * 1.02)
        ym = max(ub, np.abs(uu[:, ch]).max()); ax.set_ylabel(lbl, fontsize=FS_C); ax.set_ylim(-1.25 * ym, 1.25 * ym); ax.xaxis.set_tick_params(labelsize=14); ax.yaxis.set_tick_params(labelsize=14); ax.grid(True, alpha=0.25, which="both")
    axA.set_xlabel("Time [s]", fontsize=FS_C - 6); save(fig, "fig_manip_compare"); plt.close(fig)


paper_figs(V, U, Rr, x0s)
pickle.dump(dict(i0=i0, x0=x0s, t_van=float(V["t"][-1]), t_ura=float(U["t"]), t_ra=float(Rr["t"][-1]), all3=int(all3.sum())), open(f"{cf.DATA}/shared_x0.pkl", "wb"))
