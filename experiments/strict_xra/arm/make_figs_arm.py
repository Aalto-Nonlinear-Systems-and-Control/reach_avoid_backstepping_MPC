"""Figures for the arm: verification of Lemma 2 / Theorem 2 / Corollary 1, period sweep, paper figures
(fig_manip_vanilla/uncon/rampc, fig_manip_compare, fig_manip_cost) in the style of make_paper_figs_manipulator.ipynb."""
import os, pickle, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config_arm as cf
from config_arm import C, fk, ik, phi_np, ra, TARGET_POS, CertArm

M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0, P0, logs, van = M["X0"], M["Psi0"], M["logs"], M["van"]; n = X0.shape[1]; Y0 = fk(X0[:2])
L2 = pickle.load(open(f"{cf.DATA}/lemma2.pkl", "rb")) if os.path.exists(f"{cf.DATA}/lemma2.pkl") else None
SW = pickle.load(open(f"{cf.DATA}/sweep.pkl", "rb")) if os.path.exists(f"{cf.DATA}/sweep.pkl") else None
ura = np.array([u["status"] == "reached" for u in M["ura"]]); fv = np.array([len(v["J"]) > 0 for v in van]); vr = np.array([v["status"] == "reached" for v in van]); rr = np.array([L["status"] == "reached" for L in logs])
ura_viol = np.array([(np.asarray(u["upeak"]) > np.asarray(cf.U_MAX)).any() for u in M["ura"]])
t_ura = np.array([u["t"] for u in M["ura"]]); t_ra = np.array([L["check"]["t_entry"] if L["status"] == "reached" else np.nan for L in logs])
COL = dict(van="#e08214", ura="#2166ac", ra="#8856a7"); save = lambda fig, name, pdf=False: [fig.savefig(os.path.join(cf.FIGS, f"{name}.{e}"), dpi=200) for e in (("png", "pdf") if pdf else ("png",))] and print("saved", name)
g1 = np.linspace(-3.5, 8.2, 500); g2 = np.linspace(-3.0, 3.0, 300); G1, G2 = np.meshgrid(g1, g2); PS = np.broadcast_to(C.psi(G1, G2), G1.shape); PH = phi_np(G1, G2)
yv = np.linspace(-3.5, 8.2, 24); yw = np.linspace(-3.0, 3.0, 14); VY1, VY2 = np.meshgrid(yv, yw); kk = C.k1(VY1, VY2)
KU = np.broadcast_to(np.asarray(kk[0], float), VY1.shape).astype(float); KV = np.broadcast_to(np.asarray(kk[1], float), VY1.shape).astype(float)


def sets(ax, cert_region=True):
    if cert_region: ax.contourf(G1, G2, PS, levels=[0, cf.DELTA + C.shift, 10], colors=["#e6eef7", "#c9dcf0"], alpha=0.9)
    ax.contour(G1, G2, PS, levels=[0], colors="k", linewidths=1.0); ax.contourf(G1, G2, -PH, levels=[0, 1e9], colors=["#bfe5bf"]); ax.contour(G1, G2, PH, levels=[0], colors="g", linewidths=1.0)
    th = np.linspace(0, 2 * np.pi, 300); ax.plot(8 * np.cos(th), 8 * np.sin(th), "k:", lw=0.8); ax.plot(0, 0, "ko", ms=3)
    ax.set_aspect("equal"); ax.set_xlim(-3.5, 8.2); ax.set_ylim(-3.0, 3.0); ax.set_xlabel("$y_1$"); ax.set_ylabel("$y_2$")


# shared start: farthest (in the output space) among the states where all three methods succeed
all3 = ura & vr & rr; d2 = ((Y0 - TARGET_POS[:, None])**2).sum(0); pref = all3 & ura_viol if (all3 & ura_viol).any() else all3
i0 = int(np.where(pref)[0][np.argmax(d2[pref])]); x0s = X0[:, i0]
U = cf.ura_run(x0s, dense=True); V = van[i0]; Rr = logs[i0]
print(f"shared x0 = state {i0}, y0 {np.round(Y0[:, i0], 3)}, q0 {np.round(x0s[:2], 3)}; all-three successes {all3.sum()} (URA exceeding the bounds among them {(all3 & ura_viol).sum()}); times: vanilla {V['t'][-1]:.2f} s, URA {U['t']:.2f} s (peaks {np.round(U['upeak'], 1)}), RA-MPC {Rr['t'][-1]:.2f} s")

# ── V1 outcomes ──
fig, ax = plt.subplots(1, 4, figsize=(20, 4.2), constrained_layout=True)
for a, ok, ttl, col in zip(ax, (fv, ura, rr), ("vanilla MPC", rf"URA (scale $s = {cf.URA_SCALE:g}$)", "RA-MPC"), (COL["van"], COL["ura"], COL["ra"])):
    sets(a); a.scatter(*Y0[:, ok], s=14, c=col); a.scatter(*Y0[:, ~ok], s=18, c="0.45", marker="x", linewidths=0.9); a.set_title(f"{ttl}: {ok.sum()}/{n}", fontsize=11)
both = ura & rr; a = ax[3]; a.scatter(t_ura[both], t_ra[both], s=12, c=COL["ra"]); a.set_xscale("log"); a.set_yscale("log"); lim = [0.05, 100]; a.plot(lim, lim, "k--", lw=0.8); a.set_xlim(lim); a.set_ylim(lim); a.grid(alpha=0.3, which="both")
a.set_xlabel("time to target, URA (s)"); a.set_ylabel("time to target, RA-MPC (s)"); a.set_title(f"{both.sum()} common successes: median {np.median(t_ura[both]):.2f} s vs {np.median(t_ra[both]):.2f} s", fontsize=10)
fig.suptitle(rf"Arm, {n} initial states at rest, positions uniform in the band; terminal controller $s = {cf.SCALE:g}$, $T = {1e3*cf.T:g}$ ms, $N = {cf.N}$, bounds {cf.U_MAX}", fontsize=10)
save(fig, "verify_outcomes"); plt.close(fig)

# ── V2 Lemma 2 ──
if L2 is not None:
    TS = L2["Tstar2"]; fig, ax = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
    H = L2["traj"]["psi"]; tt = np.arange(H.shape[0]) * L2["traj"]["dt"]; a = ax[0]
    for j in range(H.shape[1]): a.plot(tt, H[:, j], lw=0.8, color=COL["ura"], alpha=0.6)
    a.axhline(cf.DELTA, color="k", ls="--", lw=1); a.set_xscale("log"); a.set_xlabel("time (s)"); a.set_ylabel(r"$\Psi(x^s(t_k))$ at the sampling instants"); a.grid(alpha=0.3)
    a.set_title(rf"(a) held torque, $T = {1e3*cf.T:g}$ ms, {H.shape[1]} states of $\mathcal{{X}}_{{RA}}^{{\delta}}$", fontsize=9.5)
    a = ax[1]; Tg = np.array([r["T"] for r in L2["onestep"]]); dm = np.array([r["dmin"] for r in L2["onestep"]]); m2 = np.array([r["M2"] for r in L2["onestep"]])
    a.plot(Tg, dm / Tg, "o-", color=COL["ura"], label="smallest one-step rate of $\\Psi$"); a.plot(Tg, L2["eta"] - 0.5 * m2 * Tg, "s--", color="0.35", ms=4, label=r"second-order bound $\eta - \frac{1}{2} M_2 T$")
    a.axhline(0, color="k", lw=0.8); a.axvline(TS, color="g", ls="--"); a.axvline(cf.T, color="r", ls=":"); a.set_xscale("log"); a.set_ylim(-2.5 * L2["eta"], 1.6 * L2["eta"])
    a.text(TS, 0.93, rf" $T^\star$ = {1e3*TS:.1f} ms", color="g", transform=a.get_xaxis_transform(), fontsize=9); a.text(cf.T, 0.93, rf"$T$ = {1e3*cf.T:g} ms ", color="r", ha="right", transform=a.get_xaxis_transform(), fontsize=9)
    a.set_xlabel("hold period $T$ (s)"); a.set_ylabel("rate of $\\Psi$ over one period"); a.set_title(rf"(b) one-step increase on $\mathcal{{X}}_{{RA}}^{{\delta}}$ ({L2['n_delta']} states)", fontsize=9.5); a.legend(fontsize=8, loc="center left", bbox_to_anchor=(0, 0.25)); a.grid(alpha=0.3, which="both")
    a = ax[2]; Tc = np.array([r["T"] for r in L2["closed_loop"]]); nn = L2["closed_loop"][0]["delta"]["n"]
    a.plot(Tc, [100 * r["delta"]["reached"] / nn for r in L2["closed_loop"]], "o-", color=COL["ura"], label=r"reach $\mathcal{X}_T$ safely, from $\mathcal{X}_{RA}^{\delta}$")
    a.plot(Tc, [100 * r["zero"]["reached"] / nn for r in L2["closed_loop"]], "s-", color=COL["ra"], label=r"reach $\mathcal{X}_T$ safely, from $\mathcal{X}_{RA}$")
    a.plot(Tc, [100 * r["delta"]["monotone"] / nn for r in L2["closed_loop"]], "^--", color="0.4", label=r"$\Psi$ increasing at every sample, from $\mathcal{X}_{RA}^{\delta}$")
    a.axvline(TS, color="g", ls="--"); a.axvline(cf.T, color="r", ls=":"); a.set_xscale("log"); a.set_ylim(-3, 104); a.set_xlabel("hold period $T$ (s)"); a.set_ylabel(f"share of {nn} initial states (%)"); a.legend(fontsize=8, loc="lower left"); a.grid(alpha=0.3, which="both")
    a.set_title("(c) sampled-data closed loop against the period", fontsize=9.5); save(fig, "verify_lemma2"); plt.close(fig)

# ── V3 Theorem 2 ──
RL = [L for L in logs if L["status"] == "reached"]; cfg = M["cfg"]
fig, axs = plt.subplots(2, 3, figsize=(16, 8.6), constrained_layout=True); ax = axs.ravel(); a = ax[0]
for L in RL: a.plot(np.arange(len(L["J"])) * cf.T, L["J"], lw=0.6, color=COL["ra"], alpha=0.35)
a.plot(np.arange(len(Rr["J"])) * cf.T, Rr["J"], lw=2.2, color="k", label="shared $x_0$"); a.set_yscale("log"); a.set_xlabel("time (s)"); a.set_ylabel(r"$J^*_k$"); a.legend(fontsize=8); a.grid(alpha=0.3, which="both")
nin = sum(1 for L in logs if L["status"] == "infeasible" and len(L["J"]) > 0); a.set_title(f"(a) optimal cost, {len(RL)} runs feasible at $k = 0$:\n{nin} infeasible solves afterwards, {len(RL)} reach $\\mathcal{{X}}_T$", fontsize=9.5)
ell = np.concatenate([L["ell"][:-1] for L in RL if len(L["J"]) > 1]); dJ = np.concatenate([-np.diff(L["J"]) for L in RL if len(L["J"]) > 1]); a = ax[1]
a.scatter(ell, dJ, s=4, c=COL["ra"], alpha=0.35); lim = [0.5 * ell.min(), ell.max() * 3]; a.plot(lim, lim, "k--", lw=1, label=r"$J^*_k - J^*_{k+1} = \ell(x_k, u_k)$"); a.axhline(ra.ELL_LOWER, color="r", ls=":", label=rf"$\ell_{{\min}}$ = {ra.ELL_LOWER:g}")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(lim); a.set_xlabel(r"stage cost $\ell(x_k, u_k)$"); a.set_ylabel(r"decrease $J^*_k - J^*_{k+1}$"); a.legend(fontsize=8, loc="upper left"); a.grid(alpha=0.3, which="both")
a.set_title(f"(b) cost decrease over {dJ.size} steps: smallest {dJ.min():.2f},\n{int((dJ < ell * (1 - 1e-4)).sum())} below $\\ell(x_k, u_k)$, {int((dJ < ra.ELL_LOWER).sum())} below $\\ell_{{\\min}}$", fontsize=9.5)
a = ax[2]
for L in RL: a.step(np.arange(len(L["N"])) * cf.T, L["N"], where="post", lw=0.6, color=COL["ra"], alpha=0.35)
a.step(np.arange(len(Rr["N"])) * cf.T, Rr["N"], where="post", lw=2.2, color="k"); a.set_xlabel("time (s)"); a.set_ylabel("horizon $N$ used at step $k$"); a.grid(alpha=0.3)
inc = sum(int(np.any(np.diff(L["N"]) > 0)) for L in RL); a.set_title(f"(c) adaptive horizon: nonincreasing in {len(RL) - inc}/{len(RL)} runs, smallest {min(L['N'].min() for L in RL)}", fontsize=9.5)
a = ax[3]; hm = np.sort([L["h_min"] for L in RL]); a.plot(np.arange(1, hm.size + 1), hm, ".", color=COL["ra"]); a.axhline(0, color="r", lw=1); a.set_yscale("symlog", linthresh=0.01)
a.set_xlabel("run (sorted)"); a.set_ylabel(r"$\min_t\, \psi(y(t))$ before entering $\mathcal{X}_T$"); a.grid(alpha=0.3, which="both"); a.set_title(f"(d) safety of the continuous closed loop:\nsmallest value {hm.min():.4f} > 0", fontsize=9.5)
a = ax[4]; te = np.array([L["check"]["t_entry"] for L in RL]); tb = np.array([L["check"]["t_bound"] for L in RL]); a.scatter(tb, te, s=10, c=COL["ra"]); lim = [0.01, tb.max() * 2]; a.plot(lim, lim, "k--", lw=1, label="entry time = bound")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(lim); a.set_ylim(lim); a.set_xlabel(r"bound $(\lfloor J^*_0/\ell_{\min}\rfloor + 1)\,T$ (s)"); a.set_ylabel("entry time into $\\mathcal{X}_T$ (s)"); a.legend(fontsize=8, loc="upper left"); a.grid(alpha=0.3, which="both")
a.set_title(f"(e) Corollary 1: the bound holds in {int((te <= tb + 1e-9).sum())}/{len(RL)} runs", fontsize=9.5)
a = ax[5]; up = np.array([np.abs(L["u"]).max(0) for L in RL]); a.scatter(up[:, 0], up[:, 1], s=10, c=COL["ra"], label="RA-MPC, peak per run"); upu = np.array([u["upeak"] for u in M["ura"]])[ura]; a.scatter(upu[:, 0], upu[:, 1], s=10, c=COL["ura"], label="URA, peak per run")
a.axvline(cf.U_MAX[0], color="k", ls="--"); a.axhline(cf.U_MAX[1], color="k", ls="--"); a.set_xlabel(r"peak $|\tau_1|$ (N m)"); a.set_ylabel(r"peak $|\tau_2|$ (N m)"); a.legend(fontsize=8); a.grid(alpha=0.3); a.set_title("(f) peak torques against the bounds", fontsize=9.5)
save(fig, "verify_theorem2"); plt.close(fig)

# ── V4 sweep ──
if SW is not None:
    rows = SW["rows"]; Tsw = np.array([r["T"] for r in rows]); ns = rows[0]["ra"].size; TS = L2["Tstar2"] if L2 else cf.T
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3), constrained_layout=True); a = ax[0]; a.axvspan(Tsw.min() * 0.8, TS, color="g", alpha=0.10, label=rf"second-order certified range $T < T^\star$ = {1e3*TS:.1f} ms")
    a.plot(Tsw, [100 * r["ra"].mean() for r in rows], "s-", color=COL["ra"], lw=2, label="RA-MPC"); a.plot(Tsw, [100 * r["fv"].mean() for r in rows], "o-", color=COL["van"], lw=1.5, label="vanilla MPC")
    a.plot(Tsw, [100 * r["zoh"].mean() for r in rows], "^-", color=COL["ura"], lw=1.5, label="terminal controller $u^c$ alone, held torque"); a.axvline(cf.T, color="r", ls=":"); a.set_xscale("log"); a.set_ylim(0, 100)
    a.set_xlabel("period $T$ (s),  $N = 0.2/T$"); a.set_ylabel(r"success rate $\rho$ (%)"); a.legend(fontsize=8, loc="lower left"); a.grid(alpha=0.3, which="both"); a.set_title(f"(a) success rate on {ns} initial states, terminal design fixed", fontsize=9.5)
    a = ax[1]; f0 = np.array([ns - r["inf0"] for r in rows]); a.plot(Tsw, [100 * r["ra"].sum() / max(f, 1) for r, f in zip(rows, f0)], "s-", color=COL["ra"], lw=2, label="RA-MPC success among the states feasible at $k = 0$ (%)")
    a.plot(Tsw, [100 * r["later"] / max(f, 1) for r, f in zip(rows, f0)], "x-", color="r", label="runs with an infeasible solve after $k = 0$ (%)"); a.plot(Tsw, [100 * r["dec_viol"] / max(r["steps"], 1) for r in rows], "d-", color="0.4", label=r"steps violating $J^*_{k+1} \leq J^*_k - \ell_{\min}$ (%)")
    a.axvline(cf.T, color="r", ls=":"); a.set_xscale("log"); a.set_ylim(-3, 104); a.set_xlabel("period $T$ (s)"); a.set_ylabel("%"); a.legend(fontsize=8, loc="center left"); a.grid(alpha=0.3, which="both"); a.set_title("(b) the guarantees of Theorem 2 against the period", fontsize=9.5)
    save(fig, "verify_sweep"); plt.close(fig)

# ── paper figures (style of make_paper_figs_manipulator.ipynb) ──
px = 1 / plt.rcParams["figure.dpi"]; FS_P = 30; C_ARM = (146 / 255, 149 / 255, 145 / 255); C_START = (1.0, 0.549, 0.0); C_END = (0.102, 0.596, 0.314); C_GRAY = (178 / 255, 187 / 255, 205 / 255, 0.45)


def draw_robot(ax, x):
    q1, q2 = float(x[0]), float(x[1]); jx, jy = 4 * np.cos(q1), 4 * np.sin(q1); ex, ey = jx + 4 * np.cos(q1 + q2), jy + 4 * np.sin(q1 + q2)
    ax.plot([0, jx], [0, jy], "o-", alpha=0.45, lw=6, color=C_ARM, zorder=6); ax.plot([jx, ex], [jy, ey], "o-", alpha=0.45, lw=6, color=C_ARM, zorder=6); ax.plot(0, 0, "ko", markersize=6, zorder=51)


def paper_figs(V, U, Rr, x0s, suffix=""):
    PAN = {"fig_manip_vanilla": (fv, fk(V["x"].T[:2]).T, False), "fig_manip_uncon": (ura, U["y"], True), "fig_manip_rampc": (rr, fk(Rr["x"].T[:2]).T, True)}
    for name, (mask, yt, field) in PAN.items():
        fig, ax = plt.subplots(figsize=(650 * px, 600 * px), layout="constrained"); fig.set_dpi(200)
        if field: ax.streamplot(yv, yw, KU, KV, color=C_GRAY, linewidth=0.9, density=0.9, arrowsize=0.9, zorder=1)
        ax.contourf(G1, G2, PS, levels=[0, np.inf], colors=["#e6f2ff"], alpha=0.3, zorder=2); ax.contour(G1, G2, PS, levels=[0], colors=["#2166ac"], linewidths=2, zorder=3)
        ax.contourf(G1, G2, PH, levels=[-np.inf, 0], colors=["#fee0d2"], alpha=1.0, zorder=4); ax.contour(G1, G2, PH, levels=[0], colors=["#d6604d"], linewidths=3, zorder=5)
        draw_robot(ax, x0s)
        ax.scatter(*Y0[:, ~mask], s=50, color="#d73027", alpha=0.85, marker="x", linewidths=3, zorder=8); ax.scatter(*Y0[:, mask], s=50, color="#2166ac", alpha=0.35, marker="o", zorder=9)
        ax.plot(yt[:, 0], yt[:, 1], color="black", lw=3, alpha=0.65, linestyle="--", zorder=10); ax.scatter([yt[0, 0]], [yt[0, 1]], s=70, c=[C_START], marker="o", linewidths=0.4, zorder=20); ax.scatter([yt[-1, 0]], [yt[-1, 1]], s=70, c=[C_END], marker="x", linewidths=3, zorder=20)
        ax.text(0.2, 1.8, r"$\mathcal{X}_S$", zorder=100, fontsize=FS_P); ax.text(5.4, 2.4, r"$\mathcal{X}_T$", zorder=100, fontsize=FS_P)
        ax.set_xlabel("$y_1$", fontsize=FS_P - 4); ax.set_ylabel("$y_2$", fontsize=FS_P - 4); ax.xaxis.set_tick_params(labelsize=22); ax.yaxis.set_tick_params(labelsize=22); ax.set_xlim(-3.5, 8.2); ax.set_ylim(-3.0, 3.2); ax.set_aspect("equal")
        save(fig, name + suffix, pdf=True); plt.close(fig)
    FS_C = 24; fig, (axR, axW, axA) = plt.subplots(3, 1, figsize=(650 * px, 600 * px), layout="constrained", gridspec_kw={"height_ratios": [0.7, 1.0, 1.0]}); fig.set_dpi(200)
    if SW is not None:
        rows = SW["rows"]; Tsw = np.array([r["T"] for r in rows]); rho = np.array([100 * r["ra"].mean() for r in rows]); axR.plot(Tsw, rho, color=COL["ra"], marker="s", lw=2.5, ms=9, alpha=0.9)
        axR.set_ylim(5 * np.floor((rho.min() - 3) / 5), 5 * np.ceil((rho.max() + 3) / 5)); axR.set_xlim(Tsw.min() * 0.8, Tsw.max() * 1.25)
    axR.axvline(cf.T, color="gray", lw=1, ls="--", alpha=0.6); axR.set_xscale("log"); axR.set_xlabel("$T$ [s]", fontsize=FS_C - 6); axR.set_ylabel(r"$\rho$ [%]", fontsize=FS_C - 2); axR.xaxis.set_tick_params(labelsize=14); axR.yaxis.set_tick_params(labelsize=14); axR.grid(True, alpha=0.3, which="both")
    tu, uu = U["tt"], U["u"]; tr_, ur = Rr["t"], Rr["u"]
    for ax, ch, ub, lbl in zip((axW, axA), (0, 1), cf.U_MAX, (r"$\tau_1$", r"$\tau_2$")):
        ax.axhspan(-ub, ub, color="#f0f0f0", alpha=0.6, zorder=0); ax.axhline(ub, color="black", lw=1.5, ls="--", alpha=0.85, zorder=5); ax.axhline(-ub, color="black", lw=1.5, ls="--", alpha=0.85, zorder=5); ax.axhline(0, color="gray", lw=0.8, ls=":", zorder=4)
        LOGT = tu[-1] > 6 * tr_[-1]; ax.plot(np.maximum(tu, 1e-3) if LOGT else tu, uu[:, ch], color=COL["ura"], lw=3, alpha=0.9, zorder=10)
        ts_ = np.append(tr_[:-1], tr_[-1]); ax.step(np.maximum(ts_, 1e-3) if LOGT else ts_, np.append(ur[:, ch], ur[-1, ch]), where="post", color=COL["ra"], lw=3, alpha=0.9, zorder=10)
        if LOGT: ax.set_xscale("log"); ax.set_xlim(cf.T / 2, tu[-1] * 1.15)
        else: ax.set_xlim(0, tu[-1] * 1.02)
        ym = max(ub, np.abs(uu[:, ch]).max()); ax.set_ylabel(lbl, fontsize=FS_C); ax.set_ylim(-1.25 * ym, 1.25 * ym); ax.xaxis.set_tick_params(labelsize=14); ax.yaxis.set_tick_params(labelsize=14); ax.grid(True, alpha=0.25, which="both")
    axA.set_xlabel("Time [s]", fontsize=FS_C - 6); save(fig, "fig_manip_compare" + suffix, pdf=True); plt.close(fig)
    tt = Rr["t"][:-1]; fig, (axJ, axD, axN) = plt.subplots(3, 1, figsize=(650 * px, 650 * px), sharex=True, layout="constrained", gridspec_kw={"height_ratios": [1.0, 1.0, 0.8]}); fig.set_dpi(200)
    axJ.plot(tt, Rr["J"], color=COL["ra"], lw=3); axJ.set_yscale("log"); axJ.set_ylabel(r"$J^*_k$", fontsize=20); axJ.grid(True, alpha=0.3, which="both")
    dJ = -np.diff(Rr["J"]); axD.plot(tt[:-1], dJ, color=COL["ra"], lw=3, label=r"$J^*_k - J^*_{k+1}$"); axD.plot(tt[:-1], Rr["ell"][:-1], "o", color="black", ms=2.5, label=r"$\ell(x_k, u_k)$"); axD.axhline(ra.ELL_LOWER, color="black", ls="--", lw=1.2, label=r"$\ell_{\min}$")
    axD.set_yscale("log"); axD.set_ylabel("decrease", fontsize=16); axD.grid(True, alpha=0.3, which="both"); axD.legend(fontsize=10, loc="upper right", ncol=3)
    axN.step(tt, Rr["N"], where="post", color=COL["ra"], lw=3); axN.set_ylabel(r"$N$", fontsize=20); axN.set_xlabel("Time [s]", fontsize=16); axN.grid(True, alpha=0.3); save(fig, "fig_manip_cost" + suffix, pdf=True); plt.close(fig)


paper_figs(V, U, Rr, x0s)
both = ura & rr; i_far = int(np.where(both)[0][np.argmax(d2[both])]); i_van = int(np.where(vr)[0][np.argmax(d2[vr])])
U2 = cf.ura_run(X0[:, i_far], dense=True); R2 = logs[i_far]; V2 = van[i_van]
print(f"far variant: URA and RA-MPC from state {i_far} y0 {np.round(Y0[:, i_far], 3)} (distance {np.sqrt(d2[i_far]):.2f}; URA {U2['t']:.2f} s, peaks {np.round(U2['upeak'], 1)}; RA-MPC {R2['t'][-1]:.2f} s), vanilla from state {i_van} (distance {np.sqrt(d2[i_van]):.2f}, {V2['t'][-1]:.2f} s)")
paper_figs(V2, U2, R2, X0[:, i_far], "_far")
pickle.dump(dict(i0=i0, x0=x0s, t_van=float(V["t"][-1]), t_ura=float(U["t"]), t_ra=float(Rr["t"][-1]), all3=int(all3.sum()), far=dict(i=i_far, i_van=i_van)), open(f"{cf.DATA}/shared_x0.pkl", "wb"))
