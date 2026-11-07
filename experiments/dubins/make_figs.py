"""Figures of Example 1 in the paper: fig_dubins_vanilla, fig_dubins_uncon, fig_dubins_rampc (Fig. 4) and fig_dubins_compare (Fig. 5). Output: figs/."""
import os, pickle, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as cf
from config import C, rm

M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); L2 = pickle.load(open(f"{cf.DATA}/lemma2.pkl", "rb")); SW = pickle.load(open(f"{cf.DATA}/sweep.pkl", "rb"))
X0, P0, logs = M["X0"], M["Psi0"], M["logs"]; n = X0.shape[1]
# URA baseline: the same k1 at cf.URA_SCALE, stored by run_main.py
URA_LIST = M["ura"]
ura = np.array([u["status"] == "reached" for u in URA_LIST]); fv = M["fv"].astype(bool); ra = np.array([L["status"] == "reached" for L in logs])
ura_viol = np.array([(np.asarray(u["upeak"]) > np.asarray(cf.U_MAX)).any() for u in URA_LIST])
TS = L2["Tstar2"]; COL = dict(van="#e08214", ura="#2166ac", ra="#8856a7")
g = np.linspace(-2.1, 2.1, 400); G1, G2 = np.meshgrid(g, g)
Zs = rm.safe_np(G1, G2, 0 * G1, 0 * G1); Zt = rm.tgt_np(G1, G2, 0 * G1, 0 * G1)
save = lambda fig, name: fig.savefig(os.path.join(cf.FIGS, f"{name}.png"), dpi=200) or print("saved", name)

# ── shared initial state: farthest state from which all three methods succeed ─────────────────────────────
van_ok = np.array([i in M["van"] and M["van"][i]["status"] == "reached" for i in range(n)])
all3 = ura & ra & fv & van_ok
pref = all3 & ura_viol if (all3 & ura_viol).any() else all3          # prefer a start where URA exceeds the bounds
d2 = ((X0[:2] - rm.TARGET_POS[:, None])**2).sum(0); i0 = int(np.where(pref)[0][np.argmax(d2[pref])]); x0s = X0[:, i0]
U = cf.ura_run(x0s, dense=True, scale=cf.URA_SCALE); V = M["van"][i0]; Rr = logs[i0]
print(f"shared x0 = state {i0}: {np.round(x0s, 4)}; all-three successes {all3.sum()} (URA exceeding the bounds among them: {(all3 & ura_viol).sum()}); "
      f"times: vanilla {V['t'][-1]:.2f} s, URA {U['t']:.2f} s (peak inputs {np.round(U['upeak'], 2)}), RA-MPC {Rr['t'][-1]:.2f} s")


rows = SW["rows"]; Tsw = np.array([r["T"] for r in rows])

# ── paper figures ─────────────────────────────────────
px = 1 / plt.rcParams["figure.dpi"]; FS_P = 30
C_FIELD, FIELD_LW = "0.55", 1.2; C_START = (1.0, 0.549, 0.0); C_END = (0.102, 0.596, 0.314)
yv = np.linspace(-2.1, 2.1, 30); VY1, VY2 = np.meshgrid(yv, yv); kk = C.k1(VY1, VY2)
KU = np.broadcast_to(np.asarray(kk[0], float), VY1.shape).astype(float); KV = np.broadcast_to(np.asarray(kk[1], float), VY1.shape).astype(float)
def paper_figs(V, U, Rr):
  PAN = {"fig_dubins_vanilla": (fv, V["x"][:, :2], False), "fig_dubins_uncon": (ura, U["x"][:, :2], True), "fig_dubins_rampc": (ra, Rr["x"][:, :2], True)}
  for name, (mask, yt, field) in PAN.items():
      fig, ax = plt.subplots(figsize=(650 * px, 600 * px), layout="constrained"); fig.set_dpi(200)
      if field:
          # field k1: opaque mid grey, drawn above the safe-set fill and below its boundary so that it survives printing
          ax.streamplot(yv, yv, KU, KV, color=C_FIELD, linewidth=FIELD_LW, density=0.75, arrowsize=1.5, zorder=1.5)
      ax.contourf(G1, G2, Zs, levels=[0, np.inf], colors=["#e6f2ff"], alpha=0.3, zorder=1); ax.contour(G1, G2, Zs, levels=[0], colors=["#2166ac"], linewidths=2, zorder=2)
      ax.contourf(G1, G2, Zt, levels=[-np.inf, 0], colors=["#fee0d2"], alpha=1.0, zorder=3); ax.contour(G1, G2, Zt, levels=[0], colors=["#d6604d"], linewidths=2, zorder=4)
      ax.scatter(*X0[:2, ~mask], s=50, color="#d73027", alpha=0.85, marker="x", linewidths=3, zorder=8)
      ax.scatter(*X0[:2, mask], s=50, color="#2166ac", alpha=0.35, marker="o", zorder=9)
      ax.plot(yt[:, 0], yt[:, 1], color="black", lw=3, alpha=0.65, linestyle="--", zorder=10)
      ax.scatter([yt[0, 0]], [yt[0, 1]], s=70, c=[C_START], marker="o", linewidths=0.4, zorder=20)
      ax.scatter([yt[-1, 0]], [yt[-1, 1]], s=70, c=[C_END], marker="x", linewidths=3, zorder=20)
      ax.text(0, 1.5, r"$\mathcal{X}_S$", zorder=100, fontsize=FS_P); ax.text(-1.9, 0, r"$\mathcal{X}_T$", zorder=100, fontsize=FS_P)
      ax.set_xlabel("$y_1$", fontsize=FS_P - 4); ax.set_ylabel("$y_2$", fontsize=FS_P - 4)
      ax.xaxis.set_tick_params(labelsize=22); ax.yaxis.set_tick_params(labelsize=22); ax.set_xlim(-2.1, 2.1); ax.set_ylim(-2.1, 2.1); ax.set_aspect("equal")
      save(fig, name); plt.close(fig)

  FS_C = 24
  fig, (axR, axW, axA) = plt.subplots(3, 1, figsize=(650 * px, 600 * px), layout="constrained", gridspec_kw={"height_ratios": [0.7, 1.0, 1.0]}); fig.set_dpi(200)
  axR.axvspan(Tsw.min() * 0.8, TS, color="g", alpha=0.12)
  axR.plot(Tsw, [100 * r["ra"].mean() for r in rows], color=COL["ra"], marker="s", lw=2.5, ms=9, alpha=0.9)
  rho = np.array([100 * r["ra"].mean() for r in rows]); axR.set_ylim(5 * np.floor((rho.min() - 3) / 5), 5 * np.ceil((rho.max() + 3) / 5))
  axR.axvline(cf.T, color="0.2", lw=1.6, ls="--", alpha=0.9, zorder=2); axR.set_xscale("log"); axR.set_xlim(Tsw.min() * 0.8, Tsw.max() * 1.25)
  # mark the success rate at the sampling period of the main run (paper: Fig. 5 / Fig. 8)
  iT = int(np.argmin(np.abs(Tsw - cf.T))); ylo, yhi = axR.get_ylim()
  axR.plot(Tsw[iT], rho[iT], marker="s", ms=14, mfc="none", mec=COL["ra"], mew=2.5, zorder=5)
  axR.annotate(rf"$\rho = {rho[iT]:.1f}\,\%$", xy=(Tsw[iT], rho[iT]), xytext=(cf.T * 1.7, ylo + 0.2 * (yhi - ylo)), fontsize=15, color=COL["ra"], ha="left", va="center",
               arrowprops=dict(arrowstyle="-", color="0.3", lw=1.0, shrinkA=0, shrinkB=7), bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9), zorder=6)
  axR.set_xlabel("$T$ [s]", fontsize=FS_C - 6); axR.set_ylabel(r"$\rho$ [%]", fontsize=FS_C - 2)
  axR.xaxis.set_tick_params(labelsize=14); axR.yaxis.set_tick_params(labelsize=14); axR.grid(True, alpha=0.3, which="both")
  tu, uu = U["tt"], U["u"]; tr_, ur = Rr["t"], Rr["u"]
  for ax, ch, ub, lbl in zip((axW, axA), (0, 1), cf.U_MAX, (r"$\omega$", r"$a$")):
      ax.axhspan(-ub, ub, color="#f0f0f0", alpha=0.6, zorder=0); ax.axhline(ub, color="black", lw=1.5, ls="--", alpha=0.85, zorder=5); ax.axhline(-ub, color="black", lw=1.5, ls="--", alpha=0.85, zorder=5)
      ax.axhline(0, color="gray", lw=0.8, ls=":", zorder=4)
      LOGT = tu[-1] > 6 * tr_[-1]                      # logarithmic time axis only when the two durations differ a lot
      ax.plot(np.maximum(tu, 1e-3) if LOGT else tu, uu[:, ch], color=COL["ura"], lw=3, alpha=0.9, zorder=10)
      ts_ = np.append(tr_[:-1], tr_[-1]); ax.step(np.maximum(ts_, 1e-3) if LOGT else ts_, np.append(ur[:, ch], ur[-1, ch]), where="post", color=COL["ra"], lw=3, alpha=0.9, zorder=10)
      if LOGT:
          ax.set_xscale("log"); ax.set_xlim(cf.T / 2, tu[-1] * 1.15)
      else:
          ax.set_xlim(0, tu[-1] * 1.02)
      ym = max(ub, np.abs(uu[:, ch]).max()); ax.set_ylabel(lbl, fontsize=FS_C); ax.set_ylim(-1.25 * ym, 1.25 * ym)
      ax.xaxis.set_tick_params(labelsize=14); ax.yaxis.set_tick_params(labelsize=14); ax.grid(True, alpha=0.25, which="both")
  axA.set_xlabel("Time [s]", fontsize=FS_C - 6)
  save(fig, "fig_dubins_compare"); plt.close(fig)


paper_figs(V, U, Rr)
pickle.dump(dict(i0=i0, x0=x0s, t_van=float(V["t"][-1]), t_ura=float(U["t"]), t_ra=float(Rr["t"][-1]), all3=int(all3.sum())), open(f"{cf.DATA}/shared_x0.pkl", "wb"))
