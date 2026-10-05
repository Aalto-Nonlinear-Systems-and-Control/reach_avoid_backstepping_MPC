"""Outcome maps for the two initial sets (S_RA inside X_RA, S_0 at rest) of setting R2."""
import os, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from eval_init import C, DELTA
import ra_mpc_dubins as rm

tag, N = (sys.argv[1], int(sys.argv[2])) if len(sys.argv) > 2 else ("3_3", 25)
g = np.linspace(-2.1, 2.1, 601); Y1, Y2 = np.meshgrid(g, g)
PS = np.broadcast_to(C.psi(Y1, Y2), Y1.shape) - C.shift
SAFE = rm.safe_np(Y1, Y2, 0 * Y1, 0 * Y1); TGT = rm.tgt_np(Y1, Y2, 0 * Y1, 0 * Y1)
PS = np.where(SAFE >= 0, PS, np.nan)
fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.6), constrained_layout=True)
for a, grp, ttl in zip(ax, ("I", "II"), (r"$S_{RA}$: velocity inside $\mathcal{X}_{RA}$", r"$S_0$: at rest, random heading")):
    f = f"logs/init2_{grp}_{tag}_N{N}.npz"
    d = np.load(f if os.path.exists(f) else f.replace("init2_", "init3_")); X = d["X0"]; ra = d["ra"]; fv = d["fv"]
    a.contourf(Y1, Y2, PS, levels=[0, DELTA, 10], colors=["#dfe9f5", "#b9d3ee"], alpha=0.9)
    a.contour(Y1, Y2, SAFE, levels=[0], colors="k", linewidths=1.0)
    a.contourf(Y1, Y2, -TGT, levels=[0, 1e9], colors=["#bfe5bf"]); a.contour(Y1, Y2, TGT, levels=[0], colors="g", linewidths=1.0)
    both = ra & fv; only = ra & ~fv; none = ~ra & ~fv; vo = fv & ~ra
    a.scatter(*X[:2, both], s=13, c="tab:green", label=f"vanilla and RA-MPC ({both.sum()})")
    a.scatter(*X[:2, only], s=13, c="tab:blue", label=f"RA-MPC only ({only.sum()})")
    if vo.any(): a.scatter(*X[:2, vo], s=22, c="tab:orange", marker="s", label=f"vanilla only ({vo.sum()})")
    lab = f"neither MPC ({none.sum()})" + (r", $u^c$ alone reaches" if grp == "I" else "")
    a.scatter(*X[:2, none], s=20, c="tab:red", marker="x", linewidths=1.0, label=lab)
    a.set_title(f"{ttl}  ({X.shape[1]} states)", fontsize=10); a.set_aspect("equal"); a.set_xlim(-2.15, 2.15); a.set_ylim(-2.15, 2.15)
    a.set_xlabel("$y_1$"); a.set_ylabel("$y_2$"); a.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
fig.suptitle(rf"R2, $(\omega_{{\max}}, a_{{\max}}) = ({tag.replace('_', ', ')})$, $N = {N}$, $T = 0.05$ s.  Light blue: $\psi \geq \delta_s/\lambda$ (positions of $\mathcal{{X}}_{{RA}}$), darker: positions of $\mathcal{{X}}_{{RA}}^{{\delta}}$", fontsize=9)
fig.savefig(f"figs/init_sets_{tag}_N{N}.png", dpi=170); print("saved", f"figs/init_sets_{tag}_N{N}.png")
