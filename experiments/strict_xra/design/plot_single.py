"""Outcome maps of URA, vanilla MPC and RA-MPC on the single initial set, plus time to the target."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from eval_init import C, DELTA
import ra_mpc_dubins as rm

r0, tag, N = (sys.argv[1], sys.argv[2], int(sys.argv[3])) if len(sys.argv) > 3 else ("0.03", "3_3", 25)
d = np.load(f"logs/single_r{r0}_{tag}_N{N}.npz"); X = d["X0"]; P = d["Psi"]
g = np.linspace(-2.1, 2.1, 601); Y1, Y2 = np.meshgrid(g, g)
SAFE = rm.safe_np(Y1, Y2, 0 * Y1, 0 * Y1); TGT = rm.tgt_np(Y1, Y2, 0 * Y1, 0 * Y1)
PS = np.where(SAFE >= 0, np.broadcast_to(C.psi(Y1, Y2), Y1.shape) - C.shift, np.nan)
fig, ax = plt.subplots(1, 4, figsize=(17, 4.4), constrained_layout=True)
for a, key, ttl, col in zip(ax, ("fv", "ura", "ra"), ("vanilla MPC", r"URA ($u^c$ alone)", "RA-MPC"), ("tab:orange", "tab:blue", "tab:purple")):
    ok = d[key].astype(bool)
    a.contourf(Y1, Y2, PS, levels=[0, DELTA, 10], colors=["#e6eef7", "#c9dcf0"], alpha=0.9)
    a.contour(Y1, Y2, SAFE, levels=[0], colors="k", linewidths=1.0)
    a.contourf(Y1, Y2, -TGT, levels=[0, 1e9], colors=["#bfe5bf"]); a.contour(Y1, Y2, TGT, levels=[0], colors="g", linewidths=1.0)
    a.scatter(*X[:2, ok], s=14, c=col); a.scatter(*X[:2, ~ok], s=18, c="0.45", marker="x", linewidths=0.9)
    a.scatter(*X[:2, P >= 0], s=42, facecolors="none", edgecolors="k", linewidths=0.4)
    a.set_title(f"{ttl}: {ok.sum()}/{X.shape[1]}", fontsize=11); a.set_aspect("equal"); a.set_xlim(-2.15, 2.15); a.set_ylim(-2.15, 2.15)
    a.set_xlabel("$y_1$"); a.set_ylabel("$y_2$")
both = d["ura"].astype(bool) & d["ra"].astype(bool); a = ax[3]
a.scatter(d["t_ura"][both], d["te"][both], s=12, c="tab:purple")
a.set_xscale("log"); a.set_yscale("log"); lim = [0.05, 200]; a.plot(lim, lim, "k--", lw=0.8); a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("time to target, URA (s)"); a.set_ylabel("time to target, RA-MPC (s)")
a.set_title(f"{both.sum()} common successes: median {np.median(d['t_ura'][both]):.1f} s vs {np.median(d['te'][both]):.2f} s", fontsize=10); a.grid(alpha=0.3, which="both")
fig.suptitle(rf"One initial set: positions uniform in $\mathcal{{X}}_S \setminus \mathcal{{X}}_T$, $\|\dot y(0) - k_1(y(0))\| \leq {r0}$; "
             rf"$(\omega_{{\max}}, a_{{\max}}) = ({tag.replace('_', ', ')})$, $N = {N}$.  Circled: the {int((P >= 0).sum())} states inside $\mathcal{{X}}_{{RA}}$; shading: positions of $\mathcal{{X}}_{{RA}}$ and $\mathcal{{X}}_{{RA}}^{{\delta}}$", fontsize=9.5)
out = f"figs/single_r{r0}_{tag}_N{N}.png"; fig.savefig(out, dpi=150); print("saved", out)
