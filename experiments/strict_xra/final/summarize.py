"""Print the numbers and markdown tables of PAPER_DUBINS.md from data/*.pkl."""
import pickle, numpy as np
import config as cf
from config import rm, C

M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); L2 = pickle.load(open(f"{cf.DATA}/lemma2.pkl", "rb")); SW = pickle.load(open(f"{cf.DATA}/sweep.pkl", "rb"))
X0, P0, logs, cfg = M["X0"], M["Psi0"], M["logs"], M["cfg"]; n = X0.shape[1]
ura = np.array([u["status"] == "reached" for u in M["ura"]]); fv = M["fv"].astype(bool); ra = np.array([L["status"] == "reached" for L in logs])
inX, inT = P0 >= 0, P0 >= cf.DELTA
print("## setting"); print({k: (round(float(v), 6) if np.isscalar(v) else v) for k, v in cfg.items()})
print(f"speeds of the initial set: [{X0[3].min():.4f}, {X0[3].max():.4f}]")
print("\n## main counts")
print("| subset | states | URA | vanilla MPC | RA-MPC |\n| --- | --- | --- | --- | --- |")
for nm, m in (("all", np.ones(n, bool)), ("inside $\\mathcal{X}_{RA}$", inX), ("inside $\\mathcal{X}_{RA}^{\\delta}$", inT), ("outside $\\mathcal{X}_{RA}$", ~inX)):
    print(f"| {nm} | {m.sum()} | {(ura & m).sum()} | {(fv & m).sum()} | {(ra & m).sum()} |")
print(f"pairwise: RA&~van {(ra & ~fv).sum()}, van&~RA {(fv & ~ra).sum()}, RA&~URA {(ra & ~ura).sum()}, URA&~RA {(ura & ~ra).sum()}, all three {(ura & fv & ra).sum()}, none {(~ura & ~fv & ~ra).sum()}")
zs = M["zoh"]["status"]; print(f"sampled-data URA (T={cf.T}): {(zs == 1).sum()} reached, same states as continuous URA: {np.array_equal(zs == 1, ura)}; URA statuses: { {s: sum(u['status'] == s for u in M['ura']) for s in set(u['status'] for u in M['ura'])} }")
RL = [L for L in logs if L["status"] == "reached"]; IL = [L for L in logs if L["status"] != "reached"]
print("RA-MPC statuses:", {s: sum(L["status"] == s for L in logs) for s in set(L["status"] for L in logs)}, "| infeasible at k=0:", sum(L["status"] == "infeasible" and len(L["J"]) == 0 for L in logs), "later:", sum(L["status"] == "infeasible" and len(L["J"]) > 0 for L in logs))
van = M["van"]; print("vanilla closed loop:", {s: sum(v["status"] == s for v in van.values()) for s in set(v["status"] for v in van.values())}, "min h", min(v["h_min"] for v in van.values() if v["status"] == "reached"))
print("\n## performance on common successes")
tu = np.array([u["t"] for u in M["ura"]]); cu = np.array([u["cost"] for u in M["ura"]]); up = np.array([u["upeak"] for u in M["ura"]])
tr = np.array([L["check"]["t_entry"] if L["status"] == "reached" else np.nan for L in logs]); cr = np.array([cf.T * np.sum(L["ell"]) if L["status"] == "reached" else np.nan for L in logs])
tv = np.full(n, np.nan); cv = np.full(n, np.nan)
for i, v in van.items():
    if v["status"] == "reached":
        tv[i] = v["t"][-1]; cv[i] = cf.T * sum(rm.stage_cost(v["x"][k], v["u"][k]) for k in range(len(v["u"])))
for nm, m in (("URA and RA-MPC", ura & ra), ("all three", ura & ra & fv)):
    print(f"{nm}: {m.sum()} states | time median: URA {np.median(tu[m]):.1f} (max {tu[m].max():.1f}), RA {np.median(tr[m]):.2f} (max {np.nanmax(tr[m]):.2f})" + (f", vanilla {np.nanmedian(tv[m]):.2f}" if nm == "all three" else "")
          + f" | cost median: URA {np.median(cu[m]):.1f}, RA {np.median(cr[m]):.1f}" + (f", vanilla {np.nanmedian(cv[m]):.1f}" if nm == "all three" else "")
          + f" | median ratio time {np.median(tu[m] / tr[m]):.1f}, cost {np.median(cu[m] / cr[m]):.1f}")
print(f"peaks: URA on its successes {np.round(up[ura].max(0), 4)}, all URA runs {np.round(up.max(0), 4)}; RA-MPC {np.round(np.max([np.abs(L['u']).max(0) for L in RL], axis=0), 4)}; vanilla {np.round(np.max([np.abs(v['u']).max(0) for v in van.values() if len(v['u'])], axis=0), 4)}")
sat = np.concatenate([np.abs(L["u"]) >= 0.98 * np.array(cf.U_MAX) for L in RL]); print(f"RA-MPC share of applied inputs within 2% of the bound: omega {100*sat[:,0].mean():.0f}%, a {100*sat[:,1].mean():.0f}%; RA entry time max {np.nanmax(tr):.2f}")
ms = np.concatenate([L["solve_ms"] for L in RL]); mv = np.concatenate([v["solve_ms"] for v in van.values()])
print(f"solve time per step [ms]: RA-MPC mean {ms.mean():.1f} median {np.median(ms):.1f} (first solve {np.mean([L['solve_ms'][0] for L in RL]):.1f}, later {np.mean(np.concatenate([L['solve_ms'][1:] for L in RL])):.1f}, n={ms.size}); vanilla mean {mv.mean():.1f} median {np.median(mv):.1f} (n={mv.size})")
print("\n## Theorem 2 / Corollary 1")
steps = sum(len(L["J"]) for L in RL); ell = np.concatenate([L["ell"][:-1] for L in RL if len(L["J"]) > 1]); dJ = np.concatenate([-np.diff(L["J"]) for L in RL if len(L["J"]) > 1])
te = np.array([L["check"]["t_entry"] for L in RL]); tb = np.array([L["check"]["t_bound"] for L in RL]); J0 = np.array([L["J"][0] for L in RL])
print(f"runs feasible at k=0: {sum(len(L['J']) > 0 for L in logs)}; reached {len(RL)}; solves {steps}; candidate kept {sum(int(np.sum(L['used_candidate'])) for L in RL)}; u^c of the candidate within bounds: {all(np.all(L['uc_in_bounds']) for L in RL)}")
print(f"decrease: {dJ.size} steps, min {dJ.min():.3f}, min ratio dJ/ell {np.min(dJ/ell):.7f}, below ell(1-1e-4): {(dJ < ell*(1-1e-4)).sum()}, below ell_min: {(dJ < rm.ELL_LOWER).sum()}; ell range [{ell.min():.2f}, {ell.max():.1f}]")
print(f"horizon: nonincreasing {sum(int(np.all(np.diff(L['N']) <= 0)) for L in RL)}/{len(RL)}, shrinks in {sum(int(L['N'].min() < cf.N) for L in RL)} runs, min {min(L['N'].min() for L in RL)}, median final {np.median([L['N'][-1] for L in RL])}")
print(f"safety: min h over all closed loops {min(L['h_min'] for L in RL):.4f}; runs with h_min<0: {sum(L['h_min'] < 0 for L in RL)}")
print(f"entry time: median {np.median(te):.2f}, max {te.max():.2f}; bound: min {tb.min():.2f}, median {np.median(tb):.0f}, max {tb.max():.0f}; holds {(te <= tb + 1e-9).sum()}/{len(RL)}; J0 range [{J0.min():.1f}, {J0.max():.3g}]; median bound/entry {np.median(tb/te):.0f}")
vf = L2["vf_check"]; print(f"condition (22) on {vf['n']} sampled states (every 8th stored): min omega*dPsi - ell = {(cfg['omega']*vf['dpsi'] - vf['ell']).min():.3f}, violations {(cfg['omega']*vf['dpsi'] < vf['ell']).sum()}, sup ell {vf['sup_ell']:.2f}, eps_big {L2['eps_big']:.3e}, omega*eps_big {cfg['omega']*L2['eps_big']:.2f}")
print("\n## Lemma 2")
print(f"samples: {L2['n_delta']} in X_RA^delta, big {L2['n_big']}; eta {L2['eta']:.4e} (lambda*delta {L2['lam_delta']:.4e}); M2 {L2['M2_ref']:.4f}; T* = 2eta/M2 = {1e3*L2['Tstar2']:.1f} ms (self-consistent {1e3*L2['Tstar2_self']:.1f} ms); paper bound T* = {L2['Tstar_paper']:.3e} s (M_g {np.round(L2['Mg'],3)}, l_uc {np.round(L2['lu'],2)}, M_v {L2['Mv']:.3f})")
print(f"u^c on X_RA^delta: {np.round(L2['umax_delta'],4)}, on X_RA: {np.round(L2['umax_0'],4)}; speeds in X_RA^delta {np.round(L2['speed_delta'],4)}")
print("| $T$ [s] | arcs | smallest increase of $\\Psi$ | increasing | end in $\\mathcal{X}_{RA}^{\\delta}$ | arc inside $\\mathcal{X}_{RA}$ | second-order bound |\n| --- | --- | --- | --- | --- | --- | --- |")
for r in L2["onestep"]:
    print(f"| {r['T']:g} | {r['n']} | {r['dmin']:+.2e} | {100*r['frac_inc']:.2f} % | {100*r['frac_inv']:.2f} % | {100*r['frac_in_xra']:.2f} % | {r['lb2']:+.2e} |")
print("| $T$ [s] | from $\\mathcal{X}_{RA}^{\\delta}$: reach | $\\Psi$ increasing at every sample | from $\\mathcal{X}_{RA}$: reach | $\\Psi$ increasing at every sample | median time |\n| --- | --- | --- | --- | --- | --- |")
for r in L2["closed_loop"]:
    d, z = r["delta"], r["zero"]; print(f"| {r['T']:g} | {d['reached']}/{d['n']} | {d['monotone']}/{d['n']} | {z['reached']}/{z['n']} | {z['monotone']}/{z['n']} | {d['t_med']:.1f} s |")
print("\n## sweep"); print(f"omega fixed at {SW['omega']:.4g}")
print("| $T$ [s] | $N$ | vanilla | sampled URA | RA-MPC | feasible at $k=0$ | infeasible later | unsafe | timeout | decrease violations | $\\min\\Delta\\Psi$ on $\\mathcal{X}_{RA}^{\\delta}$ | median entry | ms/solve |\n| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for r in SW["rows"]:
    print(f"| {r['T']:g} | {r['N']} | {r['fv'].sum()} | {r['zoh'].sum()} | {r['ra'].sum()} | {n - r['inf0']} | {r['later']} | {r['unsafe']} | {r['timeout']} | {r['dec_viol']}/{r['steps']} | {r['eps']:+.1e} | {r['t_med']:.2f} s | {r['solve_ms']:.0f} |")
base = [r for r in SW["rows"] if abs(r["T"] - cf.T) < 1e-9][0]; print("sweep base row equals main:", np.array_equal(base["ra"], ra), np.array_equal(base["fv"], fv))
for r in SW["rows"]:
    print(f"T={r['T']:g}: vanilla&~RA {(r['fv'] & ~r['ra']).sum()}, RA&~vanilla {(r['ra'] & ~r['fv']).sum()}, status {r['status']}")
try:
    S = pickle.load(open(f"{cf.DATA}/shared_x0.pkl", "rb")); print("\nshared x0:", S)
except Exception:
    pass
