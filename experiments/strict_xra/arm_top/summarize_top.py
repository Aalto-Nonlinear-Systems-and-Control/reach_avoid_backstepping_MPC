"""Numbers and tables for the arm document."""
import os, pickle, numpy as np
import config_top as cf
from config_top import C, fk, ra
M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0, P0, logs, van, cfg = M["X0"], M["Psi0"], M["logs"], M["van"], M["cfg"]; n = X0.shape[1]
ura = np.array([u["status"] == "reached" for u in M["ura"]]); uc = np.array([u["status"] == "reached" for u in M["uc"]]); fv = np.array([len(v["J"]) > 0 for v in van]); vr = np.array([v["status"] == "reached" for v in van]); rr = np.array([L["status"] == "reached" for L in logs])
up = np.array([u["upeak"] for u in M["ura"]]); viol = (up > np.asarray(cf.U_MAX)).any(1); zs = M["zoh"]["status"]
print("## setting", {k: (round(float(v), 6) if np.isscalar(v) and not isinstance(v, str) else v) for k, v in cfg.items()}, "| scale.json:", open(f"{cf.HERE}/data/scale_{cf.TASK}.json").read().replace("\n", " "))
print(f"\n## counts: URA {ura.sum()} (exceeding the bounds {viol.sum()}, among successes {(viol&ura).sum()}), terminal controller alone continuous {uc.sum()} sampled {(zs==1).sum()} | vanilla feasible at k=0 {fv.sum()}, reached {vr.sum()} | RA-MPC {rr.sum()} (infeasible at k=0 {sum(len(L['J'])==0 for L in logs)}, later {sum(L['status']=='infeasible' and len(L['J'])>0 for L in logs)}) | in X_RA {(P0>=0).sum()}, terminal set {(P0>=cf.DELTA).sum()}")
print(f"pairwise: RA&~van {(rr&~vr).sum()}, van&~RA {(vr&~rr).sum()}, RA&~URA {(rr&~ura).sum()}, URA&~RA {(ura&~rr).sum()}, all three {(ura&vr&rr).sum()}, none {(~ura&~vr&~rr).sum()}")
tu = np.array([u["t"] for u in M["ura"]]); cu = np.array([u["cost"] for u in M["ura"]]); tr = np.array([L["check"]["t_entry"] if L["status"] == "reached" else np.nan for L in logs]); cr = np.array([cf.T * np.sum(L["ell"]) if L["status"] == "reached" else np.nan for L in logs])
tv = np.array([v["t"][-1] if v["status"] == "reached" else np.nan for v in van]); cv = np.array([cf.T * np.sum(v["ell"]) if v["status"] == "reached" else np.nan for v in van])
for nm, m in (("URA and RA-MPC", ura & rr), ("all three", ura & rr & vr)):
    if m.any(): print(f"{nm}: {m.sum()} states | time med URA {np.median(tu[m]):.2f} (max {tu[m].max():.2f}), RA {np.nanmedian(tr[m]):.2f} (max {np.nanmax(tr[m]):.2f})" + (f", vanilla {np.nanmedian(tv[m]):.2f}" if nm == "all three" else "") + f" | cost med URA {np.median(cu[m]):.1f}, RA {np.nanmedian(cr[m]):.1f}" + (f", vanilla {np.nanmedian(cv[m]):.1f}" if nm == "all three" else ""))
RL = [L for L in logs if L["status"] == "reached"]; VL = [v for v in van if v["status"] == "reached"]
print(f"peaks: URA successes {np.round(up[ura].max(0),1)}, all {np.round(up.max(0),1)}; RA-MPC {np.round(np.max([np.abs(L['u']).max(0) for L in RL], axis=0),1)}; vanilla {np.round(np.max([np.abs(v['u']).max(0) for v in VL], axis=0),1)}")
sat = np.concatenate([np.abs(L["u"]) >= 0.98 * np.array(cf.U_MAX) for L in RL]); print(f"RA-MPC inputs within 2% of the bound: tau1 {100*sat[:,0].mean():.0f}%, tau2 {100*sat[:,1].mean():.0f}%; entry time median {np.nanmedian(tr):.3f} max {np.nanmax(tr):.3f}")
ms = np.concatenate([L["solve_ms"] for L in RL]); mv = np.concatenate([v["solve_ms"] for v in VL]); print(f"solve ms (parallel load): RA mean {ms.mean():.0f} (first {np.mean([L['solve_ms'][0] for L in RL]):.0f}), vanilla mean {mv.mean():.0f}")
print("\n## Theorem 2 / Corollary 1")
ell = np.concatenate([L["ell"][:-1] for L in RL if len(L["J"]) > 1]); dJ = np.concatenate([-np.diff(L["J"]) for L in RL if len(L["J"]) > 1]); te = np.array([L["check"]["t_entry"] for L in RL]); tb = np.array([L["check"]["t_bound"] for L in RL])
print(f"runs feasible at k=0 {sum(len(L['J'])>0 for L in logs)}, reached {len(RL)}, solves {sum(len(L['J']) for L in RL)}, candidate kept {sum(int(np.sum(L['used_candidate'])) for L in RL)}, u^c within bounds {all(np.all(L['uc_in_bounds']) for L in RL)}")
print(f"decrease: {dJ.size} steps, min {dJ.min():.3f}, min ratio {np.min(dJ/ell):.6f}, below ell(1-1e-4) {(dJ < ell*(1-1e-4)).sum()}, below ell_min {(dJ < ra.ELL_LOWER).sum()}; horizon nonincreasing {sum(int(np.all(np.diff(L['N'])<=0)) for L in RL)}/{len(RL)}, min N {min(L['N'].min() for L in RL)}")
print(f"safety: min psi along closed loops RA {min(L['h_min'] for L in RL):.4f}, vanilla {min(v['h_min'] for v in VL):.4f}; entry median {np.median(te):.3f} max {te.max():.3f}; bound min {tb.min():.3f} median {np.median(tb):.1f} max {tb.max():.0f}; holds {(te <= tb+1e-9).sum()}/{len(RL)}")
for f in ("lemma2.pkl", "sweep.pkl", "timing.pkl"):
    if os.path.exists(f"{cf.DATA}/{f}"):
        D = pickle.load(open(f"{cf.DATA}/{f}", "rb")); print(f"\n## {f}")
        if f == "lemma2.pkl":
            print(f"eta {D['eta']:.4e} (lambda delta {D['lam_delta']:.4e}), M2 {D['M2_ref']:.3f}, T* {1e3*D['Tstar2']:.2f} ms; u^c on terminal set {np.round(D['umax_delta'],1)}, on X_RA {np.round(D['umax_0'],1)}")
            print("| T [ms] | arcs | smallest increase | increasing | end in terminal set | arc in X_RA | second-order bound |\n| --- | --- | --- | --- | --- | --- | --- |")
            for r in D["onestep"]: print(f"| {1e3*r['T']:g} | {r['n']} | {r['dmin']:+.2e} | {100*r['frac_inc']:.2f} % | {100*r['frac_inv']:.2f} % | {100*r['frac_in_xra']:.2f} % | {r['lb2']:+.2e} |")
            print("| T [ms] | from terminal set: reach | increasing | from X_RA: reach | increasing | median time |\n| --- | --- | --- | --- | --- | --- |")
            for r in D["closed_loop"]: print(f"| {1e3*r['T']:g} | {r['delta']['reached']}/{r['delta']['n']} | {r['delta']['monotone']}/{r['delta']['n']} | {r['zero']['reached']}/{r['zero']['n']} | {r['zero']['monotone']}/{r['zero']['n']} | {r['delta']['t_med']:.2f} s |")
        elif f == "sweep.pkl":
            print("| T [ms] | N | vanilla feasible (reached) | u^c held | RA-MPC | feasible at k=0 | later infeasible | decrease viol | min dPsi on terminal set | median entry | ms/solve |\n| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
            for r in D["rows"]: print(f"| {1e3*r['T']:g} | {r['N']} | {r['fv'].sum()} ({r['van_reached'].sum()}) | {r['zoh'].sum()} | {r['ra'].sum()} | {r['ra'].size - r['inf0']} | {r['later']} | {r['dec_viol']}/{r['steps']} | {r['eps']:+.1e} | {r['t_med']:.3f} s | {r['solve_ms']:.0f} |")
        else: print(D)
