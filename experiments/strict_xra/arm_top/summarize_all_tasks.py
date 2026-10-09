"""Summary of all eight kitting tasks of the top-view arm (square layout, common torque bounds of red_fwd_sq), with the
same definitions as summarize_top.py: design and scale, comparison on the 150-state common initial set of each task,
Theorem 2 / Corollary 1 checks, Lemma 2, times and solve times. Reads data/<task>/{main,lemma2,timing}.pkl and
data/scale_<task>.json. Writes data/all_tasks_summary.json and data/all_tasks_summary.md."""
import os, json, pickle, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); os.environ.setdefault("TASK", "red_fwd_sq")
import sys; sys.path.insert(0, HERE)
from ra_mpc_top import ELL_LOWER
TASKS = [f"{n}_{d}_sq" for d in ("fwd", "rev") for n in ("red", "blue", "green", "orange")]; T = 0.01; rows = []
for task in TASKS:
    P = os.path.join(HERE, "data", task)
    if not os.path.exists(f"{P}/main.pkl"): print("missing", task); continue
    M = pickle.load(open(f"{P}/main.pkl", "rb")); logs, van, cfg = M["logs"], M["van"], M["cfg"]; U_MAX = np.asarray(cfg["U_MAX"], float)
    SC = json.load(open(os.path.join(HERE, "data", f"scale_{task}.json"))); s = SC["s"]; t1 = np.asarray(SC.get("tau1_terminal", np.asarray(SC["tau_terminal"]) / s**2), float)
    s_eps = SC.get("s_eps", 0.8 * SC["T_eps1_ms"] / 1e3 / T); s_bnd = SC.get("s_bnd", float(np.sqrt((U_MAX / (1.25 * t1)).min())))
    ura = np.array([u["status"] == "reached" for u in M["ura"]]); fv = np.array([len(v["J"]) > 0 for v in van]); vr = np.array([v["status"] == "reached" for v in van]); rr = np.array([L["status"] == "reached" for L in logs])
    viol = (np.array([u["upeak"] for u in M["ura"]]) > U_MAX).any(1); RL = [L for L in logs if L["status"] == "reached"]
    tu = np.array([u["t"] for u in M["ura"]]); tr = np.array([L["check"]["t_entry"] if L["status"] == "reached" else np.nan for L in logs]); both = ura & rr
    cu = np.array([u["cost"] for u in M["ura"]]); cr = np.array([T * np.sum(L["ell"]) if L["status"] == "reached" else np.nan for L in logs])
    ell = np.concatenate([L["ell"][:-1] for L in RL if len(L["J"]) > 1]); dJ = np.concatenate([-np.diff(L["J"]) for L in RL if len(L["J"]) > 1])
    te = np.array([L["check"]["t_entry"] for L in RL]); tb = np.array([L["check"]["t_bound"] for L in RL])
    r = dict(task=task, lam=float(open(os.path.join(HERE, "data", f"lam_{task}.txt")).read()) if os.path.exists(os.path.join(HERE, "data", f"lam_{task}.txt")) else 0.1,
             s=s, s_eps=s_eps, s_bnd=s_bnd, binding="bounds" if s_bnd < s_eps - 1e-9 else "period", T_eps1_ms=SC["T_eps1_ms"], T_eps_s_ms=SC["T_eps_s_ms"], Tstar2_ms=SC["Tstar2_ms"],
             U_MAX=U_MAX.tolist(), own_bounds=(np.ceil(1.25 * s_eps**2 * t1 * 2) / 2).tolist(), ura_scale=SC["URA_SCALE"], n=int(rr.size),
             ura=int(ura.sum()), ura_exceed=int((viol & ura).sum()), van_feasible0=int(fv.sum()), van=int(vr.sum()), ra=int(rr.sum()), ra_infeasible0=sum(len(L["J"]) == 0 for L in logs),
             ra_later_infeasible=sum(L["status"] == "infeasible" and len(L["J"]) > 0 for L in logs), ra_not_van=int((rr & ~vr).sum()), van_not_ra=int((vr & ~rr).sum()),
             ra_not_ura=int((rr & ~ura).sum()), ura_not_ra=int((ura & ~rr).sum()), solves=int(sum(len(L["J"]) for L in RL)), decrease_viol=int((dJ < ell * (1 - 1e-4)).sum()),
             below_ell_min=int((dJ < ELL_LOWER).sum()), horizon_nonincreasing=int(sum(int(np.all(np.diff(L["N"]) <= 0)) for L in RL)), uc_in_bounds=bool(all(np.all(L["uc_in_bounds"]) for L in RL)),
             min_psi_ra=float(min(L["h_min"] for L in RL)), entry_bound_holds=int((te <= tb + 1e-9).sum()), t_med_ura=float(np.median(tu[both])), t_med_ra=float(np.nanmedian(tr[both])),
             cost_med_ura=float(np.median(cu[both])), cost_med_ra=float(np.nanmedian(cr[both])))
    if os.path.exists(f"{P}/lemma2.pkl"):
        L2 = pickle.load(open(f"{P}/lemma2.pkl", "rb")); inc = [o["T"] for o in L2["onestep"] if o["frac_inc"] >= 1.0]
        r.update(lemma2_Tstar2_ms=1e3 * L2["Tstar2"], lemma2_largest_T_all_increasing_ms=1e3 * max(inc) if inc else 0.0,
                 lemma2_cl_10ms={k: f"{c[k]['reached']}/{c[k]['n']}" for c in L2["closed_loop"] if abs(c["T"] - T) < 1e-9 for k in ("delta", "zero")})
    if os.path.exists(f"{P}/timing.pkl"):
        TM = pickle.load(open(f"{P}/timing.pkl", "rb")); r.update(solve_ms_ra_median=TM["ra"]["median"], solve_ms_van_median=TM["vanilla"]["median"])
    rows.append(r)
json.dump(rows, open(os.path.join(HERE, "data", "all_tasks_summary.json"), "w"), indent=1)
L = ["# All eight kitting tasks (square layout, common torque bounds (28.5, 8.5) N m, T = 10 ms, N = 60, 150 initial states each)", "",
     "## Design and scale", "", "| task | lambda | s_eps | s_bnd | s | binding | T_eps(1) [ms] | T_eps(s) [ms] | second-order T* [ms] | own bounds would be | URA scale |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
L += [f"| {r['task']} | {r['lam']:g} | {r['s_eps']:.3f} | {r['s_bnd']:.3f} | {r['s']:.3f} | {r['binding']} | {r['T_eps1_ms']:.1f} | {r['T_eps_s_ms']:.1f} | {r['Tstar2_ms']:.2f} | {tuple(r['own_bounds'])} | {r['ura_scale']:.3f} |" for r in rows]
L += ["", "## Comparison", "", "| task | URA (exceed bounds) | vanilla feasible k=0 / reached | RA-MPC (infeasible k=0, later) | RA not van | van not RA | RA not URA | URA not RA | median time URA / RA [s] | median cost URA / RA |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
L += [f"| {r['task']} | {r['ura']} ({r['ura_exceed']}) | {r['van_feasible0']} / {r['van']} | {r['ra']} ({r['ra_infeasible0']}, {r['ra_later_infeasible']}) | {r['ra_not_van']} | {r['van_not_ra']} | {r['ra_not_ura']} | {r['ura_not_ra']} | {r['t_med_ura']:.2f} / {r['t_med_ra']:.2f} | {r['cost_med_ura']:.1f} / {r['cost_med_ra']:.1f} |" for r in rows]
L += ["", "## Theorem 2 / Corollary 1 and Lemma 2", "", "| task | solves | decrease violations | below ell_min | horizon nonincreasing | u^c in bounds | min psi RA | entry bound holds | Lemma 2: all increasing up to [ms] | T* [ms] | ZOH closed loop at 10 ms (terminal set, X_RA) | solve ms RA / vanilla |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
L += [f"| {r['task']} | {r['solves']} | {r['decrease_viol']} | {r['below_ell_min']} | {r['horizon_nonincreasing']}/{r['ra']} | {r['uc_in_bounds']} | {r['min_psi_ra']:.4f} | {r['entry_bound_holds']}/{r['ra']} | {r.get('lemma2_largest_T_all_increasing_ms', float('nan')):g} | {r.get('lemma2_Tstar2_ms', float('nan')):.2f} | {', '.join(r.get('lemma2_cl_10ms', {}).values())} | {r.get('solve_ms_ra_median', float('nan')):.0f} / {r.get('solve_ms_van_median', float('nan')):.0f} |" for r in rows]
open(os.path.join(HERE, "data", "all_tasks_summary.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
