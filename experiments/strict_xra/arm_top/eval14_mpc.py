"""RA-MPC with a terminal design from study14.py on the task's common 150-state initial set (data/<task>/main.pkl), all
else as in run_main_top.py (T, N, bounds, cost, solver). Compares with the RA-MPC runs of the scaled design (12)
stored in main.pkl. Usage: TASK=blue_fwd_sq LAM=0.1 python eval14_mpc.py <s> <kfac>
Writes data/study14/mpc_<task>_s<s>_k<kfac>.pkl."""
import os, sys, json, time, pickle, numpy as np
from collections import Counter
from joblib import Parallel, delayed
import config_top as cf
from config_top import ra
s, kfac = float(sys.argv[1]), float(sys.argv[2]); TAG = f"{cf.TASK}_s{s:g}_k{kfac:g}"; D = json.load(open(f"{cf.HERE}/data/study14/{TAG}.json"))
KF = os.path.join(cf.HERE, "k1", D["k1_file"]); MU = D["mu"]


def one(x0):
    L = ra.run_one(np.asarray(x0, float), KF, MU, 1.0, cf.T, cf.N, cf.DFRAC, None, cf.T_MAX_MPC, "cert", cf.U_MAX, cf.SAFE_MARGIN, max(1, int(round(cf.T / 0.002))))
    L["h_min"], L["entered"] = cf.safety_until_target(L["x"], L["u"], cf.T) if len(L["u"]) else (np.nan, False); return L


if __name__ == "__main__":
    t0 = time.time(); M = pickle.load(open(f"{cf.DATA}/main.pkl", "rb")); X0 = M["X0"]; n = X0.shape[1]
    logs = Parallel(n_jobs=-1, prefer="processes")(delayed(one)(X0[:, i]) for i in range(n))
    rb = np.array([L["status"] == "reached" for L in logs]); ra_ = np.array([L["status"] == "reached" for L in M["logs"]])
    RL = [L for L in logs if L["status"] == "reached"]; RA = [L for L in M["logs"] if L["status"] == "reached"]
    ell = np.concatenate([L["ell"][:-1] for L in RL if len(L["J"]) > 1]); dJ = np.concatenate([-np.diff(L["J"]) for L in RL if len(L["J"]) > 1])
    both = rb & ra_; tB = np.array([L["check"]["t_entry"] if L["status"] == "reached" else np.nan for L in logs]); tA = np.array([L["check"]["t_entry"] if L["status"] == "reached" else np.nan for L in M["logs"]])
    cB = np.array([cf.T * np.sum(L["ell"]) if L["status"] == "reached" else np.nan for L in logs]); cA = np.array([cf.T * np.sum(L["ell"]) if L["status"] == "reached" else np.nan for L in M["logs"]])
    out = dict(tag=TAG, s=s, kfac=kfac, reached=int(rb.sum()), statuses=dict(Counter(L["status"] for L in logs)), infeasible0=sum(len(L["J"]) == 0 for L in logs),
               later_infeasible=sum(L["status"] == "infeasible" and len(L["J"]) > 0 for L in logs), A_reached=int(ra_.sum()), B_not_A=int((rb & ~ra_).sum()), A_not_B=int((ra_ & ~rb).sum()),
               decrease_viol=int((dJ < ell * (1 - 1e-4)).sum()), horizon_nonincreasing=int(sum(int(np.all(np.diff(L["N"]) <= 0)) for L in RL)), uc_in_bounds=bool(all(np.all(L["uc_in_bounds"]) for L in RL)),
               min_psi=float(min(L["h_min"] for L in RL)) if RL else np.nan, t_med_B=float(np.nanmedian(tB[both])), t_med_A=float(np.nanmedian(tA[both])), cost_med_B=float(np.nanmedian(cB[both])),
               cost_med_A=float(np.nanmedian(cA[both])), solve_ms_B=float(np.mean(np.concatenate([L["solve_ms"] for L in RL]))), solve_ms_A=float(np.mean(np.concatenate([L["solve_ms"] for L in RA]))), sec=time.time() - t0)
    print(json.dumps(out)); pickle.dump(dict(summary=out, logs=logs), open(f"{cf.HERE}/data/study14/mpc_{TAG}.pkl", "wb"))
