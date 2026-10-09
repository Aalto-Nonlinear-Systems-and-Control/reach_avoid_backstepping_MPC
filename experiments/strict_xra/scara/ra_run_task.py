"""RA-MPC of one kitting task of ../arm_top from a given state (its own certificate, scale and the common bounds), for the
SCARA cycle demo; run as a subprocess because the task is fixed when config_top is imported.
Usage: python ra_run_task.py <task> <x0.npy> <out.pkl>"""
import os, sys, pickle, numpy as np
task, x0_file, out = sys.argv[1:4]; HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["TASK"] = task; os.environ["LAM"] = open(os.path.join(HERE, "..", "arm_top", "data", f"lam_{task}.txt")).read().strip()
sys.path.insert(0, os.path.join(HERE, "..", "arm_top"))
import config_top as cf
x0 = np.load(x0_file); L = cf.mpc_run(x0, "cert")
pickle.dump(dict(task=task, lam=os.environ["LAM"], status=L["status"], t=L["t"], x=L["x"], u=L["u"], J=L["J"], N=L["N"], solve_ms=L["solve_ms"], scale=cf.SCALE, U_MAX=cf.U_MAX, check=L["check"]), open(out, "wb"))
print(f"{task}: {L['status']} in {L['t'][-1]:.2f} s, {len(L['u'])} steps, solve median {np.median(L['solve_ms']):.0f} ms, scale {cf.SCALE:.3f}, bounds {cf.U_MAX}")
