"""Offline computation time of Table I: every step before the closed loop, each timed alone in one process (median of REPS
repetitions), with nothing else running. The time of the SOS design is measured by run_design.m. Writes data/timing_offline.json."""
import contextlib, io, json, os, runpy, time, numpy as np
import config as cf
from config import ra, CertArm

REPS = 3


def clock(fn):
    ts = []
    for _ in range(REPS):
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


script = lambda name: (lambda: runpy.run_path(os.path.join(cf.HERE, name)))
sos = float(open(f"{cf.DATA}/sos_time.txt").read())
steps = dict(sos_design=sos,
             controller=clock(lambda: CertArm(cf.K1_FILE, mu=cf.MU_BASE).scaled(cf.SCALE)),    # Psi and u^c from k1
             time_scale=clock(script("solve_scale.py")),                                        # s, input bounds and sampling period
             terminal_weight=clock(lambda: ra.ArmMPC(cf.C, dt=cf.T, N=cf.N, dfrac=cf.DFRAC, safe_margin=cf.SAFE_MARGIN, u_max=cf.U_MAX, nsub=cf.NSUB, terminal="cert")))   # eps, omega, model
out = dict(steps=steps, URA=steps["sos_design"] + steps["controller"], RA_MPC=sum(steps.values()), reps=REPS); print(out)
json.dump(out, open(f"{cf.DATA}/timing_offline.json", "w"), indent=1)
