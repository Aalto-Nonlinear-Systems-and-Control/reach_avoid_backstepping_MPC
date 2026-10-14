# Experiments of the journal version

Index of this folder; the repository README describes the method, the
requirements and the reproduction steps.

| Folder or file | Content |
| --- | --- |
| `solve_vanilla_k1_controller.m`, `solve_vanilla_k1_core.m` | SOS program for `k1` with `lambda` fixed (`LMIN`), `abs(k1_i) <= KMAX` on `{psi >= 0}` by the S-procedure, outputs solved in `y / SCALE`. Set through MATLAB globals by the design scripts; they shadow the copies of the letter version. |
| `design/` | Dubins certificate: `run_design.m` writes `k1/<weight>_l<lambda>.txt`; `cert.py` builds `Psi` and `u^c`; `ra_mpc_cert.py` is Algorithm 1 |
| `final/` | Example 1, Dubins car; notes in `PAPER_DUBINS.md`, design history in `REPORT.md` |
| `arm/` | workcell layout of Example 2 (`screen_racks.py` with `screen_front.py`, `screen_paths.py`, `screen_zigzag.py`, `cert_arm.py`) |
| `arm_top/` | Example 2, SCARA arm, eight kitting tasks; notes in `PAPER_ARM_TOP.md` |
| `scara/` | MuJoCo scene, stills and videos of the SCARA; notes in `PAPER_SCARA.md` |
| `franka/` | Panda-like arm in MuJoCo; notes in `PAPER_FRANKA.md` |
| `STUDY_14_DIRECT.md` | why the input-constrained program is solved in two steps; its scripts are on branch `tac-revision-archive` |

Earlier studies of this folder, including the first strict reach-avoid-set
report that this file used to hold, are on branch `tac-revision-archive`.
