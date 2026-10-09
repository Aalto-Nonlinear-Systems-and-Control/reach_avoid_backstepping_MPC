# Two-link example as a SCARA in MuJoCo: video of the red task

Setting of 2026-10-09 (scripts in `scara/`, environment `rab_franka`). The horizontal two-link arm of `arm_top/`
(task `red_fwd_sq`, lambda 0.1, s = 2.39, T = 10 ms, N = 60) is the planar part of a SCARA: both hinges turn about vertical
axes, so gravity is carried by the bearings and g = 0 holds exactly. Nothing in `arm_top/` is changed.

## 1. Scaling to a SCARA

| item | paper model (task units) | SCARA (SI) |
|---|---|---|
| link lengths, centres of mass | 4, 2 | 0.4 m, 0.2 m (ALPHA = 0.1) |
| link masses, inertias | 1, 0.02 | 10 kg, 0.002 kg m^2 (MSCALE = 10) |
| mass matrix, torques | M(q), u | 0.1 M(q), 0.1 u (C_TAU = MSCALE ALPHA^2) |
| torque bounds | (28.5, 8.5) | (2.85, 0.85) N m |
| joint trajectories, time | q(t) | the same q(t) |

The certificate, the terminal controller and the MPC map one to one (psi(y / ALPHA), ALPHA k1(y / ALPHA), output weight / ALPHA^2, input weight / C_TAU^2), so the closed loops of the paper are SCARA closed loops.
Task coordinates map to the world as X = ALPHA y2, Y = -ALPHA y1, the layout of the Franka scene.
Checks: MuJoCo's mass matrix equals C_TAU M(q) to 3e-15 and the tip equals the model's kinematics to 4e-16 m (200 random q).

## 2. Demo (`scara_demo_traj.py`, `data/demo.npz`)

Pick 2.3 s: the arm rests above the red cube (u = 0) while the Z axis lowers the gripper from 0.30 m to the cube, the fingers close and the Z axis lifts it back.
Transport: RA-MPC from rest, reached X_T in 1.40 s (140 steps, 202 ms per solve), min psi along the path 0.143.
Place 2.9 s: the tool point is centred over the slot (minimum-jerk reference, clipped computed torque, not certified), then the arm holds while the Z axis lowers the cube, the fingers open and the Z axis retracts.
The planar part runs on the paper's plant; the same inputs times C_TAU re-simulated on `model/scara_arm.xml` differ by at most 6.7e-16. Peak torques equal the bounds (2.85, 0.85) N m, never above.
Clearances: fingers to the obstacle 52 mm, carried cube to the other cubes 60 mm, fingers to the tabletop 8 mm; the links pass 0.2 m above the 0.30 m obstacle.
The Z axis, the wrist and the fingers are scripted mocap bodies (they do not enter the planar dynamics; the wrist keeps the gripper's world orientation).

## 3. Video and replay window

`figs/scara_red_task.mp4` (`render_scara.py`, 31 s): angled front view (pick and place at 0.5x, transport at 0.25x), then the fixed top camera (pick and place at 1x, transport at 0.25x) with the arm drawn as a translucent ghost, as in the paper figures, so the corridor, the trail and the cubes stay visible.
Replay window `./view_scara.sh [speed]` (default 0.5).

Reproduction: `python build_scene_scara.py`, `python scara_demo_traj.py` (about 30 s), `python render_scara.py [--stills]` (about 2 min).

## 4. Full kitting cycle (`cycle_demo_traj.py`, `data/cycle.npz`, `figs/scara_cycle.mp4`)

Eight tasks in a row: red, blue, green and orange cubes to their slots (crossing assignment), then back to the rack, 69.75 s of simulated time.
Every task uses its own certificate (k1 at lambda 0.1) and scale under the **common torque bounds** (28.5, 8.5), i.e. (2.85, 0.85) N m on the SCARA; see `PAPER_ARM_TOP.md` Section 11 for why one set of bounds works for all tasks.
Per task: pick (Z axis from 0.40 m to the cube and up to 0.30 m, the arm holding its position), RA-MPC from the actual state (certified, `ra_run_task.py` in a subprocess, logs in `data/cycle_ra/<task>.pkl`), place (centre over the slot, Z axis down, release, up to 0.40 m), transit to the next cube.
Transits are point-to-point joint-space minimum-jerk moves at 0.40 m (the fingers clear the 0.30 m obstacle by 88 mm; joint 2 is interpolated, so the tool stays between the start and end radii, away from the column), with the shortest duration whose feedforward torque stays within 90 % of the bounds: no clipping in any transit.
A first version with straight task-space transits failed on the move from the left rack to the right slot across the front of the base (272 clipped steps, 41 cm short), which put the next RA-MPC start outside X_S.

| task | scale s | RA-MPC time [s] | steps | solve median [ms] | min psi along the transport | centring [s] | transit [s] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| red_fwd_sq | 2.386 | 1.40 | 140 | 371 | 0.1429 | 0.63 | 2.25 |
| blue_fwd_sq | 1.747 | 1.62 | 162 | 398 | 0.0059 | 0.63 | 3.00 |
| green_fwd_sq | 1.623 | 1.40 | 140 | 359 | 0.3914 | 0.62 | 2.25 |
| orange_fwd_sq | 1.690 | 1.65 | 165 | 389 | 0.0058 | 0.63 | 1.00 |
| red_rev_sq | 1.175 | 1.32 | 132 | 358 | 0.4481 | 0.62 | 2.75 |
| blue_rev_sq | 1.767 | 1.61 | 161 | 372 | 0.4023 | 0.62 | 2.25 |
| green_rev_sq | 1.499 | 1.35 | 135 | 271 | 0.0093 | 0.62 | 2.50 |
| orange_rev_sq | 1.653 | 1.61 | 161 | 252 | 0.0050 | 0.62 | (last) |

psi is normalised to maximum 1; the smallest values equal the MPC's safety back-off 0.005 (the corridors are narrow). Solve times are under the parallel load of the main experiments.
Model check: every held-input period is re-simulated on `model/scara_arm.xml` from the paper plant's state at its start; the largest one-period difference is 2.2e-16. (An open-loop replay of the whole 70 s amplifies round-off by about 40 per task and reached 5e-3 rad, since the recorded inputs carry no feedback; the single red task is short enough for an open-loop replay, 6.7e-16.)

Video `figs/scara_cycle.mp4` (`render_cycle.py`, 84 s): angled front view with the top view as an inset at the top right (arm drawn as a ghost); the current task's corridor edge, target ring and transport trail; pick, place and transit at 1x, transports at 0.5x. Scene `model/scene_scara_cycle.xml` (all cubes mocap, slots dashed in the colour of their cube, homes in a lighter shade). Replay window `./view_cycle.sh [speed]` (default 1.0).

Reproduction: `python build_scene_scara.py`, `python cycle_demo_traj.py` (about 7 min), `python render_cycle.py [--stills]` (about 9 min).
