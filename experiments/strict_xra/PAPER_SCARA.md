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
