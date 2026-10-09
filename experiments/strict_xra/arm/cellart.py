"""Drawing helpers for the top-view kitting workcell: the bench fills the whole view, the robot is mounted on the bench."""
import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, Rectangle

BENCH, RACK, RACK_EDGE, FIX, FIX_EDGE = "#f4efe4", "#e4dccb", "#7d705a", "#9a9a9a", "#4d4d4d"
LINK, LINK_EDGE, JOINT, BASE, BASE_IN = "#b9bec4", "#5f6368", "#6b7075", "#3c4043", "#8a8f94"
OBJ_COLORS = ("#d62728", "#1f77b4", "#2ca02c", "#ff7f0e")


def bench(ax, view):
    (x0, x1), (y0, y1) = view; ax.add_patch(Rectangle((x0 - 1, y0 - 1), x1 - x0 + 2, y1 - y0 + 2, color=BENCH, zorder=0))


def racks(ax, X0, YS, depth, wall=0.5, opening=0.25, lw=1.4):
    for side in (-1, 1):
        xo, xb = side * (X0 - opening), side * (X0 + depth); ylo, yhi = YS[0] - wall, YS[-1] + wall
        ax.add_patch(Rectangle((min(xo, xb), ylo), abs(xb - xo), yhi - ylo, color=RACK, zorder=0.5))
        ax.plot([xb, xb], [ylo, yhi], color=RACK_EDGE, lw=2.4 * lw / 1.4, solid_capstyle="butt", zorder=1)
        for yy in [ylo] + [(YS[i] + YS[i + 1]) / 2 for i in range(len(YS) - 1)] + [yhi]: ax.plot([xo, xb], [yy, yy], color=RACK_EDGE, lw=lw, zorder=1)


def fixture(ax, w, y0, top, label=False):
    ax.add_patch(Rectangle((-w / 2, y0), w, top - y0, facecolor=FIX, edgecolor=FIX_EDGE, lw=1.0, zorder=1))
    if label: ax.text(0, (y0 + top) / 2, "fixture", ha="center", va="center", fontsize=8, color="w", zorder=2)


def objects(ax, X0, YS, which, size=0.4, alpha=1.0, outline=False, zorder=8):
    """which: list of (column side, height index, colour index)."""
    for side, j, k in which:
        p = np.array([side * X0, YS[j]]) - size / 2
        if outline: ax.add_patch(Rectangle(p, size, size, facecolor="none", edgecolor=OBJ_COLORS[k], lw=1.0, ls=":", zorder=zorder))
        else: ax.add_patch(Rectangle(p, size, size, facecolor=OBJ_COLORS[k], edgecolor="k", lw=0.6, alpha=alpha, zorder=zorder))


def _bar(ax, p, q, w, alpha, zorder):
    d = q - p; n = np.array([-d[1], d[0]]) / np.linalg.norm(d) * w / 2
    ax.add_patch(Polygon([p + n, q + n, q - n, p - n], closed=True, facecolor=LINK, edgecolor=LINK_EDGE, lw=0.9, alpha=alpha, zorder=zorder))


def robot(ax, q, l1=4.0, l2=4.0, w=0.34, alpha=0.75, zorder=6, gripper=True):
    """Two-link arm at joint angles q, base at the origin mounted on the bench."""
    q1, q2 = float(q[0]), float(q[1]); b = np.zeros(2); e = np.array([l1 * np.cos(q1), l1 * np.sin(q1)]); h = e + np.array([l2 * np.cos(q1 + q2), l2 * np.sin(q1 + q2)])
    ax.add_patch(Circle(b, 0.55, facecolor=BASE, edgecolor="k", lw=0.8, zorder=zorder - 0.2)); ax.add_patch(Circle(b, 0.30, facecolor=BASE_IN, edgecolor="k", lw=0.6, zorder=zorder + 0.3))
    _bar(ax, b, e, w, alpha, zorder); _bar(ax, e, h, w * 0.85, alpha, zorder + 0.1)
    ax.add_patch(Circle(e, 0.22, facecolor=JOINT, edgecolor="k", lw=0.7, alpha=min(1, alpha + 0.15), zorder=zorder + 0.2))
    if gripper:
        u = (h - e) / np.linalg.norm(h - e); n = np.array([-u[1], u[0]])
        for s in (1, -1): ax.plot(*np.array([h + s * 0.2 * n, h + s * 0.2 * n + 0.28 * u]).T, color=LINK_EDGE, lw=1.6, alpha=min(1, alpha + 0.15), solid_capstyle="round", zorder=zorder + 0.2)
        ax.plot(*np.array([h + 0.2 * n, h - 0.2 * n]).T, color=LINK_EDGE, lw=1.6, alpha=min(1, alpha + 0.15), zorder=zorder + 0.2)
    return h
