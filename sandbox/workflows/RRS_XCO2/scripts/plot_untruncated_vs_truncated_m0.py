#!/usr/bin/env python3
"""Compare saved untruncated and truncated m=0 Stokes I/Q."""

from pathlib import Path
import os
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_compare_m0")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "truth_map_aerosols"
BY_M = ROOT / "untruncated_stokes_by_m.dat"
TRUNC_M0 = ROOT / "truncated_m0_stokes.dat"
OUTPUT = ROOT / "untruncated_vs_truncated_m0_IQ.png"
CAPS = (8, 12, 16, 24, 32)

un = np.loadtxt(BY_M, comments="#")
un = un[un[:, 0] == 0]
x = un[:, 1]

tr = np.loadtxt(TRUNC_M0, comments="#")

def signed_cap(cap, column):
    rows = tr[tr[:, 0] == cap]
    pos = rows[rows[:, 2] == 0]
    neg = rows[(rows[:, 2] == 180) & (rows[:, 1] > 0)]
    sx = np.r_[-neg[:, 1], pos[:, 1]]
    sy = np.r_[neg[:, column], pos[:, column]]
    order = np.argsort(sx)
    return sx[order], sy[order]

fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                         constrained_layout=True)
colors = plt.cm.viridis(np.linspace(.05, .85, len(CAPS)))
for ax, un_col, tr_col, name in zip(axes, (2, 3), (3, 4), ("I", "Q")):
    ax.plot(x, un[:, un_col], "ko-", lw=2.5, ms=4.5,
            label="untruncated m=0")
    for cap, color in zip(CAPS, colors):
        xt, yt = signed_cap(cap, tr_col)
        ax.plot(xt, yt, color=color, lw=1.5, label=f"l_trunc={cap}, m=0")
    ax.set_ylabel(name)
    ax.set_xlim(-75, 75)
    ax.grid(alpha=.3)
    ax.legend(ncol=2, fontsize=8)
axes[-1].set_xlabel(
    "Signed VZA (deg): negative = RelAz 180°, positive = RelAz 0°")
fig.suptitle("O₂ A-band m=0: untruncated versus δ-BGE truncated")
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)
