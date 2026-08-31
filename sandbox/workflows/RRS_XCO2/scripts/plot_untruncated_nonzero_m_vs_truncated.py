#!/usr/bin/env python3
"""Compare saved untruncated m>=1 sums with saved truncated solutions."""

from pathlib import Path
import os
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_nonzero_m")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "truth_map_aerosols"
BY_M = ROOT / "untruncated_stokes_by_m.dat"
TRUNC = ROOT / "truncation_iqU.dat"
OUTPUT = ROOT / "untruncated_nonzero_m_vs_truncated_IQ.png"
CAPS = (8, 12, 16, 24, 32)

moments = np.loadtxt(BY_M, comments="#")
m0 = moments[moments[:, 0] == 0]
x = m0[:, 1]
m0_I, m0_Q = m0[:, 2], m0[:, 3]

dtype = [("vza", float), ("vaz", float), ("case", "U12"),
         ("I", float), ("Q", float), ("U", float)]
saved = np.loadtxt(TRUNC, dtype=dtype, comments="#")

def signed(case, component):
    pos = saved[(saved["case"] == case) & (saved["vaz"] == 0)]
    neg = saved[(saved["case"] == case) & (saved["vaz"] == 180)]
    neg = neg[neg["vza"] > 0]
    sx = np.r_[-neg["vza"], pos["vza"]]
    sy = np.r_[neg[component], pos[component]]
    order = np.argsort(sx)
    return sx[order], sy[order]

xfull, full_I = signed("untruncated", "I")
_, full_Q = signed("untruncated", "Q")
if not np.array_equal(x, xfull):
    raise ValueError("saved m=0 and full-solution VZA grids differ")
nonzero_I = full_I - m0_I
nonzero_Q = full_Q - m0_Q

fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                         constrained_layout=True)
zero = np.flatnonzero(x == 0)[0]
axes[0].plot(x, nonzero_I - nonzero_I[zero], "ko-", lw=2.2, ms=4,
             label="untruncated sum m=1…220")
axes[1].plot(x, nonzero_Q, "ko-", lw=2.2, ms=4,
             label="untruncated sum m=1…220")

colors = plt.cm.viridis(np.linspace(.05, .85, len(CAPS)))
for cap, color in zip(CAPS, colors):
    xt, It = signed(str(cap), "I")
    _, Qt = signed(str(cap), "Q")
    axes[0].plot(xt, It - It[zero], color=color, lw=1.4,
                 label=f"full truncated l={cap}")
    axes[1].plot(xt, Qt, color=color, lw=1.4,
                 label=f"full truncated l={cap}")

axes[0].set_ylabel("I − I(nadir)")
axes[0].set_title("Angular variation in I (nadir baseline removed)")
axes[1].set_ylabel("Q")
axes[1].set_title("Q (nadir is zero by symmetry)")
axes[1].set_xlabel("Signed VZA (deg): negative = RelAz 180°, positive = RelAz 0°")
for ax in axes:
    ax.set_xlim(-75, 75)
    ax.grid(alpha=.3)
    ax.legend(ncol=2, fontsize=8)
fig.suptitle("Saved untruncated nonzero Fourier sum versus saved truncated solutions")
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)
