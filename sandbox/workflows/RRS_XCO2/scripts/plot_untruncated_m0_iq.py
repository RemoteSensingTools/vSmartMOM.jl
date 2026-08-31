#!/usr/bin/env python3
"""Overplot isolated untruncated m=0:3 Stokes I/Q."""

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_untruncated_m0")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "truth_map_aerosols"
TABLE = ROOT / "untruncated_stokes_by_m.dat"
OUTPUT = ROOT / "untruncated_m0_m1_m2_m3_IQ_vs_vza.png"

data = np.loadtxt(TABLE, comments="#")
m0 = data[data[:, 0] == 0]
m1 = data[data[:, 0] == 1]
m2 = data[data[:, 0] == 2]
m3 = data[data[:, 0] == 3]
vza = m0[:, 1]

fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                         constrained_layout=True)
for ax, column, label in zip(axes, (2, 3), ("I", "Q")):
    ax.plot(vza, m0[:, column], "o-", color="tab:blue", lw=1.6, ms=4,
            label="m=0")
    ax.plot(vza, m1[:, column], "s-", color="tab:orange", lw=1.6, ms=3.5,
            label="m=1")
    ax.plot(vza, m2[:, column], "^-", color="tab:green", lw=1.6, ms=4,
            label="m=2")
    ax.plot(vza, m3[:, column], "d-", color="tab:red", lw=1.6, ms=3.5,
            label="m=3")
    ax.set_ylabel(label)
    ax.grid(alpha=0.3)
    ax.set_xlim(-75, 75)
    ax.legend()
axes[-1].set_xlabel("Signed VZA (deg): negative = RelAz 180°, positive = RelAz 0°")
fig.suptitle("Untruncated O₂ A-band: isolated Fourier orders m=0…3")
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)
