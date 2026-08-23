#!/usr/bin/env python3
"""Plot isolated and cumulative Stokes contributions for m=0:35."""

from pathlib import Path
import os
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_untruncated_by_m")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import SymLogNorm

ROOT = Path(__file__).resolve().parents[1] / "truth_map_aerosols"
TABLE = ROOT / "untruncated_stokes_by_m.dat"
OUTPUT = ROOT / "untruncated_stokes_by_m.png"

data = np.loadtxt(TABLE, comments="#")
m = np.unique(data[:, 0]).astype(int)
vza = data[data[:, 0] == m[0], 1]
cube = data[:, 2:8].reshape(len(m), len(vza), 6)
isolated, cumulative = cube[:, :, :3], cube[:, :, 3:]

fig, axes = plt.subplots(3, 2, figsize=(14, 11), constrained_layout=True)
names = ("I", "Q", "U")
for s, name in enumerate(names):
    z = isolated[:, :, s].T
    vmax = np.max(np.abs(z))
    if vmax == 0:
        mesh = axes[s, 0].pcolormesh(m, vza, z, shading="nearest",
                                     cmap="RdBu_r", vmin=-1, vmax=1)
        axes[s, 0].text(.5, .5, "identically zero in the principal plane",
                        transform=axes[s, 0].transAxes, ha="center", va="center")
    else:
        nz = np.abs(z[np.nonzero(z)])
        linthresh = max(vmax * 1e-7,
                        np.percentile(nz, 2) if nz.size else 1e-30)
        mesh = axes[s, 0].pcolormesh(m, vza, z, shading="nearest",
            cmap="RdBu_r",
            norm=SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax))
    fig.colorbar(mesh, ax=axes[s, 0], label=f"isolated Δ{name}$_m$")
    axes[s, 0].set_title(f"Isolated Fourier contribution to {name}")
    axes[s, 0].set_ylabel("Signed VZA (deg)")

    # Plot every cumulative order, colored monotonically by m; emphasize a
    # compact set of milestones without hiding intermediate trajectories.
    colors = plt.cm.viridis(np.linspace(0, 1, len(m)))
    for im, order in enumerate(m):
        lw = 1.8 if order in (0, 1, 2, 4, 8, 12, 16, 24, 35) else .45
        alpha = .95 if lw > 1 else .28
        label = f"m≤{order}" if lw > 1 else None
        axes[s, 1].plot(vza, cumulative[im, :, s], color=colors[im],
                        lw=lw, alpha=alpha, label=label)
    axes[s, 1].set_title(f"Cumulative {name}")
    axes[s, 1].set_ylabel(name)
    axes[s, 1].grid(alpha=.25)
    axes[s, 1].legend(ncol=3, fontsize=7)

for ax in axes[-1, :]:
    ax.set_xlabel("Fourier order m" if ax is axes[-1, 0] else "Signed VZA (deg)")
fig.suptitle("Untruncated O₂ A-band Stokes contributions, m=0…35")
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)
