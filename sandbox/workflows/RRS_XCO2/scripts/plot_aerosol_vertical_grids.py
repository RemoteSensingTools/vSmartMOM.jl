#!/usr/bin/env python3
"""Compare intended altitude lognormals with exact 12/16-layer AOD allocation."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "truth_map" / "aerosol_vertical_grids.dat"
OUTPUT = ROOT / "truth_map" / "aerosol_vertical_distribution_12_vs_16_layers.png"

SPECIES = (
    ("Sulfate", 1.0, 0.38, "#d62728"),
    ("Organic carbon", 1.2, 0.35, "#1f77b4"),
    ("Stratospheric sulfate", 12.0, 0.16, "#2ca02c"),
)


def lognormal_pdf(z, median, sigma):
    out = np.zeros_like(z)
    positive = z > 0
    zp = z[positive]
    out[positive] = np.exp(-0.5 * ((np.log(zp) - np.log(median)) / sigma) ** 2) / (
        zp * sigma * np.sqrt(2 * np.pi)
    )
    return out


def main():
    data = np.loadtxt(DATA)
    fig, axes = plt.subplots(1, 2, figsize=(13, 10), sharey=True,
                             constrained_layout=True)
    z = np.linspace(0.001, 62.8, 12000)

    for ax, nlayer in zip(axes, (12, 16)):
        rows = data[data[:, 0] == nlayer]
        top, bottom, center, dz = rows[:, 2], rows[:, 3], rows[:, 4], rows[:, 5]

        for boundary in np.r_[top[0], bottom]:
            ax.axhline(boundary, color="0.35", lw=0.65, alpha=0.7, zorder=0)

        # Each horizontal bar occupies one third of its layer height. The
        # offsets order sulfate, organic carbon, and stratospheric sulfate
        # from top to bottom without overlap.
        for j, (name, median, sigma, color) in enumerate(SPECIES):
            y = bottom + dz * (5 - 2 * j) / 6
            ax.barh(y, rows[:, 6 + j], height=dz / 3, left=0,
                    color=color, alpha=0.25, edgecolor=color, linewidth=0.7,
                    label=f"{name}: layer AOD" if nlayer == 12 else None)

        ax.set_xlabel("AOD at 760 nm in layer")
        ax.set_title(f"{nlayer} layers")
        ax.set_ylim(0, max(top) + 0.5)
        ax.grid(axis="x", alpha=0.15)

        curve_ax = ax.twiny()
        for name, median, sigma, color in SPECIES:
            curve_ax.plot(lognormal_pdf(z, median, sigma), z, color=color,
                          lw=4, alpha=0.62,
                          label=f"{name}: intended altitude PDF" if nlayer == 16 else None)
        curve_ax.set_xlabel("Assigned altitude-lognormal density [km$^{-1}$]")
        curve_ax.set_xlim(left=0)

        if nlayer == 12:
            ax.legend(loc="upper right", fontsize=8, frameon=False)
        else:
            curve_ax.legend(loc="center right", fontsize=8, frameon=False)

    axes[0].set_ylabel("Altitude above surface, z [km]")
    fig.suptitle(
        "Aerosol vertical profiles: intended altitude lognormals and exact model layer AOD\n"
        "Horizontal lines are layer boundaries; bars occupy one third of each layer",
        fontsize=14,
    )
    fig.savefig(OUTPUT, dpi=190)
    print(OUTPUT)


if __name__ == "__main__":
    main()
