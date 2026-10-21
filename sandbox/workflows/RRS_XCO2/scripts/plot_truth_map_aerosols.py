#!/usr/bin/env python3
"""Plot continuous and exactly layer-integrated truth-map aerosol profiles."""

from pathlib import Path
from math import pi

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "truth_map_aerosols"
DATA = OUTDIR / "aerosol_vertical_profiles.dat"

SPECIES = (
    ("sulfate", "Tropospheric sulfate", "#d62728"),
    ("organic_carbon", "Organic carbon", "#1f77b4"),
    ("utls_sulfate", "UTLS sulfate", "#2ca02c"),
)


def load_rows():
    return np.genfromtxt(DATA, comments="#", dtype=None, encoding="utf-8")


def lognormal_pdf(z, median, sigma):
    result = np.zeros_like(z)
    positive = z > 0
    zp = z[positive]
    result[positive] = np.exp(-0.5 * ((np.log(zp) - np.log(median)) / sigma) ** 2) / (
        zp * sigma * np.sqrt(2 * pi)
    )
    return result


def continuous_plot(rows):
    fig, axes = plt.subplots(1, 3, figsize=(14, 6), sharey=True, constrained_layout=True)
    z = np.linspace(0.001, 25.0, 12000)
    sample = rows[rows["f0"] == 16]
    for ax, (key, title, color) in zip(axes, SPECIES):
        row = sample[sample["f6"] == key][0]
        tau, median, sigma, mode = row["f7"], row["f8"], row["f9"], row["f10"]
        extinction = tau * lognormal_pdf(z, median, sigma)
        ax.plot(extinction, z, color=color, lw=4, alpha=0.75)
        ax.axhline(mode, color=color, ls="--", lw=1.4, label=f"mode = {mode:g} km")
        ax.axhline(median, color=color, ls=":", lw=1.4, label=f"median = {median:.3g} km")
        ax.set_title(title)
        ax.set_xlabel(r"AOD density at 760 nm [km$^{-1}$]")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel("Altitude above surface [km]")
    axes[0].set_ylim(0, 25)
    fig.suptitle("Truth-map altitude-lognormal aerosol modes")
    fig.savefig(OUTDIR / "continuous_vertical_profiles.png", dpi=200)


def layer_plot(rows):
    fig, axes = plt.subplots(1, 2, figsize=(13, 10), sharey=True, constrained_layout=True)
    for ax, nlayer in zip(axes, (12, 16)):
        subset = rows[rows["f0"] == nlayer]
        geometry = subset[subset["f6"] == "sulfate"]
        top, bottom, dz = geometry["f2"], geometry["f3"], geometry["f5"]
        for boundary in np.r_[top[0], bottom]:
            ax.axhline(boundary, color="0.35", lw=0.65, alpha=0.7, zorder=0)
        for j, (key, title, color) in enumerate(SPECIES):
            values = subset[subset["f6"] == key]
            y = values["f3"] + values["f5"] * (5 - 2 * j) / 6
            ax.barh(y, values["f12"], height=values["f5"] / 3,
                    color=color, edgecolor=color, alpha=0.32,
                    label=f"{title}: layer AOD")
        ax.set_title(f"{nlayer}-layer model")
        ax.set_xlabel("Layer-integrated AOD at 760 nm")
        ax.set_ylim(0, 25)
        ax.grid(axis="x", alpha=0.2)
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    axes[0].set_ylabel("Altitude above surface [km]")
    fig.suptitle("Exact CDF integration of aerosol modes over model layers")
    fig.savefig(OUTDIR / "layer_integrated_profiles_12_vs_16.png", dpi=200)


def main():
    rows = load_rows()
    continuous_plot(rows)
    layer_plot(rows)
    print(OUTDIR / "continuous_vertical_profiles.png")
    print(OUTDIR / "layer_integrated_profiles_12_vs_16.png")


if __name__ == "__main__":
    main()
