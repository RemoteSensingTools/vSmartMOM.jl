#!/usr/bin/env python3
"""Plot full-Mie aerosol phase-matrix diagnostics at 757.001655 nm."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "truth_map_aerosols" / "aerosol_phase_functions_757nm.dat"
OUTPUT = ROOT / "truth_map_aerosols" / "aerosol_phase_functions_757nm.png"
LABELS = {
    "sulfate": "Tropospheric sulfate",
    "organic_carbon": "Organic carbon",
    "utls_sulfate": "UTLS sulfate",
    "mixture": "Scattering-weighted mixture",
    "rayleigh": "Rayleigh reference",
}
COLORS = {
    "sulfate": "#d62728",
    "organic_carbon": "#1f77b4",
    "utls_sulfate": "#2ca02c",
    "mixture": "black",
    "rayleigh": "#777777",
}


def main():
    dtype = [("angle", float), ("species", "U20"), ("weight", float),
             ("f11", float), ("f12", float), ("dolp", float)]
    data = np.loadtxt(DATA, dtype=dtype, comments="#")
    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    for species in LABELS:
        d = data[data["species"] == species]
        width = 2.4 if species == "mixture" else 1.5
        linestyle = "--" if species == "rayleigh" else "-"
        axes[0].semilogy(d["angle"], d["f11"], color=COLORS[species],
                        linewidth=width, linestyle=linestyle, label=LABELS[species])
        axes[1].plot(d["angle"], d["f12"], color=COLORS[species],
                     linewidth=width, linestyle=linestyle)
        axes[2].plot(d["angle"], d["dolp"], color=COLORS[species],
                     linewidth=width, linestyle=linestyle)

    axes[0].set_ylabel(r"$F_{11}$")
    axes[1].set_ylabel(r"$F_{12}$")
    axes[1].set_yscale("symlog", linthresh=1e-3)
    axes[2].set_ylabel(r"$-F_{12}/F_{11}$")
    axes[2].set_xlabel("Scattering angle [deg]")
    axes[0].legend(frameon=False, ncol=2)
    axes[0].set_title(
        "Full-Mie aerosol phase functions at 757.001655 nm\n"
        "mixture weighted by scattering optical depth")
    for ax in axes:
        ax.set_xlim(0, 180)
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTPUT, dpi=220)
    print(OUTPUT)


if __name__ == "__main__":
    main()
