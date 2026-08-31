#!/usr/bin/env python3
"""Plot a matched aerosol/no-aerosol truth-map pair in a 3x3 component grid."""

from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TRUTH = ROOT / "truth_map"
NO_AEROSOL_STATE = 1
AEROSOL_STATE = 9
OUTPUT = TRUTH / "aerosol_vs_noaerosol_components_state001_state009.png"


def read_scene(index):
    with nc.Dataset(TRUTH / f"hiressim_{index:03d}.nc") as ds:
        return {name: np.asarray(var[:]) for name, var in ds.variables.items()}


def read_wavelengths():
    with nc.Dataset(TRUTH / "sim_wavelength.nc") as ds:
        return {
            "o2a": np.asarray(ds.variables["o2a_wavelength"][:]),
            "weak_co2": np.asarray(ds.variables["weak_co2_wavelength"][:]),
            "strong_co2": np.asarray(ds.variables["strong_co2_wavelength"][:]),
        }


def stokes_i(scene, variable, n_wavelength):
    """Return Stokes I independent of NetCDF C/Fortran dimension ordering."""
    values = scene[variable]
    if values.shape == (n_wavelength, 3):
        return values[:, 0]
    if values.shape == (3, n_wavelength):
        return values[0, :]
    raise ValueError(f"Unexpected {variable} shape: {values.shape}")


def main():
    clear = read_scene(NO_AEROSOL_STATE)
    aerosol = read_scene(AEROSOL_STATE)
    wavelength = read_wavelengths()

    colors = {"No aerosol": "#2166ac", "Aerosol": "#b2182b"}
    scenes = {"No aerosol": clear, "Aerosol": aerosol}
    bands = (
        ("O$_2$ A-band", "o2a"),
        ("Weak CO$_2$", "weak_co2"),
        ("Strong CO$_2$", "strong_co2"),
    )
    columns = ("Cabannes / elastic", "RRS", "Cabannes + RRS − Rayleigh")

    fig, axes = plt.subplots(3, 3, figsize=(16, 11), constrained_layout=True)
    for row, (band_title, band) in enumerate(bands):
        wl = wavelength[band]
        order = np.argsort(wl)
        for col, title in enumerate(columns):
            ax = axes[row, col]
            if row == 0:
                for label, scene in scenes.items():
                    ray = stokes_i(scene, "radiance_rayleigh_o2a", len(wl))
                    cab = stokes_i(scene, "radiance_cabannes_o2a", len(wl))
                    rrs = stokes_i(scene, "radiance_rrs_o2a", len(wl))
                    values = (cab, rrs, cab + rrs - ray)[col]
                    ax.plot(wl[order], values[order], color=colors[label],
                            lw=0.9, label=label)
            elif col == 0:
                variable = f"radiance_rayleigh_{band}"
                for label, scene in scenes.items():
                    values = stokes_i(scene, variable, len(wl))
                    ax.plot(wl[order], values[order], color=colors[label],
                            lw=0.9, label=label)
                ax.text(0.02, 0.04, "noRS elastic output\n(RRS negligible/not simulated)",
                        transform=ax.transAxes, fontsize=8, va="bottom")
            else:
                ax.axhline(0.0, color="0.35", lw=0.8)
                ax.set_ylim(-1.0e-6, 1.0e-6)
                ax.text(0.5, 0.5, "RRS not simulated\n(defined as zero)",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=10, color="0.35")

            if row == 0:
                ax.set_title(title)
            if col == 0:
                ax.set_ylabel(f"{band_title}\nTOA Stokes I\n[mW m$^{{-2}}$ sr$^{{-1}}$ (cm$^{{-1}}$)$^{{-1}}$]")
            if row == 2:
                ax.set_xlabel("Wavelength [nm]")
            ax.grid(alpha=0.2)

    axes[0, 0].legend(loc="best", frameon=False)
    fig.suptitle(
        "Matched urban scenes: state 001 (no aerosol) vs state 009 "
        "(AOD$_{760}$ = 0.28)\nSIF off, XCO$_2$ = 380 ppm, SZA = 30°, nadir",
        fontsize=14,
    )
    fig.savefig(OUTPUT, dpi=180)
    print(OUTPUT)


if __name__ == "__main__":
    main()
