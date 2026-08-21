#!/usr/bin/env python3
"""Plot ECOSTRESS mixtures and their bandwise P0--P2 representations."""
from pathlib import Path
import csv
import os

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[1]
out = root / "surface_albedos"
path = out / "surface_spectra_and_fits.csv"
rows = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
surfaces = ("urban", "rural", "desert", "forest")
bands = ("o2a", "weak_co2", "strong_co2")

fig, axes = plt.subplots(4, 3, figsize=(14, 12), constrained_layout=True)
for i, surface in enumerate(surfaces):
    for j, band in enumerate(bands):
        ax = axes[i, j]
        use = (rows["surface"] == surface) & (rows["band"] == band)
        wl = rows["wavelength_nm"][use]
        order = np.argsort(wl)
        ax.plot(wl[order], rows["ecostress_mix"][use][order], label="ECOSTRESS mixture")
        ax.plot(wl[order], rows["P0P2_fit"][use][order], "--", label="P0–P2 fit")
        ax.grid(alpha=0.25)
        ax.set_title(f"{surface.title()} — {band.replace('_', ' ').upper()}")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Lambertian reflectance")
        if i == 0 and j == 0:
            ax.legend()
fig.savefig(out / "surface_legendre_fits.png", dpi=180)
plt.close(fig)
