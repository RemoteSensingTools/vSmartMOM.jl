#!/usr/bin/env python3
"""Plot the SIF_OLD shape used by the creategrid*.jl benchmark drivers."""
from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "src" / "SIF_emission" / "sif-spectra.csv"
OUT = Path(os.environ.get("RRS_XCO2_OUTPUT_DIR", ROOT / "output"))
OUT.mkdir(parents=True, exist_ok=True)

raw = np.genfromtxt(SOURCE, delimiter=",", names=True)
wl_nm = raw["WL"]
sif_nm = raw["SIF_OLD"]
sif_nm *= 0.5 * np.pi / sif_nm.max()
nu_cm1 = 1e7 / wl_nm
# This is the exact Jacobian conversion used by load_sif_spectrum().
sif_cm1 = sif_nm * 1e7 / nu_cm1**2

np.savetxt(
    OUT / "creategrid_sif_old_shape.csv",
    np.column_stack((wl_nm, nu_cm1, sif_nm, sif_cm1)),
    delimiter=",",
    header="wavelength_nm,wavenumber_cm-1,SIF_normalized_per_nm,SIF_per_cm-1",
    comments="",
)

fig, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
axes[0].plot(wl_nm, sif_nm, lw=1.5)
axes[0].axvspan(757, 773, color="tab:orange", alpha=0.18, label="OCO O2 A band")
axes[0].set_ylabel("Normalized SIF per nm")
axes[0].legend()
axes[0].grid(alpha=0.25)

axes[1].plot(wl_nm, sif_cm1, lw=1.5)
axes[1].axvspan(757, 773, color="tab:orange", alpha=0.18)
axes[1].set_xlabel("Wavelength (nm)")
axes[1].set_ylabel("SIF passed to RT per cm$^{-1}$")
axes[1].grid(alpha=0.25)
fig.suptitle("SIF_OLD shape used by creategrid*.jl")
fig.savefig(OUT / "creategrid_sif_old_shape.png", dpi=180)
plt.close(fig)
