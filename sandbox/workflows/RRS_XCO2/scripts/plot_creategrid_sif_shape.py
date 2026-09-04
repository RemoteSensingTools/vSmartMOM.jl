#!/usr/bin/env python3
"""Plot the SIF_OLD shape with the RRS-XCO2 campaign normalization.

The historical creategrid drivers normalized this template by its peak.  The
truth campaign instead preserves the shape while imposing
``2*pi*L_lambda(760 nm) = 0.5 mW m^-2 nm^-1``.  This script mirrors that
campaign convention; it does not change the generic integrated-SIF loader.
"""
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
reference_wavelength_nm = 760.0
angular_integral_760 = 0.5
target_radiance_760 = angular_integral_760 / (2.0 * np.pi)
template_at_760 = np.interp(reference_wavelength_nm, wl_nm, sif_nm)
if not np.isfinite(template_at_760) or template_at_760 <= 0.0:
    raise RuntimeError("SIF_OLD is nonpositive or invalid at 760 nm")
sif_nm *= target_radiance_760 / template_at_760
nu_cm1 = 1e7 / wl_nm
# This is the exact Jacobian conversion used by load_sif_spectrum().
sif_cm1 = sif_nm * 1e7 / nu_cm1**2

np.savetxt(
    OUT / "creategrid_sif_old_shape.csv",
    np.column_stack((wl_nm, nu_cm1, sif_nm, sif_cm1)),
    delimiter=",",
    header=(
        "wavelength_nm,wavenumber_cm-1,"
        "SIF_radiance_mW_m-2_sr-1_nm-1,"
        "SIF_radiance_mW_m-2_sr-1_per_cm-1"
    ),
    comments="",
)

fig, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
axes[0].plot(wl_nm, sif_nm, lw=1.5)
axes[0].axvspan(757, 773, color="tab:orange", alpha=0.18, label="OCO O2 A band")
axes[0].axvline(reference_wavelength_nm, color="0.25", ls="--", lw=1.0)
axes[0].scatter(
    [reference_wavelength_nm], [target_radiance_760], color="0.15", s=24,
    zorder=3, label=r"$L_{760}=0.5/(2\pi)$",
)
axes[0].set_ylabel(r"SIF radiance (mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)")
axes[0].legend()
axes[0].grid(alpha=0.25)

axes[1].plot(wl_nm, sif_cm1, lw=1.5)
axes[1].axvspan(757, 773, color="tab:orange", alpha=0.18)
axes[1].set_xlabel("Wavelength (nm)")
axes[1].set_ylabel(r"SIF radiance (mW m$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$)")
axes[1].grid(alpha=0.25)
fig.suptitle(
    r"SIF_OLD campaign shape: $2\pi L_\lambda(760\,\mathrm{nm})=0.5$"
)
fig.savefig(OUT / "creategrid_sif_old_shape.png", dpi=180)
plt.close(fig)
