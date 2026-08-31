#!/usr/bin/env python3
"""Plot forward spectra written by run_forward.jl (requires h5py/matplotlib)."""
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[1]
names = ("o2a", "weak_co2", "strong_co2")
fig, axes = plt.subplots(3, 1, figsize=(9, 9), constrained_layout=True)

for ax, name in zip(axes, names):
    path = root / "output" / f"forward_{name}.jld2"
    with h5py.File(path, "r") as handle:
        wavelength = np.asarray(handle["wavelength_nm"]).squeeze()
        toa = np.asarray(handle["toa"])
    # Select Stokes I/view 1 regardless of whether the HDF5 reader reverses
    # Julia's array axes: the spectral axis is uniquely identified by length.
    spectral_axes = [i for i, size in enumerate(toa.shape) if size == wavelength.size]
    if len(spectral_axes) != 1:
        raise ValueError(f"Cannot identify spectral axis in TOA shape {toa.shape}")
    intensity = np.moveaxis(toa, spectral_axes[0], -1).reshape(-1, wavelength.size)[0]
    order = np.argsort(wavelength)
    ax.plot(wavelength[order], intensity[order], lw=1)
    ax.set(title=name.replace("_", " ").upper(), ylabel="TOA Stokes I")
    ax.grid(alpha=0.25)
axes[-1].set_xlabel("Wavelength (nm)")
fig.savefig(root / "output" / "forward_spectra.png", dpi=180)
