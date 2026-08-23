#!/usr/bin/env python3
"""Make one 3-band forward plot and one 3-band plot per Jacobian column."""
from pathlib import Path
import re
import os

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("RRS_XCO2_OUTPUT_DIR", ROOT / "output"))
JACOBIAN_OUT = OUT / "basis_jacobian_three_bands"
BANDS = ("o2a", "weak_co2", "strong_co2")


def shared_state_label(label):
    """Map legacy band-local gas labels onto the shared retrieval state."""
    match = re.fullmatch(r"(?:o2a|weak_co2|strong_co2)_gas_(\d+)", label)
    if not match:
        return label
    index = int(match.group(1))
    if 1 <= index <= 12:
        return f"h2o_layer_{index:02d}"
    if 13 <= index <= 24:
        return f"co2_layer_{index - 12:02d}"
    raise ValueError(f"Unexpected legacy gas index in {label}")


forward = {}
linear = {}
for band in BANDS:
    fpath = OUT / f"forward_{band}.csv"
    if fpath.exists():
        data = np.loadtxt(fpath, delimiter=",")
        forward[band] = (data[:, 0], data[:, 1])
    lpath = OUT / f"linearized_{band}.csv"
    if lpath.exists():
        with lpath.open() as stream:
            labels = [shared_state_label(label) for label in
                      stream.readline().strip().split(",")[1:]]
        data = np.loadtxt(lpath, delimiter=",", skiprows=1)
        linear[band] = (data[:, 0], labels, data[:, 1:].T)

if forward:
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), constrained_layout=True)
    for ax, band in zip(axes, BANDS):
        if band not in forward:
            ax.text(0.5, 0.5, "not available", ha="center", transform=ax.transAxes)
            continue
        wl, y = forward[band]
        order = np.argsort(wl)
        ax.plot(wl[order], y[order], lw=0.8)
        ax.set_title(band.replace("_", " ").upper())
        ax.set_ylabel("TOA Stokes I")
        ax.grid(alpha=0.25)
    axes[-1].set_xlabel("Wavelength (nm)")
    fig.savefig(OUT / "basis_forward_three_bands.png", dpi=180)
    plt.close(fig)

# Preserve retrieval-state order rather than sorting alphabetically. Labels
# first encountered in O2A, weak CO2, then strong CO2 define the plot order.
all_labels = []
for band in BANDS:
    if band in linear:
        for label in linear[band][1]:
            if label not in all_labels:
                all_labels.append(label)


def display_name(label):
    replacements = {
        "surface_pressure": "surface pressure",
        "tau_ref": "reference AOD",
        "n_real": "real refractive index",
        "n_imag": "imaginary refractive index",
        "median_radius": "median radius",
        "geometric_width": "geometric size width",
        "profile_location": "vertical-profile location",
        "profile_width": "vertical-profile width",
        "surface_P0": "surface Legendre P0 coefficient",
        "surface_P1": "surface Legendre P1 coefficient",
        "surface_P2": "surface Legendre P2 coefficient",
    }
    if label in replacements:
        return replacements[label]
    for token, text in replacements.items():
        if label.endswith("_" + token):
            species = label[: -(len(token) + 1)].replace("_", " ")
            return f"{species}: {text}"
    gas = re.fullmatch(r"(h2o|co2)_layer_(\d+)", label)
    if gas:
        return f"{gas.group(1).upper()}: profile layer {int(gas.group(2))}"
    return label.replace("_", " ")


if all_labels:
    JACOBIAN_OUT.mkdir(parents=True, exist_ok=True)
    # A regenerated state vector can contain fewer rows than an older run.
    # Remove only this script's prior row plots so stale derivatives cannot be
    # mistaken for output from the current benchmark.
    for old_plot in JACOBIAN_OUT.glob("jacobian_row_*.png"):
        old_plot.unlink()

for row_number, label in enumerate(all_labels, start=1):
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), constrained_layout=True)
    for ax, band in zip(axes, BANDS):
        ax.set_title(band.replace("_", " ").upper())
        if band not in linear or label not in linear[band][1]:
            if label in ("SIF760", "SIF_ref", "mSIF"):
                wl = forward[band][0] if band in forward else np.array([0.0, 1.0])
                ax.plot(wl, np.zeros_like(wl), lw=0.8)
            else:
                ax.text(0.5, 0.5, "parameter not present in this band",
                        ha="center", transform=ax.transAxes)
        else:
            wl, labels, jac = linear[band]
            y = jac[labels.index(label)]
            order = np.argsort(wl)
            ax.plot(wl[order], y[order], lw=0.8)
        ax.set_ylabel("dI/dx")
        ax.grid(alpha=0.25)
    axes[-1].set_xlabel("Wavelength (nm)")
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
    title = display_name(label)
    fig.suptitle(f"TOA Stokes-I Jacobian: ∂I/∂({title})")
    fig.savefig(JACOBIAN_OUT / f"jacobian_row_{row_number:03d}__{safe}.png",
                dpi=180)
    plt.close(fig)
