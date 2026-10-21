#!/usr/bin/env python3
"""Plot the archived versus coordinate-corrected aerosol A-band components.

The corrected 2026-08 aerosol truth was evaluated in 64-point core chunks.
The historical shoulder construction rebuilt each core as part of a longer
Float32 range; consequently, some retained solve nodes differ from the
canonical output nodes by one Float32 ULP.  This diagnostic reconstructs those
actual retained nodes and linearly maps the saved spectra back to the nominal
grid.  It does not overwrite the high-resolution truth files.
"""

import argparse
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import netCDF4
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCENE = ROOT / "truth_map" / "aerosol_chunked" / "hiressim_009.nc"
DEFAULT_GRID = ROOT / "truth_map" / "aerosol_chunked" / "sim_wavelength.nc"


def legacy_retained_nodes(core: np.ndarray, shoulder_cm: float = 234.0) -> np.ndarray:
    """Reproduce Julia's historical Float32 shoulder-origin range exactly."""
    step = np.float32(0.1)
    start = np.float32(np.float32(core[0]) - np.float32(shoulder_cm))
    shoulder_points = int(round(shoulder_cm / float(step)))
    return np.asarray(
        [np.float32(start + np.float32(i) * step)
         for i in range(shoulder_points, shoulder_points + len(core))],
        dtype=np.float32,
    )


def correct_legacy_chunks(values: np.ndarray, nominal_nu: np.ndarray,
                          chunk_points: int = 64) -> np.ndarray:
    corrected = np.empty_like(values, dtype=np.float64)
    nominal32 = np.asarray(nominal_nu, dtype=np.float32)
    for first in range(0, len(nominal32), chunk_points):
        indices = slice(first, min(first + chunk_points, len(nominal32)))
        target = nominal32[indices]
        actual = legacy_retained_nodes(target)
        for stokes in range(values.shape[0]):
            corrected[stokes, indices] = np.interp(
                target.astype(np.float64), actual.astype(np.float64),
                values[stokes, indices].astype(np.float64),
            )
    return corrected


def read_scene(path: Path) -> Dict[str, np.ndarray]:
    with netCDF4.Dataset(path) as dataset:
        spectra = {
            "rayleigh": np.asarray(dataset["radiance_rayleigh_o2a"][:]),
            "cabannes": np.asarray(dataset["radiance_cabannes_o2a"][:]),
            "rrs": np.asarray(dataset["radiance_rrs_o2a"][:]),
        }
    # NCDatasets.jl and netCDF4-Python expose the on-disk dimension ordering
    # oppositely for files written by the Julia workflow. Normalize to
    # (stokes, spectral) here.
    return {name: values if values.shape[0] == 3 else values.T
            for name, values in spectra.items()}


def read_grid(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    with netCDF4.Dataset(path) as dataset:
        nu = np.asarray(dataset["o2a_wavenumber"][:], dtype=np.float64)
        wavelength = np.asarray(dataset["o2a_wavelength"][:], dtype=np.float64)
    return nu, wavelength


def components(spectra: Dict[str, np.ndarray]) -> Tuple[np.ndarray, ...]:
    rayleigh = spectra["rayleigh"]
    cabannes_minus_rayleigh = spectra["cabannes"] - rayleigh
    rrs = spectra["rrs"]
    rrs_plus_cabannes_minus_rayleigh = rrs + cabannes_minus_rayleigh
    return (rayleigh, cabannes_minus_rayleigh, rrs,
            rrs_plus_cabannes_minus_rayleigh)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scene", type=Path, default=DEFAULT_SCENE)
    parser.add_argument("--grid", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--chunk-points", type=int, default=64)
    args = parser.parse_args()

    output = args.output or args.scene.with_name(
        f"{args.scene.stem}_o2_components_before_after_grid_correction.png")
    raw = read_scene(args.scene)
    nominal_nu, wavelength_nm = read_grid(args.grid)
    corrected = {
        name: correct_legacy_chunks(values, nominal_nu, args.chunk_points)
        for name, values in raw.items()
    }

    labels = (
        "Rayleigh",
        "Cabannes − Rayleigh",
        "RRS",
        "Cabannes + RRS − Rayleigh",
    )
    colors = ("#2166ac", "#b2182b", "#1b7837", "#762a83")
    raw_components = components(raw)
    corrected_components = components(corrected)
    # Retain the established truth-map plotting convention. This factor is
    # display-only and is never applied to measurements or their covariance.
    display_normalization = 0.5
    order = np.argsort(wavelength_nm)

    figure, axes = plt.subplots(
        4, 2, figsize=(14, 12), sharex=True, constrained_layout=True)
    for row, (label, color, before, after) in enumerate(zip(
            labels, colors, raw_components, corrected_components)):
        before_i = before[0, order] / display_normalization
        after_i = after[0, order] / display_normalization
        delta_i = after_i - before_i
        axes[row, 0].plot(wavelength_nm[order], before_i, color="0.55",
                          linewidth=0.9, label="archived")
        axes[row, 0].plot(wavelength_nm[order], after_i, color=color,
                          linewidth=0.8, label="grid-corrected")
        axes[row, 1].plot(wavelength_nm[order], delta_i, color=color,
                          linewidth=0.8)
        axes[row, 1].axhline(0.0, color="0.3", linewidth=0.6)
        axes[row, 0].set_ylabel(f"{label}\nI / 0.5")
        axes[row, 1].set_ylabel("after − before\nI / 0.5")
        axes[row, 0].grid(alpha=0.2)
        axes[row, 1].grid(alpha=0.2)
        axes[row, 1].text(
            0.01, 0.94,
            f"max |Δ| = {np.max(np.abs(delta_i)):.3e}",
            transform=axes[row, 1].transAxes, va="top", fontsize=9)

    axes[0, 0].legend(loc="best", frameon=False)
    axes[-1, 0].set_xlabel("Wavelength (nm)")
    axes[-1, 1].set_xlabel("Wavelength (nm)")
    figure.suptitle(
        f"{args.scene.stem}: aerosol O₂ A-band legacy-grid correction\n"
        "64-point cores, ±234 cm⁻¹ Raman shoulders; Stokes I",
        fontsize=14,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
