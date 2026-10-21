#!/usr/bin/env python3
"""Read and combine truth-map OCO noise for plotting-only noise clouds.

The covariance products are stored in raw measurement units
(`mW m-2 sr-1 nm-1`).  Existing truth-map figures divide the OCO analyzer
output by its unpolarized throughput M11=0.5 and express spectral density per
cm-1.  This module applies the same two plotting-only conversions to the noise
standard deviation.  It never changes the measurement vector or covariance
files used by the retrieval.

The truth-map difference figures represent one noisy, effect-present
measurement compared with a deterministic reference simulation. Their
zero-centered cloud therefore uses the noise of the effect-present OCO
radiance alone. `independent_noise` remains available for diagnostics that
really do subtract statistically independent measurements.
"""

from pathlib import Path

from netCDF4 import Dataset
import numpy as np


PLOT_ANALYZER_M11 = 0.5
BAND_INDEX = {"o2a": 1, "weak_co2": 2, "strong_co2": 3}


def per_nm_to_plot_per_cm(values, wavelength_nm):
    """Apply the inverse density Jacobian and plotting M11 normalization."""
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    values = np.asarray(values, dtype=float)
    return values * wavelength_nm**2 / 1.0e7 / PLOT_ANALYZER_M11


def read_oco_noise(state, noise_root, band="o2a"):
    """Return corrected/uncorrected 1-sigma noise on one synthetic OCO grid."""
    if band not in BAND_INDEX:
        raise ValueError(f"unknown synthetic OCO band: {band}")
    path = Path(noise_root) / f"OCO2noise_{state:03d}.nc"
    if not path.is_file():
        raise FileNotFoundError(f"missing synthetic OCO noise file: {path}")

    with Dataset(path) as ds:
        if int(getattr(ds, "noise_covariance_complete", 0)) != 1:
            raise RuntimeError(f"incomplete synthetic OCO noise file: {path}")
        stored_state = int(ds.getncattr("state_index"))
        if stored_state != state:
            raise RuntimeError(
                f"noise file {path} stores state {stored_state}, expected {state}"
            )
        band_index = np.asarray(ds["band_index"][:], dtype=int)
        selected = band_index == BAND_INDEX[band]
        wavelength = np.asarray(ds["wavelength"][:], dtype=float)[selected]
        corrected = np.asarray(
            ds["noise_std_corrected"][:], dtype=float
        )[selected]
        uncorrected = np.asarray(
            ds["noise_std_uncorrected"][:], dtype=float
        )[selected]

    if wavelength.size == 0:
        raise RuntimeError(f"noise file {path} has no samples for band {band}")
    if not (
        np.isfinite(wavelength).all()
        and np.isfinite(corrected).all()
        and np.isfinite(uncorrected).all()
        and np.all(corrected > 0)
        and np.all(uncorrected > 0)
    ):
        raise RuntimeError(f"invalid noise samples in {path} for band {band}")

    return {
        "wavelength": wavelength,
        "corrected": per_nm_to_plot_per_cm(corrected, wavelength),
        "uncorrected": per_nm_to_plot_per_cm(uncorrected, wavelength),
    }


def independent_noise(*standard_deviations):
    """Propagate independent one-sigma terms through a signed difference."""
    if not standard_deviations:
        raise ValueError("at least one standard deviation is required")
    arrays = [np.asarray(sigma, dtype=float) for sigma in standard_deviations]
    shape = arrays[0].shape
    if any(array.shape != shape for array in arrays[1:]):
        raise ValueError("noise arrays must share a spectral grid")
    return np.sqrt(sum(array**2 for array in arrays))


def assert_same_grid(first, second, *, label="OCO noise"):
    """Reject silent interpolation between supposedly identical OCO grids."""
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if first.shape != second.shape or not np.allclose(first, second, rtol=0, atol=1e-10):
        raise ValueError(f"{label} wavelength grid does not match the OCO radiance grid")


def single_scene_o2_clouds(oco, noise):
    """Noise-cloud definitions for the four standard O2 component panels."""
    assert_same_grid(oco["wavelength"], noise["wavelength"])
    zero = np.zeros_like(noise["corrected"])
    return (
        (
            np.asarray(oco["rayleigh"], dtype=float),
            noise["corrected"],
            "corrected OCO radiance ±1σ",
        ),
        (
            zero,
            noise["uncorrected"],
            "±1σ uncorrected-measurement noise scale",
        ),
        (
            zero,
            noise["uncorrected"],
            "±1σ uncorrected-measurement noise scale",
        ),
        (
            zero,
            noise["uncorrected"],
            "±1σ uncorrected OCO noise",
        ),
    )


def effect_scene_o2_difference_clouds(noisy_scene):
    """Noise clouds for an effect-present scene minus a model reference.

    Only the effect-present truth radiance is treated as noisy. The reference
    (no aerosol or no SIF) is a deterministic forward-model spectrum. The
    Rayleigh panel therefore uses the corrected measurement covariance and
    the other three panels use the uncorrected measurement covariance.
    """
    zero = np.zeros_like(noisy_scene["corrected"])
    return (
        (zero, noisy_scene["corrected"], "±1σ effect-scene corrected noise"),
        (
            zero,
            noisy_scene["uncorrected"],
            "±1σ effect-scene uncorrected noise",
        ),
        (
            zero,
            noisy_scene["uncorrected"],
            "±1σ effect-scene uncorrected noise",
        ),
        (
            zero,
            noisy_scene["uncorrected"],
            "±1σ effect-scene uncorrected noise",
        ),
    )
