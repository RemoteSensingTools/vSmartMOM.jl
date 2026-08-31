#!/usr/bin/env python3
"""Plot matched synthetic-OCO CO2-band differences inside noise clouds.

The default comparison is the canonical urban, aerosol-free, SIF-off nadir
scene at 420 ppm minus the otherwise identical 400 ppm scene. Other truth-map
scene classes and CO2 pairs can be selected from the command line.

As in the existing aerosol and SIF difference plots, the high-XCO2 scene is
treated as one noisy measurement and the low-XCO2 scene as a deterministic
forward-model reference. The zero-centered cloud is therefore +/- the noise
standard deviation of the high-XCO2 scene, not the quadrature sum of two
independent measurement noises.
"""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_xco2_scene_differences")

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

from oco_noise_clouds import (
    assert_same_grid,
    per_nm_to_plot_per_cm,
    read_oco_noise,
)


TRUTH_ROOT = Path(__file__).resolve().parents[1] / "truth_map"
DEFAULT_OCO_ROOT = TRUTH_ROOT / "OCO_radiances"
DEFAULT_NOISE_ROOT = DEFAULT_OCO_ROOT / "noise_covariances"
DEFAULT_OUTPUT_ROOT = DEFAULT_OCO_ROOT / "xco2_differences"

BANDS = (
    ("weak_co2", "Weak CO₂ band"),
    ("strong_co2", "Strong CO₂ band"),
)
MATCH_ATTRIBUTES = (
    "surface",
    "aerosol_case",
    "sif_case",
    "psurf_hpa",
    "sza_deg",
    "vza_deg",
    "atmospheric_layers",
)


def scene_metadata(path):
    """Return the metadata needed to select and validate one OCO scene."""
    with Dataset(path) as dataset:
        if int(getattr(dataset, "instrument_processing_complete", 0)) != 1:
            raise RuntimeError(f"incomplete synthetic OCO file: {path}")
        metadata = {
            "state": int(dataset.getncattr("state_index")),
            "xco2_ppm": float(dataset.getncattr("xco2_ppm")),
        }
        metadata.update(
            {name: dataset.getncattr(name) for name in MATCH_ATTRIBUTES}
        )
    metadata["path"] = path
    return metadata


def find_scene(oco_root, *, surface, aerosol, sif, xco2_ppm):
    """Find exactly one truth-map scene with the requested physical state."""
    matches = []
    for path in sorted(oco_root.glob("OCO2sims_*.nc")):
        metadata = scene_metadata(path)
        if (
            metadata["surface"] == surface
            and metadata["aerosol_case"] == aerosol
            and metadata["sif_case"] == sif
            and np.isclose(metadata["xco2_ppm"], xco2_ppm, rtol=0, atol=1e-10)
        ):
            matches.append(metadata)
    if len(matches) != 1:
        raise RuntimeError(
            "expected exactly one matching OCO scene for "
            f"surface={surface}, aerosol={aerosol}, SIF={sif}, "
            f"XCO2={xco2_ppm:g} ppm; found {len(matches)}"
        )
    return matches[0]


def assert_only_xco2_differs(low, high):
    """Reject a comparison in which any non-CO2 scene attribute changes."""
    mismatches = []
    for name in MATCH_ATTRIBUTES:
        low_value = low[name]
        high_value = high[name]
        if isinstance(low_value, (float, np.floating)):
            same = np.isclose(low_value, high_value, rtol=0, atol=1e-10)
        else:
            same = low_value == high_value
        if not same:
            mismatches.append(f"{name}: {low_value!r} != {high_value!r}")
    if mismatches:
        raise RuntimeError(
            "selected states differ in parameters other than XCO2:\n  "
            + "\n  ".join(mismatches)
        )


def read_band(scene, band):
    """Read plotting-normalized Rayleigh/noRS OCO radiance for one band."""
    with Dataset(scene["path"]) as dataset:
        wavelength = np.asarray(dataset[f"{band}_wavelength"][:], dtype=float)
        radiance_per_nm = np.asarray(
            dataset[f"I_OCO_rayleigh_{band}"][:], dtype=float
        )
    if not np.isfinite(wavelength).all() or not np.isfinite(radiance_per_nm).all():
        raise RuntimeError(
            f"non-finite {band} data in synthetic OCO state {scene['state']:03d}"
        )
    return wavelength, per_nm_to_plot_per_cm(radiance_per_nm, wavelength)


def safe_token(value):
    return f"{value:g}".replace("-", "m").replace(".", "p")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", default="urban")
    parser.add_argument("--aerosol", default="none")
    parser.add_argument("--sif", default="off")
    parser.add_argument("--xco2-low", type=float, default=400.0)
    parser.add_argument("--xco2-high", type=float, default=420.0)
    parser.add_argument("--oco-root", type=Path, default=DEFAULT_OCO_ROOT)
    parser.add_argument("--noise-root", type=Path, default=DEFAULT_NOISE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.xco2_high == args.xco2_low:
        raise ValueError("the two XCO2 values must differ")

    oco_root = args.oco_root.resolve()
    noise_root = args.noise_root.resolve()
    output_root = args.output_root.resolve()

    low = find_scene(
        oco_root,
        surface=args.surface,
        aerosol=args.aerosol,
        sif=args.sif,
        xco2_ppm=args.xco2_low,
    )
    high = find_scene(
        oco_root,
        surface=args.surface,
        aerosol=args.aerosol,
        sif=args.sif,
        xco2_ppm=args.xco2_high,
    )
    assert_only_xco2_differs(low, high)

    if args.output is None:
        output_root.mkdir(parents=True, exist_ok=True)
        output = output_root / (
            f"state{high['state']:03d}_minus{low['state']:03d}_"
            f"xco2_{safe_token(args.xco2_high)}_minus_"
            f"{safe_token(args.xco2_low)}ppm_co2_noise_clouds.png"
        )
    else:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(13.5, 8.5), constrained_layout=True)
    metric_lines = []
    for axis, (band, band_title) in zip(axes, BANDS):
        wavelength_low, radiance_low = read_band(low, band)
        wavelength_high, radiance_high = read_band(high, band)
        assert_same_grid(
            wavelength_low,
            wavelength_high,
            label=f"matched {band} OCO radiance",
        )

        noise = read_oco_noise(high["state"], noise_root, band)
        assert_same_grid(wavelength_high, noise["wavelength"])
        sigma = noise["corrected"]
        difference = radiance_high - radiance_low
        normalized = difference / sigma
        normalized_rms = float(np.sqrt(np.mean(normalized**2)))
        maximum_absolute = float(np.max(np.abs(normalized)))
        fraction_above_one = float(np.mean(np.abs(normalized) > 1.0))
        metric_lines.append(
            f"{band}: RMS(DeltaI/sigma)={normalized_rms:.3f}, "
            f"max={maximum_absolute:.3f}, "
            f"fraction |DeltaI|>sigma={fraction_above_one:.3f}"
        )

        order = np.argsort(wavelength_high)
        wavelength_plot = wavelength_high[order]
        difference_plot = difference[order]
        sigma_plot = sigma[order]
        axis.fill_between(
            wavelength_plot,
            -sigma_plot,
            sigma_plot,
            color="0.50",
            alpha=0.30,
            linewidth=0,
            label=f"±1σ noise of the {args.xco2_high:g}-ppm measurement",
        )
        axis.plot(
            wavelength_plot,
            difference_plot,
            color="#b2182b",
            linewidth=1.0,
            label=(
                f"ΔI = I({args.xco2_high:g} ppm) − "
                f"I({args.xco2_low:g} ppm)"
            ),
        )
        axis.axhline(0.0, color="black", linewidth=0.65, alpha=0.75)
        axis.text(
            0.012,
            0.96,
            f"RMS(ΔI/σ) = {normalized_rms:.2f}\n"
            f"max |ΔI|/σ = {maximum_absolute:.2f}\n"
            f"samples with |ΔI| > σ = {100*fraction_above_one:.1f}%",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
        )
        # Center the panel title so Matplotlib's scientific-notation exponent
        # at the upper-left of the y axis cannot overlap it.
        axis.set_title(band_title, loc="center", pad=8)
        axis.set_xlabel("Wavelength [nm]")
        axis.set_ylabel("Δ OCO radiance\n[mW m⁻² sr⁻¹ (cm⁻¹)⁻¹]")
        axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axis.grid(alpha=0.22)
        axis.legend(loc="lower right", frameon=False, fontsize=8.5)

    fig.suptitle(
        f"Synthetic OCO CO₂ response: {args.xco2_high:g} − "
        f"{args.xco2_low:g} ppm\n"
        f"states {high['state']:03d} − {low['state']:03d}; "
        f"surface={high['surface']}; aerosol={high['aerosol_case']}; "
        f"SIF={high['sif_case']}; pₛ={float(high['psurf_hpa']):g} hPa; "
        f"SZA={float(high['sza_deg']):g}°; VZA={float(high['vza_deg']):g}°",
        fontsize=13,
    )
    fig.savefig(output, dpi=220)
    plt.close(fig)

    print(output)
    for line in metric_lines:
        print(line)


if __name__ == "__main__":
    main()
