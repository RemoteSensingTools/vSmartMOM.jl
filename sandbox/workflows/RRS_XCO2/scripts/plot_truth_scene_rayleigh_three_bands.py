#!/usr/bin/env python3
"""Plot three-band Rayleigh/noRS spectra with synthetic OCO overlays."""

import argparse
from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_truth_rayleigh_three_bands")

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

from oco_noise_clouds import (
    assert_same_grid,
    per_nm_to_plot_per_cm,
    read_oco_noise,
)


TRUTH_ROOT = Path(__file__).resolve().parents[1] / "truth_map"
DEFAULT_SCENE_ROOTS = (TRUTH_ROOT, TRUTH_ROOT / "aerosol_chunked")
DEFAULT_OUTPUT = TRUTH_ROOT / "rayleigh_three_bands"
DEFAULT_OCO_ROOT = TRUTH_ROOT / "OCO_radiances"
DEFAULT_NOISE_ROOT = DEFAULT_OCO_ROOT / "noise_covariances"

BANDS = (
    ("o2a", "O₂ A-band", "radiance_rayleigh_o2a"),
    ("weak_co2", "Weak CO₂ band", "radiance_rayleigh_weak_co2"),
    ("strong_co2", "Strong CO₂ band", "radiance_rayleigh_strong_co2"),
)


def stokes_i(values):
    values = np.ma.asarray(values)
    if values.ndim != 2:
        raise ValueError(f"expected a two-dimensional Stokes spectrum, got {values.shape}")
    if values.shape[0] == 3:
        return values[0, :]
    if values.shape[1] == 3:
        return values[:, 0]
    raise ValueError(f"cannot locate the three-element Stokes axis in {values.shape}")


def read_wavelengths(scene_root):
    with Dataset(scene_root / "sim_wavelength.nc") as ds:
        return {
            band: np.asarray(ds[f"{band}_wavelength"][:], dtype=float)
            for band, _, _ in BANDS
        }


def read_oco_rayleigh(state, oco_root):
    path = oco_root / f"OCO2sims_{state:03d}.nc"
    if not path.is_file():
        raise FileNotFoundError(f"missing synthetic OCO radiance file: {path}")
    result = {}
    with Dataset(path) as ds:
        for band, _, _ in BANDS:
            wavelength = np.asarray(ds[f"{band}_wavelength"][:], dtype=float)
            radiance_per_nm = np.asarray(
                ds[f"I_OCO_rayleigh_{band}"][:], dtype=float
            )
            radiance_per_cm = per_nm_to_plot_per_cm(
                radiance_per_nm, wavelength
            )
            result[band] = (wavelength, radiance_per_cm)
    return result


def scene_complete(scene):
    with Dataset(scene) as ds:
        complete = int(getattr(ds, "simulation_complete", 0)) == 1 or int(
            getattr(ds, "chunked_simulation_complete", 0)
        ) == 1
        if not complete:
            return False
        for _, _, variable_name in BANDS:
            values = ds[variable_name][:]
            if np.ma.count(values) != values.size:
                return False
            if not np.isfinite(np.asarray(values)).all():
                return False
    return True


def discover_scenes(scene_roots):
    scenes = {}
    for scene_root in scene_roots:
        if not scene_root.is_dir():
            continue
        for scene in sorted(scene_root.glob("hiressim_*.nc")):
            if not scene_complete(scene):
                continue
            with Dataset(scene) as ds:
                state = int(ds.getncattr("state_index"))
            if state in scenes:
                raise RuntimeError(
                    f"duplicate state {state:03d}: {scenes[state]} and {scene}"
                )
            scenes[state] = scene
    return [scenes[state] for state in sorted(scenes)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scene-root",
        type=Path,
        action="append",
        help="Scene directory; repeat to combine roots (default: clear + aerosol)",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--oco-root", type=Path, default=DEFAULT_OCO_ROOT)
    parser.add_argument(
        "--noise-root", type=Path, default=DEFAULT_NOISE_ROOT,
        help="Directory containing OCO2noise_NNN.nc files",
    )
    args = parser.parse_args()
    scene_roots = tuple(
        path.resolve() for path in (args.scene_root or DEFAULT_SCENE_ROOTS)
    )
    output_root = args.output.resolve()
    oco_root = args.oco_root.resolve()
    noise_root = args.noise_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    wavelength_cache = {}
    for scene_root in scene_roots:
        wavelength_cache[scene_root] = read_wavelengths(scene_root)

    scenes = discover_scenes(scene_roots)
    if not scenes:
        raise RuntimeError(f"no completed truth scenes found under {scene_roots}")

    generated = 0
    for scene in scenes:
        scene_root = scene.parent.resolve()
        wavelengths = wavelength_cache[scene_root]
        with Dataset(scene) as ds:
            state = int(ds.getncattr("state_index"))
            surface = ds.getncattr("surface")
            aerosol = ds.getncattr("aerosol_case")
            sif = ds.getncattr("sif_case")
            xco2 = float(ds.getncattr("xco2_ppm"))
            spectra = {
                band: stokes_i(ds[varname][:])
                for band, _, varname in BANDS
            }

        oco = read_oco_rayleigh(state, oco_root)
        fig, axes = plt.subplots(3, 1, figsize=(13, 11), constrained_layout=True)
        for panel_index, (ax, (band, title, _)) in enumerate(zip(axes, BANDS)):
            wavelength = wavelengths[band]
            radiance = spectra[band]
            order = np.argsort(wavelength)
            ax.plot(
                wavelength[order],
                radiance[order],
                color="#1f77b4",
                linewidth=0.8,
                label="High-resolution Stokes I",
            )
            oco_wavelength, oco_radiance = oco[band]
            noise = read_oco_noise(state, noise_root, band)
            assert_same_grid(oco_wavelength, noise["wavelength"])
            oco_order = np.argsort(oco_wavelength)
            ax.fill_between(
                oco_wavelength[oco_order],
                (oco_radiance - noise["corrected"])[oco_order],
                (oco_radiance + noise["corrected"])[oco_order],
                color="0.45",
                alpha=0.28,
                linewidth=0,
                zorder=1,
                label="corrected OCO radiance ±1σ",
            )
            ax.plot(
                oco_wavelength[oco_order],
                oco_radiance[oco_order],
                color="black",
                linewidth=1.6,
                zorder=5,
                label="OCO analyzer + Gaussian ILS; divided by M₁₁=0.5",
            )
            ax.set_title(title, loc="left")
            ax.set_xlabel("Wavelength [nm]")
            ax.set_ylabel("Stokes I\n[mW m⁻² sr⁻¹ (cm⁻¹)⁻¹]")
            ax.grid(alpha=0.22)
            ax.legend(loc="best", frameon=False, fontsize=8)

        fig.suptitle(
            f"Rayleigh/noRS truth scene {state:03d}: {surface}; "
            f"aerosol={aerosol}; SIF={sif}; XCO₂={xco2:g} ppm",
            fontsize=14,
        )
        output = output_root / f"state{state:03d}_rayleigh_three_bands.png"
        fig.savefig(output, dpi=220)
        plt.close(fig)
        print(output)
        generated += 1

    print(
        f"Generated {generated} Rayleigh/noRS three-band plots in {output_root}"
    )


if __name__ == "__main__":
    main()
