#!/usr/bin/env python3
"""Plot matched no-SIF O2 A-band differences between surface types."""

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_truth_surface_differences")

from netCDF4 import Dataset
import numpy as np

from plot_truth_state_o2_components import (
    ROOT,
    WAVELENGTHS,
    component_panels,
    draw_panels,
    o2_is_complete,
    read_scene,
)


SURFACE_PAIRS = (
    ("urban", "rural"),
    ("urban", "desert"),
    ("urban", "forest"),
    ("rural", "desert"),
    ("rural", "forest"),
    ("desert", "forest"),
)


def main():
    with Dataset(WAVELENGTHS) as ds:
        wavelength = np.asarray(ds["o2a_wavelength"][:])

    scenes = {}
    for path in sorted(ROOT.glob("hiressim_*.nc")):
        if not o2_is_complete(path):
            continue
        item = read_scene(path)
        if item["sif_case"] != "off":
            continue
        key = (item["surface"], float(item["xco2"]), item["aerosol"])
        if key in scenes:
            raise RuntimeError(f"duplicate completed no-SIF scene for {key}")
        scenes[key] = item

    output_dir = ROOT / "surface_differences_no_sif"
    output_dir.mkdir(exist_ok=True)
    generated = 0

    aerosol_cases = sorted({key[2] for key in scenes})
    xco2_values = sorted({key[1] for key in scenes})
    for first_surface, second_surface in SURFACE_PAIRS:
        for aerosol in aerosol_cases:
            for xco2 in xco2_values:
                first = scenes.get((first_surface, xco2, aerosol))
                second = scenes.get((second_surface, xco2, aerosol))
                if first is None or second is None:
                    continue

                difference = {
                    name: first[name] - second[name]
                    for name in ("rayleigh", "cabannes", "rrs")
                }
                xco2_label = f"{xco2:g}"
                output = output_dir / (
                    f"{first_surface}_minus_{second_surface}_"
                    f"xco2_{xco2_label}ppm_{aerosol}_o2_components.png"
                )
                draw_panels(
                    wavelength,
                    component_panels(difference),
                    f"O₂ A-band no-SIF surface difference: "
                    f"{first_surface} − {second_surface}; "
                    f"XCO₂={xco2_label} ppm; aerosol={aerosol}",
                    output,
                    "Δ Stokes I\n[mW m⁻² sr⁻¹ (cm⁻¹)⁻¹]",
                    explicit_scientific=True,
                )
                generated += 1

    if generated == 0:
        raise RuntimeError("no matched completed no-SIF surface pairs found")
    print(f"Generated {generated} matched surface-difference plots in {output_dir}")


if __name__ == "__main__":
    main()
