#!/usr/bin/env python3
"""Plot aerosol-minus-no-aerosol O2 A-band components for every surface."""

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from plot_truth_state_o2_components import (
    NOISE_ROOT,
    component_panels,
    draw_panels,
    o2_is_complete,
    read_oco_scene,
    read_scene,
)
from oco_noise_clouds import effect_scene_o2_difference_clouds, read_oco_noise


SURFACE_ORDER = ("urban", "rural", "desert", "forest")


def completed_scenes(truth_root):
    paths = list(truth_root.glob("hiressim_*.nc"))
    paths += list((truth_root / "aerosol_chunked").glob("hiressim_*.nc"))
    return [read_scene(path) for path in sorted(paths) if o2_is_complete(path)]


def main():
    default_root = Path(__file__).resolve().parents[1] / "truth_map"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-root", type=Path, default=default_root)
    parser.add_argument("--xco2", type=float, default=380.0)
    parser.add_argument(
        "--noise-root",
        type=Path,
        default=NOISE_ROOT,
        help="Directory containing OCO2noise_NNN.nc files",
    )
    args = parser.parse_args()
    truth_root = args.truth_root.resolve()
    noise_root = args.noise_root.resolve()

    with Dataset(truth_root / "sim_wavelength.nc") as dataset:
        wavelength = np.asarray(dataset["o2a_wavelength"][:])

    scenes = completed_scenes(truth_root)
    selected = {}
    for scene in scenes:
        if scene["sif_case"] != "off" or float(scene["xco2"]) != args.xco2:
            continue
        aerosol_key = "none" if scene["aerosol"] == "none" else "aerosol"
        key = (scene["surface"], aerosol_key)
        if key in selected:
            raise RuntimeError(f"duplicate scene for {key} at XCO2={args.xco2:g}")
        selected[key] = scene

    outputs = []
    for surface in SURFACE_ORDER:
        clear = selected.get((surface, "none"))
        aerosol = selected.get((surface, "aerosol"))
        if clear is None or aerosol is None:
            raise RuntimeError(
                f"missing matched {surface} pair at XCO2={args.xco2:g} ppm"
            )
        difference = {
            name: aerosol[name] - clear[name]
            for name in ("rayleigh", "cabannes", "rrs")
        }
        oco_clear = read_oco_scene(clear["state"], truth_root / "OCO_radiances")
        oco_aerosol = read_oco_scene(
            aerosol["state"], truth_root / "OCO_radiances"
        )
        noise_aerosol = read_oco_noise(
            aerosol["state"], noise_root, "o2a"
        )
        oco_difference = {
            name: oco_aerosol[name] - oco_clear[name]
            for name in ("rayleigh", "cabannes", "rrs")
        }
        output = truth_root / (
            f"state{aerosol['state']:03d}_minus{clear['state']:03d}_"
            "aerosol_difference_o2_components.png"
        )
        draw_panels(
            wavelength,
            component_panels(difference),
            f"O₂ A-band aerosol effect: state {aerosol['state']:03d} − "
            f"state {clear['state']:03d} ({surface}, SIF off, "
            f"XCO₂={args.xco2:g} ppm; AOD₇₆₀=0.28 minus no aerosol)",
            output,
            "Δ Stokes I\n[mW m⁻² sr⁻¹ (cm⁻¹)⁻¹]",
            explicit_scientific=True,
            oco_overlay=(
                oco_aerosol["wavelength"],
                component_panels(oco_difference),
            ),
            noise_cloud=(
                noise_aerosol["wavelength"],
                effect_scene_o2_difference_clouds(noise_aerosol),
            ),
        )
        outputs.append(output)

    print(f"Generated {len(outputs)} matched aerosol-effect plots")


if __name__ == "__main__":
    main()
