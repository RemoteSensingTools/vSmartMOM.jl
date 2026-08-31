#!/usr/bin/env python3
"""Plot O2 A-band components with synthetic OCO-radiance overlays."""

import argparse
from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_truth_o2_components")

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from netCDF4 import Dataset
import numpy as np

from oco_noise_clouds import (
    effect_scene_o2_difference_clouds,
    per_nm_to_plot_per_cm,
    read_oco_noise,
    single_scene_o2_clouds,
)

TRUTH_ROOT = Path(__file__).resolve().parents[1] / "truth_map"
ROOT = TRUTH_ROOT / "aerosol_chunked"
OCO_ROOT = TRUTH_ROOT / "OCO_radiances"
NOISE_ROOT = OCO_ROOT / "noise_covariances"


def stokes_i(a):
    a = np.asarray(a)
    if a.shape[0] == 3:
        return a[0]
    if a.shape[-1] == 3:
        return a[:, 0]
    raise ValueError(f"cannot locate Stokes dimension in {a.shape}")


def o2_is_complete(scene):
    """Require every value in all three O2 products to be written and finite."""
    with Dataset(scene) as ds:
        for name in ("radiance_rayleigh_o2a", "radiance_cabannes_o2a",
                     "radiance_rrs_o2a"):
            values = ds[name][:]
            if np.ma.count(values) != values.size:
                return False
            if not np.isfinite(np.asarray(values)).all():
                return False
    return True


def read_scene(scene):
    with Dataset(scene) as ds:
        return {
            "rayleigh": stokes_i(ds["radiance_rayleigh_o2a"][:]),
            "cabannes": stokes_i(ds["radiance_cabannes_o2a"][:]),
            "rrs": stokes_i(ds["radiance_rrs_o2a"][:]),
            "state": int(ds.getncattr("state_index")),
            "surface": ds.getncattr("surface"),
            "sif_case": ds.getncattr("sif_case"),
            "xco2": ds.getncattr("xco2_ppm"),
            "aerosol": ds.getncattr("aerosol_case"),
            "sza": float(ds.getncattr("sza_deg")),
            "vza": float(ds.getncattr("vza_deg")),
            "relative_azimuth": float(ds.getncattr("relative_azimuth_deg")),
        }


def read_oco_scene(state, oco_root=OCO_ROOT):
    """Read OCO intensities, normalize by M11, and restore per-cm-1 units."""
    path = Path(oco_root) / f"OCO2sims_{state:03d}.nc"
    if not path.is_file():
        raise FileNotFoundError(f"missing synthetic OCO radiance file: {path}")
    with Dataset(path) as ds:
        wavelength = np.asarray(ds["o2a_wavelength"][:], dtype=float)
        data = {
            component: per_nm_to_plot_per_cm(
                ds[f"I_OCO_{component}_o2a"][:], wavelength
            )
            for component in ("rayleigh", "cabannes", "rrs")
        }
    data["wavelength"] = wavelength
    return data


def component_panels(data):
    rayleigh = data["rayleigh"]
    cabannes = data["cabannes"]
    rrs = data["rrs"]
    return (
        (rayleigh, "Rayleigh (elastic noRS)", "#1f77b4"),
        (cabannes - rayleigh, "Cabannes − Rayleigh", "#d62728"),
        (rrs, "Rotational Raman contribution", "#2ca02c"),
        (cabannes + rrs - rayleigh,
         "Cabannes + RRS − Rayleigh", "#9467bd"),
    )


def draw_panels(wavelength, panels, suptitle, output, ylabel,
                explicit_scientific=False, oco_overlay=None,
                noise_cloud=None):
    order = np.argsort(wavelength)
    wavelength = wavelength[order]
    fig, axes = plt.subplots(4, 1, figsize=(13, 13), sharex=True,
                             constrained_layout=True)
    overlay_panels = None
    if oco_overlay is not None:
        overlay_wavelength, overlay_panel_data = oco_overlay
        overlay_order = np.argsort(overlay_wavelength)
        overlay_wavelength = np.asarray(overlay_wavelength)[overlay_order]
        overlay_panels = list(overlay_panel_data)

    cloud_panels = None
    if noise_cloud is not None:
        cloud_wavelength, cloud_panel_data = noise_cloud
        cloud_order = np.argsort(cloud_wavelength)
        cloud_wavelength = np.asarray(cloud_wavelength)[cloud_order]
        cloud_panels = list(cloud_panel_data)

    for panel_index, (ax, (y, title, color)) in enumerate(zip(axes, panels)):
        ax.plot(wavelength, y[order], color=color, linewidth=0.9)
        cloud_center = None
        cloud_sigma = None
        if cloud_panels is not None:
            center, sigma, cloud_label = cloud_panels[panel_index]
            cloud_center = np.asarray(center, dtype=float)[cloud_order]
            cloud_sigma = np.asarray(sigma, dtype=float)[cloud_order]
            ax.fill_between(
                cloud_wavelength,
                cloud_center - cloud_sigma,
                cloud_center + cloud_sigma,
                color="0.45",
                alpha=0.28,
                linewidth=0,
                zorder=1,
                label=cloud_label,
            )
        if overlay_panels is not None:
            overlay_y = np.asarray(overlay_panels[panel_index][0])
            overlay_y_sorted = overlay_y[overlay_order]
            ax.plot(
                overlay_wavelength,
                overlay_y_sorted,
                color="black",
                linewidth=1.6,
                zorder=5,
                label=(
                    "OCO analyzer + Gaussian ILS; divided by M₁₁=0.5"
                    if panel_index == 0 else None
                ),
            )
            if cloud_center is not None:
                contrast = np.abs(overlay_y_sorted - cloud_center)
                significance = contrast / cloud_sigma
                if np.any(contrast > 10 * np.finfo(float).eps):
                    ax.text(
                        0.01,
                        0.96,
                        f"max |signal|/σ = {np.max(significance):.2f}; "
                        f"samples >1σ = {100*np.mean(significance > 1):.1f}%",
                        transform=ax.transAxes,
                        fontsize=8,
                        va="top",
                        color="0.2",
                        bbox={
                            "facecolor": "white",
                            "edgecolor": "none",
                            "alpha": 0.72,
                            "pad": 2.0,
                        },
                    )
        ax.axhline(0.0, color="black", linewidth=0.5, alpha=0.35)
        ax.set_title(title, loc="left")
        ax.set_ylabel(ylabel)
        if explicit_scientific:
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))
        ax.grid(alpha=0.22)
        if cloud_panels is not None or (overlay_panels is not None and panel_index == 0):
            ax.legend(loc="best", frameon=False, fontsize=8)
    axes[-1].set_xlabel("Wavelength [nm]")
    fig.suptitle(suptitle, fontsize=15)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(output)


def plot_scene(scene, wavelength, output_root=ROOT, oco_root=OCO_ROOT,
               noise_root=NOISE_ROOT):
    data = read_scene(scene)
    oco = read_oco_scene(data["state"], oco_root)
    noise = read_oco_noise(data["state"], noise_root, "o2a")
    state = data["state"]
    output = Path(output_root) / f"state{state:03d}_o2_components.png"
    draw_panels(
        wavelength,
        component_panels(data),
        f"O₂ A-band truth state {state:03d}: {data['surface']} surface, "
        f"SIF={data['sif_case']}, XCO₂={data['xco2']} ppm, "
        f"aerosol={data['aerosol']}\n"
        f"Geometry: SZA={data['sza']:.3f}°, "
        f"VZA={data['vza']:.3f}°, "
        f"Δφ={data['relative_azimuth']:.1f}°",
        output,
        "Stokes I\n[mW m⁻² sr⁻¹ (cm⁻¹)⁻¹]",
        oco_overlay=(oco["wavelength"], component_panels(oco)),
        noise_cloud=(noise["wavelength"], single_scene_o2_clouds(oco, noise)),
    )


def plot_sif_difference(off, on, wavelength, output_root=ROOT,
                        oco_root=OCO_ROOT, noise_root=NOISE_ROOT):
    difference = {
        name: on[name] - off[name]
        for name in ("rayleigh", "cabannes", "rrs")
    }
    oco_off = read_oco_scene(off["state"], oco_root)
    oco_on = read_oco_scene(on["state"], oco_root)
    noise_on = read_oco_noise(on["state"], noise_root, "o2a")
    oco_difference = {
        name: oco_on[name] - oco_off[name]
        for name in ("rayleigh", "cabannes", "rrs")
    }
    output = Path(output_root) / (
        f"state{on['state']:03d}_minus{off['state']:03d}_"
        "sif_difference_o2_components.png"
    )
    draw_panels(
        wavelength,
        component_panels(difference),
        f"O₂ A-band SIF difference: state {on['state']:03d} − "
        f"{off['state']:03d} (SIF on − off), {on['surface']} surface, "
        f"XCO₂={on['xco2']} ppm, aerosol={on['aerosol']}",
        output,
        "Δ Stokes I\n[mW m⁻² sr⁻¹ (cm⁻¹)⁻¹]",
        explicit_scientific=True,
        oco_overlay=(oco_on["wavelength"], component_panels(oco_difference)),
        noise_cloud=(
            noise_on["wavelength"],
            effect_scene_o2_difference_clouds(noise_on),
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scene-root", type=Path, default=ROOT)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--oco-root", type=Path, default=OCO_ROOT)
    parser.add_argument("--noise-root", type=Path, default=NOISE_ROOT)
    args = parser.parse_args()
    scene_root = args.scene_root.resolve()
    output_root = (args.output_root or scene_root).resolve()
    oco_root = args.oco_root.resolve()
    noise_root = args.noise_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    with Dataset(scene_root / "sim_wavelength.nc") as ds:
        wavelength = np.asarray(ds["o2a_wavelength"][:])

    completed = [scene for scene in sorted(scene_root.glob("hiressim_*.nc"))
                 if o2_is_complete(scene)]
    if not completed:
        raise RuntimeError(f"no completed O2 scene files found in {scene_root}")
    data = [read_scene(scene) for scene in completed]
    for scene in completed:
        plot_scene(scene, wavelength, output_root, oco_root, noise_root)
    print(f"Generated {len(completed)} completed-state plots")

    pairs = {}
    for item in data:
        key = (item["surface"], item["xco2"], item["aerosol"])
        pairs.setdefault(key, {})[item["sif_case"]] = item
    npairs = 0
    for cases in pairs.values():
        if "off" not in cases or len(cases) < 2:
            continue
        on_cases = [item for name, item in cases.items() if name != "off"]
        for on in on_cases:
            plot_sif_difference(
                cases["off"], on, wavelength, output_root, oco_root,
                noise_root,
            )
            npairs += 1
    print(f"Generated {npairs} matched SIF-on-minus-off plots")


if __name__ == "__main__":
    main()
