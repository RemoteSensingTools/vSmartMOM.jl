#!/usr/bin/env python3
"""Plot O2 A-band components for every completed aerosol truth state."""

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_truth_o2_components")

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from netCDF4 import Dataset
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "truth_map" / "aerosol_chunked"
WAVELENGTHS = ROOT / "sim_wavelength.nc"


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
        }


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
                explicit_scientific=False):
    order = np.argsort(wavelength)
    wavelength = wavelength[order]
    fig, axes = plt.subplots(4, 1, figsize=(13, 13), sharex=True,
                             constrained_layout=True)
    for ax, (y, title, color) in zip(axes, panels):
        ax.plot(wavelength, y[order], color=color, linewidth=0.9)
        ax.axhline(0.0, color="black", linewidth=0.5, alpha=0.35)
        ax.set_title(title, loc="left")
        ax.set_ylabel(ylabel)
        if explicit_scientific:
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))
        ax.grid(alpha=0.22)
    axes[-1].set_xlabel("Wavelength [nm]")
    fig.suptitle(suptitle, fontsize=15)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(output)


def plot_scene(scene, wavelength):
    data = read_scene(scene)
    state = data["state"]
    output = ROOT / f"state{state:03d}_o2_components.png"
    draw_panels(
        wavelength,
        component_panels(data),
        f"O₂ A-band truth state {state:03d}: {data['surface']} surface, "
        f"SIF={data['sif_case']}, XCO₂={data['xco2']} ppm, "
        f"aerosol={data['aerosol']}",
        output,
        "Stokes I\n[mW m⁻² sr⁻¹ (cm⁻¹)⁻¹]",
    )


def plot_sif_difference(off, on, wavelength):
    difference = {
        name: on[name] - off[name]
        for name in ("rayleigh", "cabannes", "rrs")
    }
    output = ROOT / (
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
    )


def main():
    with Dataset(WAVELENGTHS) as ds:
        wavelength = np.asarray(ds["o2a_wavelength"][:])

    completed = [scene for scene in sorted(ROOT.glob("hiressim_*.nc"))
                 if o2_is_complete(scene)]
    if not completed:
        raise RuntimeError(f"no completed O2 scene files found in {ROOT}")
    data = [read_scene(scene) for scene in completed]
    for scene in completed:
        plot_scene(scene, wavelength)
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
            plot_sif_difference(cases["off"], on, wavelength)
            npairs += 1
    print(f"Generated {npairs} matched SIF-on-minus-off plots")


if __name__ == "__main__":
    main()
