#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ALBEDOS = [0.0, 0.3, 1.0]


def label_float(x: float) -> str:
    return f"{x:.1f}".replace(".", "p")


def read_pair(outdir: Path, version: str, albedo: float):
    tag = label_float(albedo)
    rrs = np.loadtxt(outdir / f"{version}_flat_noabs_alb{tag}_rrs.dat")
    nors = np.loadtxt(outdir / f"{version}_flat_noabs_alb{tag}_nors.dat")
    wn = rrs[:, 0]
    if not np.allclose(wn, nors[:, 0]):
        raise ValueError(f"wavenumber mismatch for {version} albedo {albedo}")
    f0 = rrs[:, 7]
    red = (rrs[:, 1] - nors[:, 1]) / f0
    green = rrs[:, 4] / f0
    blue = (rrs[:, 1] + rrs[:, 4] - nors[:, 1]) / f0
    return wn, red, green, blue


def read_meta(outdir: Path, version: str, albedo: float = 0.0) -> dict[str, str]:
    path = outdir / f"{version}_flat_noabs_alb{label_float(albedo)}_metadata.txt"
    meta: dict[str, str] = {}
    if not path.exists():
        return meta
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        meta[key.strip()] = value.strip()
    return meta


def fmt_meta_value(meta: dict[str, str], *keys: str) -> str:
    for key in keys:
        if key in meta:
            try:
                return f"{float(meta[key]):.8f}"
            except ValueError:
                return meta[key]
    return "n/a"


def column_title(outdir: Path, version: str) -> str:
    meta = read_meta(outdir, version, 0.0)
    if version == "current":
        return (
            "current\n"
            + r"$\rho_\mathrm{Cab}$="
            + fmt_meta_value(meta, "rho_cab_internal")
            + ", "
            + r"$\rho_\mathrm{Rayl}$="
            + fmt_meta_value(meta, "rho_ray_internal")
        )
    return (
        "June 2024\n"
        + r"$\rho_\mathrm{Rayl}$="
        + fmt_meta_value(meta, "rho_ray_model_greek", "rho_ray_yaml", "rho_ray_internal")
        + ", "
        + r"$\rho_\mathrm{Cab}$="
        + fmt_meta_value(meta, "rho_cab_model_greek", "rho_cab_internal")
    )


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: 06_plot_flat_noabs_comparison.py OUTDIR", file=sys.stderr)
        return 2

    outdir = Path(sys.argv[1]).expanduser().resolve()
    fig, axes = plt.subplots(
        len(ALBEDOS), 3, figsize=(14.5, 8.2), sharex=True, constrained_layout=True
    )
    colors = {"red": "#c43c35", "green": "#1b8a5a", "blue": "#2369bd"}

    for row, albedo in enumerate(ALBEDOS):
        cur = read_pair(outdir, "current", albedo)
        old = read_pair(outdir, "june2024", albedo)
        wn = cur[0]
        if not np.allclose(wn, old[0]):
            raise ValueError(f"current/june2024 grid mismatch for albedo {albedo}")

        for col, (title, data) in enumerate(
            [
                (column_title(outdir, "current"), cur),
                (column_title(outdir, "june2024"), old),
                ("current - June 2024", (wn, cur[1] - old[1], cur[2] - old[2], cur[3] - old[3])),
            ]
        ):
            ax = axes[row, col]
            ax.plot(wn, data[1], color=colors["red"], lw=1.1, label=r"$I_\mathrm{Cab}-I_\mathrm{Rayl}$")
            ax.plot(wn, data[2], color=colors["green"], lw=1.1, label=r"$I_\mathrm{RRS}$")
            ax.plot(wn, data[3], color=colors["blue"], lw=1.1, label=r"$\Delta I$")
            ax.axhline(0.0, color="0.35", lw=0.6, alpha=0.55)
            if row == 0:
                ax.set_title(title)
            if col == 0:
                ax.set_ylabel(rf"$\rho={albedo:g}$" + "\n" + r"radiance / $F_0$")
            if row == len(ALBEDOS) - 1:
                ax.set_xlabel(r"wavenumber [cm$^{-1}$]")
            ax.grid(True, color="0.88", lw=0.5)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Flat Sun, No Absorption, 761-764 nm Output Window", y=1.03)

    pdf = outdir / "flat_noabs_761_764_current_vs_june2024_rgb.pdf"
    png = outdir / "flat_noabs_761_764_current_vs_june2024_rgb.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=180)

    summary = outdir / "flat_noabs_761_764_current_minus_june2024_summary.txt"
    with summary.open("w") as f:
        for albedo in ALBEDOS:
            cur = read_pair(outdir, "current", albedo)
            old = read_pair(outdir, "june2024", albedo)
            for name, idx in [("red", 1), ("green", 2), ("blue", 3)]:
                diff = cur[idx] - old[idx]
                f.write(
                    f"albedo={albedo:g} curve={name} "
                    f"max_abs_diff={np.max(np.abs(diff)):.12e} "
                    f"mean_abs_diff={np.mean(np.abs(diff)):.12e} "
                    f"min_diff={np.min(diff):.12e} max_diff={np.max(diff):.12e}\n"
                )
    print(pdf)
    print(png)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
