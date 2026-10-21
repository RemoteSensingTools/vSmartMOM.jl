#!/usr/bin/env python3
"""Plot a terminal retrieval fit and its bandwise residuals versus noise."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


BAND_LABELS = ("O$_2$ A band", "Weak CO$_2$ band", "Strong CO$_2$ band")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("retrieval", type=Path, help="Completed retrieval NetCDF")
    parser.add_argument("--output", type=Path, help="Output PNG path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output or args.retrieval.with_name(
        args.retrieval.stem + "_fit_and_residuals.png"
    )

    with Dataset(args.retrieval) as dataset:
        if int(dataset.getncattr("retrieval_complete")) != 1:
            raise RuntimeError(f"retrieval is not marked complete: {args.retrieval}")
        wavelength = np.asarray(dataset["wavelength"][:], dtype=float)
        measurement = np.asarray(dataset["measurement_perturbed"][:], dtype=float)
        noise = np.asarray(dataset["noise_standard_deviation"][:], dtype=float)
        forward = np.asarray(dataset["final_forward_model"][:], dtype=float)
        residual = np.asarray(dataset["final_residual"][:], dtype=float)
        band_index = np.asarray(dataset["band_index"][:], dtype=int)
        chi_squared = np.asarray(
            dataset["final_band_reduced_chi_squared"][:], dtype=float
        )
        retrieval_class = str(dataset.getncattr("measurement_class"))
        state_index = int(dataset.getncattr("truth_state_index"))
        perturbation_index = int(dataset.getncattr("perturbation_index"))

    if not (
        wavelength.shape == measurement.shape == noise.shape == forward.shape
        == residual.shape == band_index.shape
    ):
        raise RuntimeError("measurement-space arrays do not share one shape")
    if not all(
        np.all(np.isfinite(values))
        for values in (wavelength, measurement, noise, forward, residual)
    ):
        raise RuntimeError("retrieval fit contains non-finite values")
    if not np.all(noise > 0):
        raise RuntimeError("noise standard deviations must be positive")
    if not np.allclose(residual, forward - measurement, rtol=1e-12, atol=1e-12):
        raise RuntimeError("stored final residual is inconsistent with F(x) - y")

    fig, axes = plt.subplots(3, 2, figsize=(15, 10), constrained_layout=True)
    for iband, (axis_row, label) in enumerate(zip(axes, BAND_LABELS), start=1):
        fit_axis, residual_axis = axis_row
        selected = band_index == iband
        fit_axis.fill_between(
            wavelength[selected],
            measurement[selected] - noise[selected],
            measurement[selected] + noise[selected],
            color="0.65", alpha=0.45, linewidth=0,
            label="Measurement noise cloud ($y \\pm 1\\sigma$)", zorder=1,
        )
        fit_axis.plot(
            wavelength[selected], measurement[selected],
            color="black", linewidth=0, marker=".", markersize=2.0,
            alpha=0.65, label="Measurement $y$", zorder=3,
        )
        fit_axis.plot(
            wavelength[selected], forward[selected],
            color="#d62728", linewidth=1.25, label="Terminal $F(x)$", zorder=2,
        )
        fit_axis.set_title(
            f"{label}: fit    reduced $\\chi^2$ = {chi_squared[iband - 1]:.2f}"
        )
        fit_axis.set_ylabel("Radiance\n(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)")

        residual_axis.fill_between(
            wavelength[selected], -noise[selected], noise[selected],
            color="0.65", alpha=0.55, linewidth=0,
            label="$\\pm1\\sigma$ noise cloud", zorder=1,
        )
        residual_axis.axhline(0, color="0.25", linewidth=0.7, zorder=2)
        residual_axis.plot(
            wavelength[selected], residual[selected],
            color="#1f77b4", linewidth=0.8, marker=".", markersize=1.5,
            alpha=0.8, label="Residual $F(x)-y$", zorder=3,
        )
        residual_axis.set_title(f"{label}: residual")
        residual_axis.set_ylabel(
            "$F(x)-y$\n(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)"
        )

        for axis in axis_row:
            axis.grid(alpha=0.2)
            axis.margins(x=0)

    axes[0, 0].legend(loc="best", frameon=False)
    axes[0, 1].legend(loc="best", frameon=False)
    axes[-1, 0].set_xlabel("Wavelength (nm)")
    axes[-1, 1].set_xlabel("Wavelength (nm)")
    fig.suptitle(
        f"{retrieval_class.capitalize()} retrieval: state {state_index:03d}, "
        f"perturbation {perturbation_index:02d}",
        fontsize=14,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
