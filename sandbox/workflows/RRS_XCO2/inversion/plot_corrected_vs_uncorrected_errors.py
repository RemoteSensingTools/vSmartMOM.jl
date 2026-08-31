#!/usr/bin/env python3
"""Scatter corrected versus uncorrected retrieval errors for one truth state.

Every marker represents one matched noise perturbation.  The horizontal axis is
``x_corrected - x_truth`` and the vertical axis is
``x_uncorrected - x_truth``.  Layer-resolved CO2 is replaced by the saved
dry-air-column XCO2 diagnostic, and logarithmic aerosol retrieval coordinates
are converted back to physical AOD and height before differencing.

The 3-by-7 panel placement is intentionally fixed to support direct comparison
between truth states.  The two unused cells in columns 1 and 2 are left blank.
"""

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from netCDF4 import Dataset

from plot_retrieval_state_convergence import (
    SCENE_COMPONENTS,
    TRUTH_TABLE,
    aerosol_aod760,
    physical_state,
    table_row,
    truth_values,
)


HERE = Path(__file__).resolve().parent
WORKFLOW_ROOT = HERE.parent
DATA_ROOT = Path(os.environ.get("RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
INVERSION_ROOT = DATA_ROOT / "inversion"

# key, panel title, row, column, physical plotting unit
PANELS = (
    ("XCO2", "XCO$_2$", 0, 0, "ppm"),
    ("psurf", "$p_s$", 1, 0, "hPa"),
    ("SIF760", "SIF$_{760}$ (reference)", 0, 1,
     "mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$"),
    ("mSIF", "$m_{SIF}$", 1, 1, "native"),
    ("sulfate_aod760", "Sulphate AOD$_{760}$", 0, 2, "1"),
    ("organic_carbon_aod760", "Organic-carbon AOD$_{760}$", 1, 2, "1"),
    ("utls_sulfate_aod760", "UTLS sulphate AOD$_{760}$", 2, 2, "1"),
    ("sulfate_z0", "Sulphate $z_0$", 0, 3, "km"),
    ("organic_carbon_z0", "Organic-carbon $z_0$", 1, 3, "km"),
    ("utls_sulfate_z0", "UTLS sulphate $z_0$", 2, 3, "km"),
    ("o2a_surface_P0", "$\\rho_0$: O$_2$ A", 0, 4, "1"),
    ("weak_co2_surface_P0", "$\\rho_0$: weak CO$_2$", 1, 4, "1"),
    ("strong_co2_surface_P0", "$\\rho_0$: strong CO$_2$", 2, 4, "1"),
    ("o2a_surface_P1", "$\\rho_1$: O$_2$ A", 0, 5, "1"),
    ("weak_co2_surface_P1", "$\\rho_1$: weak CO$_2$", 1, 5, "1"),
    ("strong_co2_surface_P1", "$\\rho_1$: strong CO$_2$", 2, 5, "1"),
    ("o2a_surface_P2", "$\\rho_2$: O$_2$ A", 0, 6, "1"),
    ("weak_co2_surface_P2", "$\\rho_2$: weak CO$_2$", 1, 6, "1"),
    ("strong_co2_surface_P2", "$\\rho_2$: strong CO$_2$", 2, 6, "1"),
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "state", type=int, help="Truth-state index, for example 1 for state001"
    )
    parser.add_argument(
        "--inversion-root", type=Path, default=INVERSION_ROOT,
        help="Directory containing corrected/ and uncorrected/",
    )
    parser.add_argument("--output", type=Path, help="Output PNG")
    parser.add_argument(
        "--table-output", type=Path,
        help="Whitespace table containing every plotted value",
    )
    return parser.parse_args()


def completed_retrievals(directory, state_index, expected_class):
    records = {}
    pattern = f"retrieval_state{state_index:03d}_perturbation*.nc"
    for path in sorted(directory.glob(pattern)):
        try:
            with Dataset(path) as dataset:
                if int(dataset.getncattr("retrieval_complete")) != 1:
                    continue
                if int(dataset.getncattr("truth_state_index")) != state_index:
                    continue
                retrieval_class = str(dataset.getncattr("measurement_class"))
                if retrieval_class != expected_class:
                    raise RuntimeError(
                        f"{path} declares class {retrieval_class!r}, expected "
                        f"{expected_class!r}"
                    )
                perturbation = int(dataset.getncattr("perturbation_index"))
                names = str(dataset.getncattr("parameter_names")).split()
                final_state = np.asarray(dataset["final_state"][:], dtype=float)
                prior_state = np.asarray(
                    dataset["a_priori_state"][:], dtype=float
                )
                final_xco2 = float(dataset["XCO2"][:])
                prior_xco2 = float(dataset["a_priori_XCO2"][:])
                metadata = {
                    "surface": str(dataset.getncattr("surface")),
                    "aerosol_case": str(dataset.getncattr("aerosol_case")),
                    "truth_xco2": float(dataset.getncattr("truth_xco2_ppm")),
                }
        except (OSError, RuntimeError, AttributeError) as error:
            raise RuntimeError(f"could not read completed retrieval {path}") from error
        if perturbation in records:
            raise RuntimeError(
                f"duplicate {expected_class} perturbation {perturbation:02d}"
            )
        records[perturbation] = {
            "path": path,
            "final": physical_state(names, final_state, final_xco2),
            "prior": physical_state(names, prior_state, prior_xco2),
            "metadata": metadata,
        }
    return records


def matched_records(inversion_root, state_index):
    corrected = completed_retrievals(
        inversion_root / "corrected", state_index, "corrected"
    )
    uncorrected = completed_retrievals(
        inversion_root / "uncorrected", state_index, "uncorrected"
    )
    perturbations = sorted(set(corrected) & set(uncorrected))
    if not perturbations:
        raise RuntimeError(
            f"state {state_index:03d} has no matched complete retrievals"
        )
    return corrected, uncorrected, perturbations


def check_pair_metadata(corrected, uncorrected, perturbations):
    reference = corrected[perturbations[0]]["metadata"]
    for perturbation in perturbations:
        for retrieval_class, records in (
            ("corrected", corrected), ("uncorrected", uncorrected)
        ):
            if records[perturbation]["metadata"] != reference:
                raise RuntimeError(
                    f"state metadata changed in {retrieval_class} perturbation "
                    f"{perturbation:02d}"
                )
    return reference


def write_table(path, perturbations, truth, corrected, uncorrected):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# Paired corrected and uncorrected terminal retrieval errors.\n"
        )
        stream.write(
            "# parameter units perturbation truth corrected uncorrected "
            "corrected_minus_truth uncorrected_minus_truth\n"
        )
        for key, _, _, _, unit in PANELS:
            for perturbation in perturbations:
                corrected_value = corrected[perturbation]["final"][key]
                uncorrected_value = uncorrected[perturbation]["final"][key]
                stream.write(
                    f"{key} {unit.replace(' ', '_')} {perturbation:02d} "
                    f"{truth[key]:.12g} {corrected_value:.12g} "
                    f"{uncorrected_value:.12g} "
                    f"{corrected_value - truth[key]:.12g} "
                    f"{uncorrected_value - truth[key]:.12g}\n"
                )


def symmetric_limit(x_values, y_values):
    maximum = float(
        np.max(np.abs(np.concatenate([x_values, y_values, np.asarray([0.0])])))
    )
    if maximum == 0.0:
        maximum = 1.0
    return 1.16 * maximum


def main():
    args = parse_args()
    corrected, uncorrected, perturbations = matched_records(
        args.inversion_root, args.state
    )
    metadata = check_pair_metadata(corrected, uncorrected, perturbations)
    truth_row = table_row(TRUTH_TABLE, args.state)
    truth = truth_values(
        truth_row,
        aerosol_aod760(SCENE_COMPONENTS, metadata["aerosol_case"]),
        corrected[perturbations[0]]["prior"],
    )

    output = args.output or args.inversion_root / (
        f"state{args.state:03d}_corrected_vs_uncorrected_parameter_errors.png"
    )
    table_output = args.table_output or output.with_suffix(".dat")
    write_table(table_output, perturbations, truth, corrected, uncorrected)

    colors = plt.get_cmap("tab10")(np.arange(len(perturbations)) % 10)
    figure, axes = plt.subplots(3, 7, figsize=(23.5, 10.5))
    for axis in axes.ravel():
        axis.set_visible(False)

    for key, title, row, column, unit in PANELS:
        axis = axes[row, column]
        axis.set_visible(True)
        x_values = np.asarray([
            corrected[perturbation]["final"][key] - truth[key]
            for perturbation in perturbations
        ])
        y_values = np.asarray([
            uncorrected[perturbation]["final"][key] - truth[key]
            for perturbation in perturbations
        ])
        limit = symmetric_limit(x_values, y_values)
        axis.plot(
            [-limit, limit], [-limit, limit], color="0.35", linestyle="--",
            linewidth=1.0, zorder=1,
        )
        axis.axhline(0.0, color="0.78", linewidth=0.8, zorder=0)
        axis.axvline(0.0, color="0.78", linewidth=0.8, zorder=0)
        axis.scatter(
            x_values, y_values, c=colors, s=49, edgecolors="white",
            linewidths=0.7, zorder=3,
        )
        axis.set_xlim(-limit, limit)
        axis.set_ylim(-limit, limit)
        axis.set_aspect("equal", adjustable="box")
        unit_title = "" if unit == "1" else f" [{unit}]"
        axis.set_title(title + unit_title, fontsize=10.5, pad=5)
        axis.ticklabel_format(
            axis="both", style="sci", scilimits=(-3, 3), useMathText=True
        )
        axis.tick_params(labelsize=8)
        axis.grid(alpha=0.16)

    legend_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=7,
            markerfacecolor=color, markeredgecolor="white",
            label=f"Perturbation {perturbation:02d}",
        )
        for perturbation, color in zip(perturbations, colors)
    ]
    figure.legend(
        handles=legend_handles, loc="lower center", ncol=min(10, len(legend_handles)),
        frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.005),
    )
    figure.text(
        0.5, 0.055,
        "Corrected retrieval error: $x_{corr}-x_{truth}$ (panel units)",
        ha="center", va="center", fontsize=13,
    )
    figure.text(
        0.018, 0.49,
        "Uncorrected retrieval error: $x_{uncorr}-x_{truth}$ (panel units)",
        ha="center", va="center", rotation="vertical", fontsize=13,
    )
    figure.suptitle(
        f"State {args.state:03d}: corrected versus uncorrected terminal-state "
        f"errors\n{metadata['surface']}, {metadata['aerosol_case']}, "
        f"truth XCO$_2$={metadata['truth_xco2']:.0f} ppm; "
        f"{len(perturbations)} matched perturbations",
        fontsize=16, y=0.975,
    )
    figure.subplots_adjust(
        left=0.06, right=0.985, bottom=0.12, top=0.87,
        hspace=0.40, wspace=0.43,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    print(output)
    print(table_output)
    print("matched perturbations:", " ".join(f"{p:02d}" for p in perturbations))


if __name__ == "__main__":
    main()
