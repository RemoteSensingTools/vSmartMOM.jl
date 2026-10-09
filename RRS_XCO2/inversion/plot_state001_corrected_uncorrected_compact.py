#!/usr/bin/env python3
"""Compare compact terminal states for one matched truth-state ensemble.

The compact state replaces the twelve retrieved layer CO2 VMRs with the
bottom-layer VMR and derived dry-air-column XCO2 diagnostic.  Every other
retrieval quantity is retained.  Aerosol optical depths and profile-center
heights are transformed from logarithmic retrieval coordinates back to
physical values.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from plot_retrieval_state_convergence import (
    SCENE_COMPONENTS,
    TRUTH_TABLE,
    aerosol_aod760,
    co2_truth_label,
    physical_state,
    table_row,
    truth_values,
)
from retrieval_plot_schema import RetrievalPlotSchema


HERE = Path(__file__).resolve().parent

PARAMETERS = (
    ("bottom_co2", "Bottom-layer CO$_2$ (ppm)"),
    ("XCO2", "XCO$_2$ (ppm)"),
    ("psurf", "Surface pressure (hPa)"),
    ("sulfate_aod760", "Sulphate AOD$_{760}$"),
    ("organic_carbon_aod760", "Organic-carbon AOD$_{760}$"),
    ("utls_sulfate_aod760", "UTLS sulphate AOD$_{760}$"),
    ("sulfate_z0", "Sulphate height (km)"),
    ("organic_carbon_z0", "Organic-carbon height (km)"),
    ("utls_sulfate_z0", "UTLS sulphate height (km)"),
    ("o2a_surface_P0", "O$_2$ A surface $P_0$"),
    ("o2a_surface_P1", "O$_2$ A surface $P_1$"),
    ("o2a_surface_P2", "O$_2$ A surface $P_2$"),
    ("weak_co2_surface_P0", "Weak CO$_2$ surface $P_0$"),
    ("weak_co2_surface_P1", "Weak CO$_2$ surface $P_1$"),
    ("weak_co2_surface_P2", "Weak CO$_2$ surface $P_2$"),
    ("strong_co2_surface_P0", "Strong CO$_2$ surface $P_0$"),
    ("strong_co2_surface_P1", "Strong CO$_2$ surface $P_1$"),
    ("strong_co2_surface_P2", "Strong CO$_2$ surface $P_2$"),
    ("SIF760", "SIF760 (mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)"),
    ("mSIF", "mSIF (native units)"),
)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("state", type=int, nargs="?", default=1)
    parser.add_argument("--inversion-root", type=Path, default=HERE)
    parser.add_argument("--truth-table", type=Path, default=TRUTH_TABLE)
    parser.add_argument(
        "--scene-components", type=Path, default=SCENE_COMPONENTS,
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--table-output", type=Path)
    return parser.parse_args()


def load_class(inversion_root, retrieval_class, state_index):
    records = []
    reference_schema = None
    for path in sorted(
        (inversion_root / retrieval_class).glob(
            "retrieval_state%03d_perturbation*.nc" % state_index
        )
    ):
        with Dataset(path) as dataset:
            if int(dataset.getncattr("retrieval_complete")) != 1:
                continue
            if int(dataset.getncattr("truth_state_index")) != state_index:
                continue
            declared_class = str(dataset.getncattr("measurement_class"))
            if declared_class != retrieval_class:
                raise RuntimeError(
                    "%s declares class %r, expected %r" %
                    (path, declared_class, retrieval_class)
                )
            schema = RetrievalPlotSchema.from_dataset(dataset, path)
            if reference_schema is None:
                reference_schema = schema
            elif schema.identity != reference_schema.identity:
                raise RuntimeError(
                    "%s uses a different retrieval-state schema from %s" %
                    (path, reference_schema.source)
                )
            perturbation = int(dataset.getncattr("perturbation_index"))
            final_state = schema.expand_absolute_state(
                np.asarray(dataset["final_state"][:], dtype=float)
            )
            prior_state = schema.expand_absolute_state(
                np.asarray(dataset["a_priori_state"][:], dtype=float)
            )
            final = physical_state(
                schema.canonical_names, final_state, float(dataset["XCO2"][:])
            )
            prior = physical_state(
                schema.canonical_names, prior_state,
                float(dataset["a_priori_XCO2"][:]),
            )
        records.append((perturbation, final, prior))
    if not records:
        raise RuntimeError(
            "no completed state-%03d %s retrievals" %
            (state_index, retrieval_class)
        )
    return records, reference_schema


def write_table(path, state_index, truth, corrected, uncorrected, schema):
    corrected_by_pert = {pert: state for pert, state, _ in corrected}
    uncorrected_by_pert = {pert: state for pert, state, _ in uncorrected}
    perturbations = sorted(set(corrected_by_pert) & set(uncorrected_by_pert))
    noisy = [p for p in perturbations if p != 11]
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# State %03d compact terminal-state comparison. Layer CO2 is "
            "summarized by bottom-layer VMR and XCO2.\n" % state_index
        )
        stream.write("# Retrieval state schema: %s\n" % schema.description())
        stream.write("# Schema identity: %r\n" % (schema.identity,))
        stream.write(
            "# Standard deviations are sample standard deviations across matched "
            "noise perturbations 01--10; perturbation 11 is excluded.\n"
        )
        stream.write(
            "# parameter truth corrected_mean corrected_sd uncorrected_mean "
            "uncorrected_sd paired_difference_mean paired_difference_sd "
            "corrected_noiseless uncorrected_noiseless\n"
        )
        for key, _ in PARAMETERS:
            if noisy:
                corr = np.array([corrected_by_pert[p][key] for p in noisy])
                uncorr = np.array([uncorrected_by_pert[p][key] for p in noisy])
                delta = uncorr - corr
                ddof = 1 if len(noisy) > 1 else 0
                statistics = (
                    corr.mean(), corr.std(ddof=ddof),
                    uncorr.mean(), uncorr.std(ddof=ddof),
                    delta.mean(), delta.std(ddof=ddof),
                )
            else:
                statistics = (np.nan,) * 6
            corrected_noiseless = (
                corrected_by_pert[11][key] if 11 in corrected_by_pert else np.nan
            )
            uncorrected_noiseless = (
                uncorrected_by_pert[11][key]
                if 11 in uncorrected_by_pert else np.nan
            )
            stream.write(
                f"{key} {truth[key]:.12e} "
                + " ".join("%.12e" % value for value in statistics)
                + " %.12e %.12e" % (
                    corrected_noiseless, uncorrected_noiseless
                )
                + "\n"
            )


def main():
    args = parse_args()
    corrected, corrected_schema = load_class(
        args.inversion_root, "corrected", args.state
    )
    uncorrected, uncorrected_schema = load_class(
        args.inversion_root, "uncorrected", args.state
    )
    if corrected_schema.identity != uncorrected_schema.identity:
        raise RuntimeError(
            "corrected and uncorrected retrievals use different "
            "retrieval-state schemas"
        )
    corrected_by_pert = {pert: state for pert, state, _ in corrected}
    uncorrected_by_pert = {pert: state for pert, state, _ in uncorrected}
    perturbations = sorted(set(corrected_by_pert) & set(uncorrected_by_pert))
    if not perturbations:
        raise RuntimeError("no matched corrected/uncorrected perturbations")
    prior = corrected[0][2]
    truth_row = table_row(args.truth_table, args.state)
    corrected_schema.validate_truth_row(truth_row, str(args.truth_table))
    aerosol_case = truth_row["aerosol_case"]
    truth = truth_values(
        truth_row, aerosol_aod760(args.scene_components, aerosol_case), prior,
        schema=corrected_schema,
    )
    output = args.output or args.inversion_root / (
        "state%03d_corrected_vs_uncorrected_compact_state.png" % args.state
    )
    table_output = args.table_output or output.with_suffix(".dat")

    fig, axes = plt.subplots(6, 4, figsize=(18, 23))
    axes = axes.ravel()
    colors = {"corrected": "#1f77b4", "uncorrected": "#d62728"}
    markers = {"corrected": "o", "uncorrected": "s"}

    for axis, (key, label) in zip(axes, PARAMETERS):
        for retrieval_class, states in (
            ("corrected", corrected_by_pert),
            ("uncorrected", uncorrected_by_pert),
        ):
            values = np.array([states[p][key] for p in perturbations])
            noisy = np.asarray(perturbations) != 11
            if np.any(noisy):
                axis.plot(
                    np.asarray(perturbations)[noisy], values[noisy],
                    color=colors[retrieval_class],
                    marker=markers[retrieval_class], markersize=4.5,
                    linewidth=1.4, label=retrieval_class.capitalize(),
                )
                mean = values[noisy].mean()
                ddof = 1 if np.count_nonzero(noisy) > 1 else 0
                sigma = values[noisy].std(ddof=ddof)
                axis.axhspan(
                    mean - sigma, mean + sigma,
                    color=colors[retrieval_class], alpha=0.07, linewidth=0,
                )
            if np.any(~noisy):
                axis.scatter(
                    np.asarray(perturbations)[~noisy], values[~noisy],
                    color=colors[retrieval_class], marker="s", s=48,
                    edgecolors="white", linewidths=0.7, zorder=4,
                    label=(retrieval_class.capitalize() + " noiseless"),
                )

        if np.isfinite(truth[key]):
            axis.axhline(truth[key], color="black", linewidth=1.4, label="Truth")
        axis.axhline(
            prior[key], color="0.45", linestyle=":", linewidth=1.5, label="Prior"
        )
        axis.set_title(label, fontsize=11)
        axis.set_xticks((1, 3, 5, 7, 9))
        axis.grid(alpha=0.22)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))

    for axis in axes[:len(PARAMETERS)]:
        axis.set_xlabel("Noise perturbation")
    handles, labels = axes[0].get_legend_handles_labels()
    legend_axis = axes[len(PARAMETERS)]
    legend_axis.axis("off")
    legend_axis.legend(handles, labels, frameon=False, loc="center", fontsize=12)
    legend_axis.text(
        0.5,
        0.24,
        "Shading: ensemble mean $\\pm1\\sigma$\n"
        "Aerosol heights are nominal profile centers;\n"
        "they are non-identifiable when truth AOD = 0",
        ha="center",
        va="center",
        fontsize=10,
        transform=legend_axis.transAxes,
    )
    for axis in axes[len(PARAMETERS) + 1:]:
        axis.axis("off")
    fig.suptitle(
        "State %03d: corrected versus uncorrected terminal retrieval states\n"
        "%s; %s\n"
        "Layer CO$_2$ summarized by bottom-layer VMR and XCO$_2$" %
        (
            args.state, co2_truth_label(truth_row),
            corrected_schema.description(),
        ),
        fontsize=16,
        y=0.985,
    )
    fig.subplots_adjust(
        left=0.06, right=0.985, bottom=0.045, top=0.935, hspace=0.52, wspace=0.30
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    table_output.parent.mkdir(parents=True, exist_ok=True)
    write_table(
        table_output, args.state, truth, corrected, uncorrected,
        corrected_schema,
    )
    print(output)
    print(table_output)


if __name__ == "__main__":
    main()
