#!/usr/bin/env python3
"""Plot truth displacements for every completed corrected/uncorrected pair.

The horizontal coordinate is ``x_corrected - x_truth`` and the vertical
coordinate is ``x_uncorrected - x_truth``.  Layer-resolved CO2 is represented
by both the retrieved bottom-layer VMR and saved dry-air-column XCO2
diagnostic, while logarithmic aerosol coordinates are transformed back to
physical AOD and height.

Only retrievals with ``retrieval_complete == 1`` in both classes are paired.
Noise perturbations 01--10 are circles and the noiseless retrieval 11 is a
square.  Colors identify truth states.  This aggregate diagnostic deliberately
does not draw the optional gain/noise overlay used by the single-state plot.
"""

import argparse
from datetime import datetime
from pathlib import Path
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from netCDF4 import Dataset

from co2_plot_metadata import bottom_layer_mapping_from_truth_table
from plot_corrected_vs_uncorrected_errors import PANELS, symmetric_limit
from plot_retrieval_state_convergence import (
    SCENE_COMPONENTS,
    TRUTH_TABLE,
    aerosol_aod760,
    physical_state,
    table_row,
    truth_values,
)
from retrieval_plot_schema import RetrievalPlotSchema


HERE = Path(__file__).resolve().parent
RETRIEVAL_NAME = re.compile(
    r"^retrieval_state([0-9]+)_perturbation([0-9]+)\.nc$"
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inversion-root", type=Path, default=HERE,
        help="Directory containing corrected/ and uncorrected/",
    )
    parser.add_argument(
        "--output", type=Path,
        help=(
            "Output PNG (default: all_completed_corrected_vs_uncorrected_"
            "truth_displacements.png under the inversion root)"
        ),
    )
    parser.add_argument(
        "--table-output", type=Path,
        help="Whitespace table containing every plotted value",
    )
    parser.add_argument(
        "--truth-table", type=Path, default=TRUTH_TABLE,
        help="Campaign true_states.dat (default: full-column truth map)",
    )
    parser.add_argument(
        "--scene-components", type=Path, default=SCENE_COMPONENTS,
        help="Campaign scene_components.dat (default: full-column truth map)",
    )
    return parser.parse_args()


def candidate_states(inversion_root):
    """Return the frozen set of state indices visible at scan time."""
    states = set()
    for retrieval_class in ("corrected", "uncorrected"):
        directory = inversion_root / retrieval_class
        for path in directory.glob("retrieval_state*_perturbation*.nc"):
            match = RETRIEVAL_NAME.match(path.name)
            if match is not None:
                states.add(int(match.group(1)))
    return sorted(states)


def completed_retrievals(directory, state_index, expected_class):
    """Load complete files, while safely ignoring files still being written."""
    records = {}
    pattern = "retrieval_state%03d_perturbation*.nc" % state_index
    for path in sorted(directory.glob(pattern)):
        try:
            dataset = Dataset(str(path))
        except OSError as error:
            warnings.warn("skipping unreadable in-progress file %s: %s" % (
                path, error,
            ))
            continue

        with dataset:
            try:
                complete = int(dataset.getncattr("retrieval_complete"))
            except AttributeError:
                warnings.warn(
                    "skipping file without retrieval_complete attribute: %s"
                    % path
                )
                continue
            if complete != 1:
                continue

            try:
                file_state = int(dataset.getncattr("truth_state_index"))
                retrieval_class = str(
                    dataset.getncattr("measurement_class")
                )
                perturbation = int(dataset.getncattr("perturbation_index"))
                schema = RetrievalPlotSchema.from_dataset(dataset, path)
                final_state = schema.expand_absolute_state(
                    np.asarray(dataset["final_state"][:], dtype=float)
                )
                prior_state = schema.expand_absolute_state(
                    np.asarray(dataset["a_priori_state"][:], dtype=float)
                )
                final_xco2 = float(dataset["XCO2"][:])
                prior_xco2 = float(dataset["a_priori_XCO2"][:])
                metadata = {
                    "surface": str(dataset.getncattr("surface")),
                    "aerosol_case": str(
                        dataset.getncattr("aerosol_case")
                    ),
                    "truth_xco2": float(
                        dataset.getncattr("truth_xco2_ppm")
                    ),
                    "truth_bottom_co2": (
                        float(dataset.getncattr("truth_bottom_co2_ppm"))
                        if "truth_bottom_co2_ppm" in dataset.ncattrs() else None
                    ),
                    "schema_identity": schema.identity,
                    "schema_description": schema.description(),
                    "state_model": schema.state_model,
                    "sif_mode": schema.sif_mode,
                    "active_state_count": len(schema.active_names),
                }
                converged = bool(dataset.getncattr("converged"))
            except (AttributeError, IndexError, KeyError, ValueError) as error:
                raise RuntimeError(
                    "completed retrieval is malformed: %s" % path
                ) from error

        if file_state != state_index:
            raise RuntimeError(
                "%s declares truth state %03d, expected %03d"
                % (path, file_state, state_index)
            )
        if retrieval_class != expected_class:
            raise RuntimeError(
                "%s declares class %r, expected %r"
                % (path, retrieval_class, expected_class)
            )
        if perturbation in records:
            raise RuntimeError(
                "duplicate %s state %03d perturbation %02d"
                % (expected_class, state_index, perturbation)
            )
        records[perturbation] = {
            "path": path,
            "final": physical_state(
                schema.canonical_names, final_state, final_xco2
            ),
            "prior": physical_state(
                schema.canonical_names, prior_state, prior_xco2
            ),
            "metadata": metadata,
            "converged": converged,
            "schema": schema,
        }
    return records


def check_pair_metadata(corrected, uncorrected, perturbations, state_index):
    reference = corrected[perturbations[0]]["metadata"]
    for perturbation in perturbations:
        for retrieval_class, records in (
            ("corrected", corrected), ("uncorrected", uncorrected)
        ):
            if records[perturbation]["metadata"] != reference:
                raise RuntimeError(
                    "state %03d metadata changed in %s perturbation %02d"
                    % (state_index, retrieval_class, perturbation)
                )
    return reference


def collect_pairs(inversion_root, states, truth_table=TRUTH_TABLE,
                  scene_components=SCENE_COMPONENTS):
    """Freeze and flatten all matched, completed class pairs."""
    pairs = []
    unmatched = []
    for state_index in states:
        corrected = completed_retrievals(
            inversion_root / "corrected", state_index, "corrected"
        )
        uncorrected = completed_retrievals(
            inversion_root / "uncorrected", state_index, "uncorrected"
        )
        perturbations = sorted(set(corrected) & set(uncorrected))
        corrected_only = sorted(set(corrected) - set(uncorrected))
        uncorrected_only = sorted(set(uncorrected) - set(corrected))
        if corrected_only or uncorrected_only:
            unmatched.append((state_index, corrected_only, uncorrected_only))
        if not perturbations:
            continue

        metadata = check_pair_metadata(
            corrected, uncorrected, perturbations, state_index
        )
        truth_row = table_row(truth_table, state_index)
        corrected[perturbations[0]]["schema"].validate_truth_row(
            truth_row, str(truth_table)
        )
        truth = truth_values(
            truth_row,
            aerosol_aod760(
                scene_components, metadata["aerosol_case"]
            ),
            corrected[perturbations[0]]["prior"],
            schema=corrected[perturbations[0]]["schema"],
        )
        for perturbation in perturbations:
            pairs.append({
                "state": state_index,
                "perturbation": perturbation,
                "noise_case": (
                    "noiseless" if perturbation == 11 else "perturbed"
                ),
                "truth": truth,
                "corrected": corrected[perturbation],
                "uncorrected": uncorrected[perturbation],
                "metadata": metadata,
            })
    return pairs, unmatched


def write_table(path, pairs, snapshot_time):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# Unified corrected-versus-uncorrected truth displacements for "
            "matched completed retrievals.\n"
        )
        stream.write("# Snapshot %s\n" % snapshot_time)
        schema_descriptions = sorted(set(
            pair["metadata"]["schema_description"] for pair in pairs
        ))
        stream.write(
            "# Retrieval state schema(s): %s\n" %
            " | ".join(schema_descriptions)
        )
        stream.write(
            "# state perturbation noise_case surface aerosol_case "
            "truth_bottom_co2_ppm truth_xco2_ppm "
            "state_model sif_mode active_state_count "
            "parameter units truth corrected_retrieved uncorrected_retrieved "
            "corrected_minus_truth uncorrected_minus_truth corrected_converged "
            "uncorrected_converged corrected_file uncorrected_file\n"
        )
        for pair in pairs:
            for key, _, _, _, unit in PANELS:
                truth = pair["truth"][key]
                corrected = pair["corrected"]["final"][key]
                uncorrected = pair["uncorrected"]["final"][key]
                stream.write(
                    "%03d %02d %s %s %s %s %.12g %s %s %d %s %s %.12g %.12g %.12g "
                    "%.12g %.12g %d %d %s %s\n" % (
                        pair["state"],
                        pair["perturbation"],
                        pair["noise_case"],
                        pair["metadata"]["surface"],
                        pair["metadata"]["aerosol_case"],
                        (
                            "NA" if pair["metadata"]["truth_bottom_co2"] is None
                            else "%.12g" % pair["metadata"]["truth_bottom_co2"]
                        ),
                        pair["metadata"]["truth_xco2"],
                        pair["metadata"]["state_model"],
                        pair["metadata"]["sif_mode"],
                        pair["metadata"]["active_state_count"],
                        key,
                        unit.replace(" ", "_"),
                        truth,
                        corrected,
                        uncorrected,
                        corrected - truth,
                        uncorrected - truth,
                        int(pair["corrected"]["converged"]),
                        int(pair["uncorrected"]["converged"]),
                        pair["corrected"]["path"],
                        pair["uncorrected"]["path"],
                    )
                )


def state_labels(pairs, states):
    labels = {}
    for state_index in states:
        pair = next(pair for pair in pairs if pair["state"] == state_index)
        metadata = pair["metadata"]
        aerosol = metadata["aerosol_case"].replace("_", " ")
        if metadata["truth_bottom_co2"] is None:
            co2_label = "XCO2 %.0f ppm" % metadata["truth_xco2"]
        else:
            co2_label = "bottom %.0f ppm (XCO2 %.3f)" % (
                metadata["truth_bottom_co2"], metadata["truth_xco2"]
            )
        labels[state_index] = "%03d: %s, %s, %s" % (
            state_index, metadata["surface"], aerosol, co2_label,
        )
    return labels


def make_plot(output, pairs, snapshot_time, mapping_note=None):
    states = sorted(set(pair["state"] for pair in pairs))
    if len(states) <= 10:
        color_values = plt.get_cmap("tab10")(
            np.linspace(0.0, 1.0, 10)
        )[:len(states)]
    else:
        color_values = plt.get_cmap("tab20")(
            np.linspace(0.0, 1.0, max(20, len(states)))
        )[:len(states)]
    colors = dict(zip(states, color_values))

    figure, axes = plt.subplots(3, 7, figsize=(23.5, 11.5))
    for axis in axes.ravel():
        axis.set_visible(False)

    for key, title, row, column, unit in PANELS:
        axis = axes[row, column]
        axis.set_visible(True)
        x_values = np.asarray([
            pair["corrected"]["final"][key] - pair["truth"][key]
            for pair in pairs
        ])
        y_values = np.asarray([
            pair["uncorrected"]["final"][key] - pair["truth"][key]
            for pair in pairs
        ])
        limit = symmetric_limit(x_values, y_values)
        axis.plot(
            [-limit, limit], [-limit, limit], color="0.35", linestyle="--",
            linewidth=1.0, zorder=1,
        )
        axis.axhline(0.0, color="0.78", linewidth=0.8, zorder=0)
        axis.axvline(0.0, color="0.78", linewidth=0.8, zorder=0)

        for state_index in states:
            state_pairs = [
                pair for pair in pairs if pair["state"] == state_index
            ]
            noisy = [
                pair for pair in state_pairs if pair["perturbation"] != 11
            ]
            noiseless = [
                pair for pair in state_pairs if pair["perturbation"] == 11
            ]
            if noisy:
                axis.scatter(
                    [pair["corrected"]["final"][key] - pair["truth"][key]
                     for pair in noisy],
                    [pair["uncorrected"]["final"][key] - pair["truth"][key]
                     for pair in noisy],
                    color=[colors[state_index]], marker="o", s=39, alpha=0.76,
                    edgecolors="white", linewidths=0.45, zorder=3,
                )
            if noiseless:
                axis.scatter(
                    [pair["corrected"]["final"][key] - pair["truth"][key]
                     for pair in noiseless],
                    [pair["uncorrected"]["final"][key] - pair["truth"][key]
                     for pair in noiseless],
                    color=[colors[state_index]], marker="s", s=72,
                    edgecolors="black", linewidths=0.85, zorder=5,
                )

        axis.set_xlim(-limit, limit)
        axis.set_ylim(-limit, limit)
        axis.set_aspect("equal", adjustable="box")
        unit_title = "" if unit == "1" else " [%s]" % unit
        axis.set_title(title + unit_title, fontsize=10.5, pad=5)
        axis.ticklabel_format(
            axis="both", style="sci", scilimits=(-3, 3), useMathText=True
        )
        axis.tick_params(labelsize=8)
        axis.grid(alpha=0.16)

    labels = state_labels(pairs, states)
    legend_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=7,
            markerfacecolor="0.35", markeredgecolor="white",
            label="Perturbed measurement (01--10)",
        ),
        Line2D(
            [0], [0], marker="s", linestyle="none", markersize=7,
            markerfacecolor="0.35", markeredgecolor="black",
            label="Noiseless measurement (11)",
        ),
    ]
    legend_handles.extend([
        Patch(facecolor=colors[state_index], edgecolor="none",
              label=labels[state_index])
        for state_index in states
    ])
    legend_columns = 4 if len(legend_handles) <= 16 else 6
    figure.legend(
        handles=legend_handles, loc="lower center", ncol=legend_columns,
        frameon=False, fontsize=8.2, bbox_to_anchor=(0.5, 0.008),
        handlelength=1.25, columnspacing=1.5,
    )

    noisy_count = sum(pair["perturbation"] != 11 for pair in pairs)
    noiseless_count = len(pairs) - noisy_count
    schema_descriptions = sorted(set(
        pair["metadata"]["schema_description"] for pair in pairs
    ))
    schema_summary = " | ".join(schema_descriptions)
    figure.text(
        0.5, 0.135,
        "Corrected truth displacement: $x_{corr}-x_{truth}$ (panel units)",
        ha="center", va="center", fontsize=13,
    )
    figure.text(
        0.018, 0.51,
        "Uncorrected truth displacement: $x_{uncorr}-x_{truth}$ (panel units)",
        ha="center", va="center", rotation="vertical", fontsize=13,
    )
    figure.suptitle(
        "All matched completed retrievals: corrected versus uncorrected "
        "truth displacement\n"
        "%d states; %d perturbed pairs; %d noiseless pairs; snapshot %s\n%s" % (
            len(states), noisy_count, noiseless_count, snapshot_time,
            schema_summary,
        ),
        fontsize=16, y=0.976,
    )
    if mapping_note:
        figure.text(
            0.5, 0.895, mapping_note, ha="center", va="center",
            fontsize=9.0, color="0.30",
        )
    figure.subplots_adjust(
        left=0.06, right=0.985, bottom=0.22,
        top=0.85 if mapping_note else 0.87,
        hspace=0.40, wspace=0.43,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(str(output), dpi=180, bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    output = args.output or args.inversion_root / (
        "all_completed_corrected_vs_uncorrected_truth_displacements.png"
    )
    table_output = args.table_output or output.with_suffix(".dat")
    snapshot_time = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")

    states = candidate_states(args.inversion_root)
    pairs, unmatched = collect_pairs(
        args.inversion_root, states, args.truth_table, args.scene_components
    )
    if not pairs:
        raise RuntimeError("no matched completed corrected/uncorrected pairs")

    write_table(table_output, pairs, snapshot_time)
    mapping_note = bottom_layer_mapping_from_truth_table(args.truth_table)
    make_plot(output, pairs, snapshot_time, mapping_note=mapping_note)

    paired_states = sorted(set(pair["state"] for pair in pairs))
    noisy_count = sum(pair["perturbation"] != 11 for pair in pairs)
    noiseless_count = len(pairs) - noisy_count
    print(output)
    print(table_output)
    print("paired states:", " ".join("%03d" % state for state in paired_states))
    print("matched pairs: %d (%d perturbed, %d noiseless)" % (
        len(pairs), noisy_count, noiseless_count,
    ))
    for state_index, corrected_only, uncorrected_only in unmatched:
        print(
            "excluded unmatched state %03d: corrected-only [%s], "
            "uncorrected-only [%s]" % (
                state_index,
                " ".join("%02d" % value for value in corrected_only),
                " ".join("%02d" % value for value in uncorrected_only),
            )
        )


if __name__ == "__main__":
    main()
