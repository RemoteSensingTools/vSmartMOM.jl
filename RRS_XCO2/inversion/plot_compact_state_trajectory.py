#!/usr/bin/env python3
"""Compare compact trajectories using bottom-layer CO2 and column XCO2."""

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
    table_row,
    truth_values as build_truth_values,
)
from retrieval_plot_schema import RetrievalPlotSchema


PARAMETERS = (
    ("bottom_co2", "Bottom-layer CO$_2$ (ppm)"),
    ("XCO2", "XCO$_2$ (ppm)"),
    ("psurf", "$p_s$ (hPa)"),
    ("sulfate_aod760", "Sulfate AOD$_{760}$"),
    ("organic_carbon_aod760", "OC AOD$_{760}$"),
    ("utls_sulfate_aod760", "UTLS AOD$_{760}$"),
    ("sulfate_z0", "Sulfate $z_0$ (km)"),
    ("organic_carbon_z0", "OC $z_0$ (km)"),
    ("utls_sulfate_z0", "UTLS $z_0$ (km)"),
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
    ("mSIF", "mSIF (native)"),
)
COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("retrievals", type=Path, nargs="+", help="Retrieval NetCDFs")
    parser.add_argument("--output", type=Path, help="Output comparison PNG")
    parser.add_argument("--table-output", type=Path, help="Output whitespace table")
    parser.add_argument("--truth-table", type=Path, default=TRUTH_TABLE)
    parser.add_argument(
        "--scene-components", type=Path, default=SCENE_COMPONENTS,
    )
    return parser.parse_args()


def physical_values(names, states, xco2):
    index = {name: i for i, name in enumerate(names)}
    values = {}
    for name, _ in PARAMETERS:
        if name == "XCO2":
            values[name] = np.asarray(xco2, dtype=float)
        elif name == "bottom_co2":
            values[name] = (
                1.0e6 * np.asarray(
                    states[:, index["co2_vmr_layer16"]], dtype=float
                )
            )
        elif name.endswith("_aod760") or name.endswith("_z0"):
            values[name] = np.exp(
                np.asarray(states[:, index["ln_" + name]], dtype=float)
            )
        else:
            data = np.asarray(states[:, index[name]], dtype=float)
            if name == "SIF760":
                data = data * 1.0e7 / 760.0**2
            values[name] = data
    return values


def load_retrieval(path):
    with Dataset(path) as dataset:
        if int(dataset.getncattr("retrieval_complete")) != 1:
            raise RuntimeError(f"retrieval is not marked complete: {path}")
        schema = RetrievalPlotSchema.from_dataset(dataset, path)
        active_states = np.asarray(dataset["state_at_trial"][:], dtype=float)
        active_final = np.asarray(dataset["final_state"][:], dtype=float)
        active_prior = np.asarray(dataset["a_priori_state"][:], dtype=float)
        trials = np.asarray(dataset["trial_index"][:], dtype=int)
        iterations = np.asarray(dataset["iteration_index"][:], dtype=int)
        accepted = np.asarray(dataset["trial_accepted"][:], dtype=int).astype(bool)
        xco2 = np.asarray(dataset["XCO2_at_trial"][:], dtype=float)
        final_xco2 = float(dataset["XCO2"][:])
        prior_xco2 = float(dataset["a_priori_XCO2"][:])
        metadata = {
            "state_index": int(dataset.getncattr("truth_state_index")),
            "perturbation": int(dataset.getncattr("perturbation_index")),
            "retrieval_class": str(dataset.getncattr("measurement_class")),
            "truth_xco2": float(dataset.getncattr("truth_xco2_ppm")),
            "aerosol_case": str(dataset.getncattr("aerosol_case")),
        }

    terminal = np.zeros(len(trials), dtype=bool)
    if np.array_equal(active_final, active_states[-1]):
        terminal[-1] = True
    else:
        active_states = np.vstack([active_states, active_final])
        trials = np.append(trials, trials[-1] + 1)
        iterations = np.append(iterations, iterations[-1] + 1)
        accepted = np.append(accepted, True)
        terminal = np.append(terminal, True)
        xco2 = np.append(xco2, final_xco2)

    names = schema.canonical_names
    states = schema.expand_absolute_state(active_states)
    prior = schema.expand_absolute_state(active_prior)
    values = physical_values(names, states, xco2)
    prior_values = physical_values(names, prior.reshape(1, -1), [prior_xco2])
    return {
        "path": path,
        "names": names,
        "states": states,
        "trials": trials,
        "iterations": iterations,
        "accepted": accepted,
        "terminal": terminal,
        "values": values,
        "prior": {key: value[0] for key, value in prior_values.items()},
        "schema_identity": schema.identity,
        "schema_description": schema.description(),
        "schema": schema,
        **metadata,
    }


def truth_values(record, truth_table, scene_components):
    row = table_row(truth_table, record["state_index"])
    truth = build_truth_values(
        row,
        aerosol_aod760(scene_components, record["aerosol_case"]),
        record["prior"],
        schema=record["schema"],
    )
    return truth, row


def default_stem(records):
    state = records[0]["state_index"]
    perturbations = "_".join(f"{r['perturbation']:02d}" for r in records)
    return f"retrieval_state{state:03d}_perturbations{perturbations}_compact_state_trajectory"


def write_table(path, records):
    columns = [name for name, _ in PARAMETERS]
    with path.open("w", encoding="utf-8") as stream:
        stream.write("# Compact state trajectory; layer-resolved CO2 is replaced by XCO2.\n")
        stream.write("# AOD760 is dimensionless; heights are km; SIF760 is per nm.\n")
        stream.write("# State schema: %s\n" % records[0]["schema_description"])
        stream.write(
            "# perturbation trial iteration accepted terminal " + " ".join(columns) + "\n"
        )
        for record in records:
            for row in range(len(record["trials"])):
                metadata = [
                    record["perturbation"], record["trials"][row],
                    record["iterations"][row], int(record["accepted"][row]),
                    int(record["terminal"][row]),
                ]
                values = [record["values"][name][row] for name in columns]
                stream.write(
                    " ".join(str(value) for value in metadata) + " " +
                    " ".join(f"{value:.12g}" for value in values) + "\n"
                )


def main():
    args = parse_args()
    records = [load_retrieval(path) for path in args.retrievals]
    states = {record["state_index"] for record in records}
    classes = {record["retrieval_class"] for record in records}
    if len(states) != 1 or len(classes) != 1:
        raise RuntimeError("all retrievals must share a truth state and measurement class")
    schemas = {record["schema_identity"] for record in records}
    if len(schemas) != 1:
        raise RuntimeError("all retrievals must share one retrieval-state schema")
    stem = default_stem(records)
    output = args.output or records[0]["path"].with_name(stem + ".png")
    table_output = args.table_output or records[0]["path"].with_name(stem + ".dat")
    truth, truth_row = truth_values(
        records[0], args.truth_table, args.scene_components
    )
    records[0]["schema"].validate_truth_row(
        truth_row, str(args.truth_table)
    )
    write_table(table_output, records)

    fig, axes = plt.subplots(6, 4, figsize=(18, 23))
    fig.subplots_adjust(
        left=0.055, right=0.985, bottom=0.045, top=0.925,
        hspace=0.55, wspace=0.30,
    )
    axes = axes.ravel()
    for panel, (name, label) in enumerate(PARAMETERS):
        axis = axes[panel]
        axis.axhline(truth[name], color="black", linewidth=1.4, label="Truth")
        axis.axhline(
            records[0]["prior"][name], color="0.5", linestyle=":", linewidth=1.4,
            label="Prior",
        )
        for color, record in zip(COLORS, records):
            values = record["values"][name]
            accepted = record["accepted"]
            axis.plot(
                record["trials"][accepted], values[accepted], color=color,
                marker="o", linewidth=1.6,
                label=f"Perturbation {record['perturbation']:02d}",
            )
            rejected = ~accepted
            if np.any(rejected):
                axis.scatter(
                    record["trials"][rejected], values[rejected], color=color,
                    marker="x", s=55, linewidths=1.8, zorder=4,
                )
            terminal = record["terminal"]
            axis.scatter(
                record["trials"][terminal], values[terminal], facecolors="none",
                edgecolors=color, marker="o", s=120, linewidths=2.0, zorder=5,
            )
        logarithmic = name.endswith("aod760") and truth[name] > 0
        if logarithmic:
            axis.set_yscale("log")
        axis.set_title(label, fontsize=10)
        axis.set_xlabel("Evaluated trial / terminal point")
        axis.grid(alpha=0.22)
        if not logarithmic:
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))

    legend_axis = axes[len(PARAMETERS)]
    legend_axis.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    legend_axis.legend(handles, labels, frameon=False, loc="center", fontsize=11)
    for axis in axes[len(PARAMETERS) + 1:]:
        axis.axis("off")
    record = records[0]
    fig.suptitle(
        f"{record['retrieval_class'].capitalize()} retrieval: state "
        f"{record['state_index']:03d}; {co2_truth_label(truth_row)}\n"
        "Compact state progression; x = rejected, "
        "open circle = terminal state\n"
        f"{record['schema_description']}",
        fontsize=16, y=0.975,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)
    print(table_output)


if __name__ == "__main__":
    main()
