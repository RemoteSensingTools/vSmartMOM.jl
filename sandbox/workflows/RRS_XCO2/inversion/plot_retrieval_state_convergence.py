#!/usr/bin/env python3
"""Plot a compact retrieval-state path together with OE convergence diagnostics.

Layer-resolved CO2 is represented by the saved dry-air-column XCO2 diagnostic.
Logarithmic aerosol coordinates are transformed back to AOD and height.  The
terminal state, which is evaluated once after the final accepted LM proposal,
is shown explicitly after the trial-state history.
"""

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


HERE = Path(__file__).resolve().parent
WORKFLOW_ROOT = HERE.parent
DATA_ROOT = Path(os.environ.get("RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
TRUTH_TABLE = DATA_ROOT / "truth_map" / "true_states.dat"
SCENE_COMPONENTS = DATA_ROOT / "truth_map" / "scene_components.dat"

PARAMETERS = (
    ("XCO2", "XCO$_2$ (ppm)"),
    ("psurf", "$p_s$ (hPa)"),
    ("sulfate_aod760", "Sulphate AOD$_{760}$"),
    ("organic_carbon_aod760", "OC AOD$_{760}$"),
    ("utls_sulfate_aod760", "UTLS AOD$_{760}$"),
    ("sulfate_z0", "Sulphate $z_0$ (km)"),
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
BAND_LABELS = ("O$_2$ A", "Weak CO$_2$", "Strong CO$_2$")
BAND_COLORS = ("#1f77b4", "#2ca02c", "#d62728")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("retrieval", type=Path, help="Completed retrieval NetCDF")
    parser.add_argument("--output", type=Path, help="Output PNG")
    parser.add_argument(
        "--table-output", type=Path,
        help="Whitespace table of the compact physical state trajectory",
    )
    return parser.parse_args()


def table_row(path, state_index):
    names = None
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("# index "):
                names = stripped[2:].split()
            elif stripped and not stripped.startswith("#"):
                values = stripped.split()
                if names is not None and int(values[0]) == state_index:
                    return dict(zip(names, values))
    raise RuntimeError(f"truth state {state_index} was not found in {path}")


def aerosol_aod760(path, aerosol_case):
    active = False
    names = None
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped == "[AEROSOL_CASES]":
                active = True
                continue
            if active and stripped.startswith("["):
                break
            if not active:
                continue
            if stripped.startswith("# case"):
                names = stripped[2:].split()
            elif stripped and not stripped.startswith("#") and names is not None:
                row = dict(zip(names, stripped.split()))
                if row["case"] == aerosol_case:
                    return {
                        "sulfate_aod760": float(row["sulfate_AOD760"]),
                        "organic_carbon_aod760": float(row["organic_AOD760"]),
                        "utls_sulfate_aod760": float(row["utls_sulfate_AOD760"]),
                    }
    raise RuntimeError(f"aerosol case {aerosol_case!r} was not found in {path}")


def physical_state(names, state, xco2):
    native = dict(zip(names, state))
    values = {"XCO2": float(xco2), "psurf": native["psurf"]}
    for species in ("sulfate", "organic_carbon", "utls_sulfate"):
        values[f"{species}_aod760"] = np.exp(native[f"ln_{species}_aod760"])
        values[f"{species}_z0"] = np.exp(native[f"ln_{species}_z0"])
    for band in ("o2a", "weak_co2", "strong_co2"):
        for order in range(3):
            key = f"{band}_surface_P{order}"
            values[key] = native[key]
    values["SIF760"] = native["SIF760"] * 1.0e7 / 760.0**2
    values["mSIF"] = native["mSIF"]
    return values


def truth_values(row, aerosol_truth, prior):
    truth = dict(prior)
    truth.update(aerosol_truth)
    truth.update({
        "XCO2": float(row["xco2_ppm"]),
        "psurf": float(row["psurf_hpa"]),
        "SIF760": float(row["SIF760"]) * 1.0e7 / 760.0**2,
        "mSIF": float(row["mSIF"]),
    })
    for state_band, truth_band in (
        ("o2a", "o2a"), ("weak_co2", "weak"), ("strong_co2", "strong")
    ):
        for order in range(3):
            truth[f"{state_band}_surface_P{order}"] = float(
                row[f"{truth_band}_P{order}"]
            )
    return truth


def write_state_table(path, truth, prior_state, trial_states, final_state,
                      prior_native, first_trial_native, trials):
    """Write parameters as rows and distinct retrieval states as columns."""
    merge_prior = np.array_equal(prior_native, first_trial_native)
    if merge_prior:
        state_labels = [f"apriori_trial{trials[0]:02d}"]
        state_values = [trial_states[0]]
        for trial, values in zip(trials[1:], trial_states[1:]):
            state_labels.append(f"trial{trial:02d}")
            state_values.append(values)
    else:
        state_labels = ["apriori"] + [f"trial{trial:02d}" for trial in trials]
        state_values = [prior_state] + trial_states
    state_labels.append("terminal")
    state_values.append(final_state)
    units = {
        "XCO2": "ppm",
        "psurf": "hPa",
        "sulfate_aod760": "1",
        "organic_carbon_aod760": "1",
        "utls_sulfate_aod760": "1",
        "sulfate_z0": "km",
        "organic_carbon_z0": "km",
        "utls_sulfate_z0": "km",
        "SIF760": "mW_m-2_sr-1_nm-1",
        "mSIF": "native",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# Compact physical state trajectory; layer CO2 is replaced by XCO2.\n"
        )
        stream.write(
            "# The a priori and trial 01 are merged only when their native state "
            "vectors are exactly identical.\n"
        )
        stream.write(
            "# parameter units truth " + " ".join(state_labels) + "\n"
        )
        for key, _ in PARAMETERS:
            unit = units.get(key, "1")
            values = " ".join(f"{state[key]:.12g}" for state in state_values)
            stream.write(f"{key} {unit} {truth[key]:.12g} {values}\n")


def main():
    args = parse_args()
    output = args.output or args.retrieval.with_name(
        args.retrieval.stem + "_state_and_convergence.png"
    )
    table_output = args.table_output or args.retrieval.with_name(
        args.retrieval.stem + "_compact_state_trajectory.dat"
    )
    with Dataset(args.retrieval) as dataset:
        if int(dataset.getncattr("retrieval_complete")) != 1:
            raise RuntimeError(f"retrieval is not complete: {args.retrieval}")
        names = str(dataset.getncattr("parameter_names")).split()
        states = np.asarray(dataset["state_at_trial"][:], dtype=float)
        final_state = np.asarray(dataset["final_state"][:], dtype=float)
        prior_state = np.asarray(dataset["a_priori_state"][:], dtype=float)
        xco2 = np.asarray(dataset["XCO2_at_trial"][:], dtype=float)
        final_xco2 = float(dataset["XCO2"][:])
        prior_xco2 = float(dataset["a_priori_XCO2"][:])
        trials = np.asarray(dataset["trial_index"][:], dtype=int)
        accepted = np.asarray(dataset["trial_accepted"][:], dtype=int).astype(bool)
        divergent = np.asarray(dataset["trial_divergent"][:], dtype=int).astype(bool)
        gamma = np.asarray(dataset["gamma"][:], dtype=float)
        ratio = np.asarray(dataset["reduction_ratio"][:], dtype=float)
        d_sigma = np.asarray(dataset["d_sigma_sq_scaled"][:], dtype=float)
        chi = np.asarray(dataset["band_reduced_chi_squared"][:], dtype=float)
        final_chi = np.asarray(
            dataset["final_band_reduced_chi_squared"][:], dtype=float
        )
        cost = np.asarray(dataset["total_cost"][:], dtype=float)
        seconds = np.asarray(dataset["evaluation_seconds"][:], dtype=float)
        state_index = int(dataset.getncattr("truth_state_index"))
        perturbation = int(dataset.getncattr("perturbation_index"))
        retrieval_class = str(dataset.getncattr("measurement_class"))
        aerosol_case = str(dataset.getncattr("aerosol_case"))
        convergence_threshold = float(dataset.getncattr("convergence_threshold"))
        fit_threshold = float(dataset.getncattr("maximum_band_chi_squared"))
        converged = bool(dataset.getncattr("converged"))

    physical_trials = [
        physical_state(names, state, column_xco2)
        for state, column_xco2 in zip(states, xco2)
    ]
    physical_final = physical_state(names, final_state, final_xco2)
    physical_prior = physical_state(names, prior_state, prior_xco2)
    points = np.append(trials, trials[-1] + 1)
    row = table_row(TRUTH_TABLE, state_index)
    truth = truth_values(
        row, aerosol_aod760(SCENE_COMPONENTS, aerosol_case), physical_prior
    )
    write_state_table(
        table_output, truth, physical_prior, physical_trials, physical_final,
        prior_state, states[0], trials,
    )

    fig, axes = plt.subplots(5, 5, figsize=(21, 19))
    axes = axes.ravel()
    for axis, (key, label) in zip(axes, PARAMETERS):
        values = np.asarray(
            [record[key] for record in physical_trials] + [physical_final[key]]
        )
        axis.axhline(truth[key], color="black", linewidth=1.35, label="Truth")
        axis.axhline(
            physical_prior[key], color="0.5", linestyle=":", linewidth=1.35,
            label="A priori",
        )
        axis.plot(points, values, color="#1f77b4", marker="o", linewidth=1.6)
        rejected = ~accepted
        if np.any(rejected):
            axis.scatter(
                trials[rejected], values[:-1][rejected], color="#d62728",
                marker="x", s=60, linewidths=1.8, zorder=4,
            )
        axis.scatter(
            points[-1], values[-1], facecolors="none", edgecolors="#1f77b4",
            marker="o", s=110, linewidths=2.0, zorder=5,
        )
        axis.set_title(label, fontsize=10)
        axis.set_xticks(points)
        axis.grid(alpha=0.22)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))

    axis = axes[19]
    for column, (label, color) in enumerate(zip(BAND_LABELS, BAND_COLORS)):
        axis.plot(
            points, np.append(chi[:, column], final_chi[column]), marker="o",
            color=color, linewidth=1.5, label=label,
        )
    axis.axhline(fit_threshold, color="black", linestyle="--", linewidth=1.1)
    axis.set_yscale("log")
    axis.set_title("Reduced $\\chi^2$ by band")
    axis.legend(frameon=False, fontsize=8)

    axis = axes[20]
    axis.plot(trials, d_sigma, color="#9467bd", marker="o", linewidth=1.6)
    axis.axhline(
        convergence_threshold, color="black", linestyle="--", linewidth=1.1,
        label=f"threshold={convergence_threshold:g}",
    )
    axis.set_yscale("log")
    axis.set_title("Posterior-normalized state step")
    axis.set_ylabel("$d_\\sigma^2/N_{eff}$")
    axis.legend(frameon=False, fontsize=8)

    axis = axes[21]
    axis.plot(trials, gamma, color="#ff7f0e", marker="o", label="$\\gamma$")
    finite = np.isfinite(ratio)
    twin = axis.twinx()
    twin.plot(
        trials[finite], ratio[finite], color="#2ca02c", marker="s",
        label="reduction ratio",
    )
    axis.set_title("LM damping and reduction ratio")
    axis.set_ylabel("$\\gamma$")
    twin.set_ylabel("actual / forecast reduction")

    axis = axes[22]
    axis.plot(trials, cost, color="#8c564b", marker="o", linewidth=1.6)
    axis.set_yscale("log")
    axis.set_title("OE total cost")

    axis = axes[23]
    axis.bar(trials, seconds, color="#17becf", alpha=0.8)
    axis.set_title("Forward + Jacobian evaluation time")
    axis.set_ylabel("seconds")

    for axis in axes[19:24]:
        axis.set_xticks(trials if axis is not axes[19] else points)
        axis.set_xlabel("Trial" if axis is not axes[19] else "Trial / terminal")
        axis.grid(alpha=0.22)

    axes[24].axis("off")
    axes[24].plot([], [], color="black", linewidth=1.35, label="Truth")
    axes[24].plot([], [], color="0.5", linestyle=":", label="A priori")
    axes[24].plot([], [], color="#1f77b4", marker="o", label="Evaluated state")
    axes[24].scatter(
        [], [], facecolors="none", edgecolors="#1f77b4", marker="o", s=90,
        label="Terminal state",
    )
    axes[24].legend(frameon=False, loc="upper center", fontsize=10)
    axes[24].text(
        0.5, 0.32,
        f"converged = {converged}\n"
        f"final $d_\\sigma^2/N_{{eff}}$ = {d_sigma[-1]:.3g}\n"
        f"final band $\\chi_r^2$ = "
        + ", ".join(f"{value:.3f}" for value in final_chi),
        ha="center", va="center", transform=axes[24].transAxes, fontsize=10,
    )

    fig.suptitle(
        f"{retrieval_class.capitalize()} retrieval: state {state_index:03d}, "
        f"perturbation {perturbation:02d}\n"
        "Compact state evolution and convergence diagnostics",
        fontsize=16, y=0.992,
    )
    fig.subplots_adjust(
        left=0.05, right=0.975, bottom=0.04, top=0.95, hspace=0.55, wspace=0.34
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)
    print(table_output)


if __name__ == "__main__":
    main()
