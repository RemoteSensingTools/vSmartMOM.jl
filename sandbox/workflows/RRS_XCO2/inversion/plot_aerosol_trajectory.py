#!/usr/bin/env python3
"""Plot aerosol-height and AOD evolution across evaluated OE trial states."""

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
DEFAULT_COMPONENTS = DATA_ROOT / "truth_map" / "scene_components.dat"
SPECIES = (
    ("sulfate", "Sulfate"),
    ("organic_carbon", "Organic carbon"),
    ("utls_sulfate", "UTLS sulfate"),
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("retrieval", type=Path, help="Completed retrieval NetCDF")
    parser.add_argument(
        "--scene-components", type=Path, default=DEFAULT_COMPONENTS,
        help="Truth-map scene_components.dat file",
    )
    parser.add_argument("--output", type=Path, help="Output PNG")
    return parser.parse_args()


def read_truth_aod760(path, aerosol_case):
    in_cases = False
    names = None
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped == "[AEROSOL_CASES]":
                in_cases = True
                continue
            if in_cases and stripped.startswith("["):
                break
            if not in_cases:
                continue
            if stripped.startswith("# case"):
                names = stripped[2:].split()
            elif stripped and not stripped.startswith("#"):
                if names is None:
                    raise RuntimeError(f"AEROSOL_CASES header not found in {path}")
                row = dict(zip(names, stripped.split()))
                if row["case"] == aerosol_case:
                    return {
                        "sulfate": float(row["sulfate_AOD760"]),
                        "organic_carbon": float(row["organic_AOD760"]),
                        "utls_sulfate": float(row["utls_sulfate_AOD760"]),
                    }
    raise RuntimeError(f"aerosol case {aerosol_case!r} not found in {path}")


def main():
    args = parse_args()
    output = args.output or args.retrieval.with_name(
        args.retrieval.stem + "_aerosol_trajectory.png"
    )

    with Dataset(args.retrieval) as dataset:
        if int(dataset.getncattr("retrieval_complete")) != 1:
            raise RuntimeError(f"retrieval is not marked complete: {args.retrieval}")
        names = str(dataset.getncattr("parameter_names")).split()
        index = {name: i for i, name in enumerate(names)}
        states = np.asarray(dataset["state_at_trial"][:], float)
        prior = np.asarray(dataset["a_priori_state"][:], float)
        accepted = np.asarray(dataset["trial_accepted"][:], int).astype(bool)
        trials = np.asarray(dataset["trial_index"][:], int)
        state_index = int(dataset.getncattr("truth_state_index"))
        perturbation_index = int(dataset.getncattr("perturbation_index"))
        retrieval_class = str(dataset.getncattr("measurement_class"))
        aerosol_case = str(dataset.getncattr("aerosol_case"))

    truth_aod = read_truth_aod760(args.scene_components, aerosol_case)
    accepted_trials = trials[accepted]
    rejected_trials = trials[~accepted]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), constrained_layout=True)

    for column, (key, label) in enumerate(SPECIES):
        height_name = f"ln_{key}_z0"
        aod_name = f"ln_{key}_aod760"
        height = np.exp(states[:, index[height_name]])
        aod = np.exp(states[:, index[aod_name]])
        prior_height = np.exp(prior[index[height_name]])
        prior_aod = np.exp(prior[index[aod_name]])

        for axis, values, truth, prior_value, ylabel in (
            (axes[0, column], height, prior_height, prior_height, "Height $z_0$ (km)"),
            (axes[1, column], aod, truth_aod[key], prior_aod, "AOD$_{760}$"),
        ):
            axis.axhline(
                truth, color="black", linewidth=1.6, label="Truth" if column == 0 else None,
            )
            axis.axhline(
                prior_value, color="0.5", linestyle=":", linewidth=1.6,
                label="Prior" if column == 0 else None,
            )
            axis.plot(
                accepted_trials, values[accepted], color="#1f77b4", marker="o",
                linewidth=1.8, label="Accepted state" if column == 0 else None,
            )
            if rejected_trials.size:
                axis.scatter(
                    rejected_trials, values[~accepted], color="#d62728", marker="x",
                    s=65, linewidths=2.0,
                    label="Rejected trial" if column == 0 else None,
                    zorder=4,
                )
            axis.set_xticks(trials)
            axis.grid(alpha=0.22)
            axis.set_ylabel(ylabel)
            if ylabel == "AOD$_{760}$" and truth > 0:
                axis.set_yscale("log")
        axes[0, column].set_title(label)
        axes[1, column].set_xlabel("Evaluated LM trial")

    axes[0, 0].legend(frameon=False, loc="best")
    fig.suptitle(
        f"{retrieval_class.capitalize()} retrieval: state {state_index:03d}, "
        f"perturbation {perturbation_index:02d}\n"
        "Aerosol height and optical-depth trajectory",
        fontsize=15,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
