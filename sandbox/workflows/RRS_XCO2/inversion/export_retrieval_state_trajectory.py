#!/usr/bin/env python3
"""Export every evaluated retrieval state in physical, human-readable units."""

import argparse
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
DEFAULT_TRUTH_TABLE = DATA_ROOT / "truth_map" / "true_states.dat"
BAND_SURFACE_NAMES = {
    "o2a": "o2a",
    "weak_co2": "weak",
    "strong_co2": "strong",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("retrieval", type=Path, help="Completed retrieval NetCDF")
    parser.add_argument("--output", type=Path, help="Output whitespace table")
    parser.add_argument(
        "--truth-table", type=Path, default=DEFAULT_TRUTH_TABLE,
        help="Truth-map state table",
    )
    parser.add_argument(
        "--surface-output", type=Path,
        help="Output truth/prior/retrieved surface comparison",
    )
    return parser.parse_args()


def read_truth_row(path, state_index):
    names = None
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("# index "):
                names = stripped[2:].split()
            elif stripped and not stripped.startswith("#"):
                if names is None:
                    raise RuntimeError(f"truth header was not found in {path}")
                values = stripped.split()
                if len(values) != len(names):
                    raise RuntimeError(f"malformed truth row in {path}")
                row = dict(zip(names, values))
                if int(row["index"]) == state_index:
                    return row
    raise RuntimeError(f"truth state {state_index} was not found in {path}")


def physical_state(names, state):
    output_names = []
    output_values = []
    for name, value in zip(names, state):
        if name.startswith("co2_vmr_"):
            output_names.append(name.replace("co2_vmr_", "co2_ppm_"))
            output_values.append(value * 1.0e6)
        elif name.startswith("ln_") and name.endswith("_aod760"):
            output_names.append(name[3:])
            output_values.append(np.exp(value))
        elif name.startswith("ln_") and name.endswith("_z0"):
            output_names.append(name[3:] + "_km")
            output_values.append(np.exp(value))
        elif name == "SIF760":
            output_names.extend(["SIF760_per_cm1", "SIF760_per_nm"])
            output_values.extend([value, value * 1.0e7 / 760.0**2])
        elif name == "mSIF":
            output_names.append("mSIF_per_cm1_squared")
            output_values.append(value)
        else:
            output_names.append(name)
            output_values.append(value)
    return output_names, np.asarray(output_values, dtype=float)


def main():
    args = parse_args()
    output = args.output or args.retrieval.with_name(
        args.retrieval.stem + "_state_trajectory.dat"
    )
    surface_output = args.surface_output or args.retrieval.with_name(
        args.retrieval.stem + "_surface_truth_vs_retrieved.dat"
    )

    with Dataset(args.retrieval) as dataset:
        if int(dataset.getncattr("retrieval_complete")) != 1:
            raise RuntimeError(f"retrieval is not marked complete: {args.retrieval}")
        names = str(dataset.getncattr("parameter_names")).split()
        state_index = int(dataset.getncattr("truth_state_index"))
        prior = np.asarray(dataset["a_priori_state"][:], dtype=float)
        final = np.asarray(dataset["final_state"][:], dtype=float)
        states = np.asarray(dataset["state_at_trial"][:], dtype=float)
        trial = np.asarray(dataset["trial_index"][:], dtype=int)
        iteration = np.asarray(dataset["iteration_index"][:], dtype=int)
        accepted = np.asarray(dataset["trial_accepted"][:], dtype=int)
        divergent = np.asarray(dataset["trial_divergent"][:], dtype=int)
        gamma = np.asarray(dataset["gamma"][:], dtype=float)
        d_sigma = np.asarray(dataset["d_sigma_sq_scaled"][:], dtype=float)
        chi_squared = np.asarray(dataset["band_reduced_chi_squared"][:], dtype=float)
        final_chi_squared = np.asarray(
            dataset["final_band_reduced_chi_squared"][:], dtype=float
        )

    if states.shape != (len(trial), len(names)):
        raise RuntimeError("state trajectory dimensions disagree with metadata")
    physical_names, physical_prior = physical_state(names, prior)
    physical_states = np.vstack(
        [physical_state(names, state)[1] for state in states]
    )

    metadata_names = [
        "trial", "iteration", "terminal", "accepted", "divergent", "gamma",
        "d_sigma_sq_scaled", "chi2_o2a", "chi2_weak_co2", "chi2_strong_co2",
    ]
    prior_metadata = np.array(
        [0, 0, 0, 1, 0, np.nan, np.nan, np.nan, np.nan, np.nan], dtype=float
    )
    terminal = np.zeros(len(trial), dtype=int)
    final_is_last_trial = np.array_equal(final, states[-1])
    if final_is_last_trial:
        terminal[-1] = 1
    trial_metadata = np.column_stack(
        [trial, iteration, terminal, accepted, divergent, gamma, d_sigma, chi_squared]
    )
    rows = [
        np.concatenate([prior_metadata, physical_prior]),
        *np.column_stack([trial_metadata, physical_states]),
    ]
    if not final_is_last_trial:
        _, physical_final = physical_state(names, final)
        final_metadata = np.array(
            [trial[-1] + 1, iteration[-1] + 1, 1, 1, 0, np.nan, np.nan,
             *final_chi_squared],
            dtype=float,
        )
        rows.append(np.concatenate([final_metadata, physical_final]))
    table = np.vstack(rows)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        stream.write("# Retrieval state trajectory in physical coordinates.\n")
        stream.write("# Row trial=0 is the a priori; subsequent rows are every evaluated LM trial.\n")
        stream.write(
            "# terminal=1 identifies the state used for the saved final products. "
            "A converged retrieval may require one extra terminal evaluation after "
            "the last LM trial.\n"
        )
        stream.write("# CO2 is ppm; AOD760 is dimensionless; aerosol z0 is km.\n")
        stream.write("# SIF760 is reported in both native per-cm^-1 and per-nm units.\n")
        stream.write("# " + " ".join(metadata_names + physical_names) + "\n")
        np.savetxt(stream, table, fmt="%.12g")

    truth = read_truth_row(args.truth_table, state_index)
    state_lookup = dict(zip(names, final))
    prior_lookup = dict(zip(names, prior))
    surface_output.parent.mkdir(parents=True, exist_ok=True)
    with surface_output.open("w", encoding="utf-8") as stream:
        stream.write("# Truth, prior, and terminal surface Legendre coefficients.\n")
        stream.write(f"# truth_state_index={state_index} source={args.truth_table}\n")
        stream.write(
            "# band coefficient truth prior retrieved prior_minus_truth "
            "retrieved_minus_truth retrieved_percent_error\n"
        )
        for state_band, truth_band in BAND_SURFACE_NAMES.items():
            for coefficient in ("P0", "P1", "P2"):
                state_name = f"{state_band}_surface_{coefficient}"
                truth_name = f"{truth_band}_{coefficient}"
                true_value = float(truth[truth_name])
                prior_value = prior_lookup[state_name]
                retrieved_value = state_lookup[state_name]
                percent_error = 100.0 * (retrieved_value - true_value) / true_value
                stream.write(
                    f"{state_band} {coefficient} {true_value:.12g} "
                    f"{prior_value:.12g} {retrieved_value:.12g} "
                    f"{prior_value - true_value:.12g} "
                    f"{retrieved_value - true_value:.12g} {percent_error:.8g}\n"
                )
    print(output)
    print(surface_output)


if __name__ == "__main__":
    main()
