#!/usr/bin/env python3
"""Plot truth, prior, and retrieved surface Legendre-moment contributions."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


DEFAULT_TRUTH_TABLE = Path(__file__).resolve().parents[1] / "truth_map" / "true_states.dat"
BANDS = (
    ("o2a", "o2a", "O$_2$ A band", 773.0, 757.0),
    ("weak_co2", "weak", "Weak CO$_2$ band", 1622.0, 1589.0),
    ("strong_co2", "strong", "Strong CO$_2$ band", 2084.0, 2042.0),
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("retrieval", type=Path, help="Completed retrieval NetCDF")
    parser.add_argument("--truth-table", type=Path, default=DEFAULT_TRUTH_TABLE)
    parser.add_argument("--output", type=Path, help="Output PNG")
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
                row = dict(zip(names, values))
                if int(row["index"]) == state_index:
                    return row
    raise RuntimeError(f"truth state {state_index} was not found in {path}")


def main():
    args = parse_args()
    output = args.output or args.retrieval.with_name(
        args.retrieval.stem + "_surface_legendre_moments.png"
    )

    with Dataset(args.retrieval) as dataset:
        if int(dataset.getncattr("retrieval_complete")) != 1:
            raise RuntimeError(f"retrieval is not marked complete: {args.retrieval}")
        names = str(dataset.getncattr("parameter_names")).split()
        state_index = int(dataset.getncattr("truth_state_index"))
        perturbation_index = int(dataset.getncattr("perturbation_index"))
        retrieval_class = str(dataset.getncattr("measurement_class"))
        prior = dict(zip(names, np.asarray(dataset["a_priori_state"][:], float)))
        retrieved = dict(zip(names, np.asarray(dataset["final_state"][:], float)))

    truth = read_truth_row(args.truth_table, state_index)
    fig, axes = plt.subplots(3, 3, figsize=(15, 10), constrained_layout=True)
    moment_labels = ("$c_0P_0(x)$", "$c_1P_1(x)$", "$c_2P_2(x)$")

    for row, (state_band, truth_band, band_label, long_nm, short_nm) in enumerate(BANDS):
        wavenumber = np.linspace(1.0e7 / long_nm, 1.0e7 / short_nm, 1000)
        x = 2.0 * (wavenumber - wavenumber[0]) / (
            wavenumber[-1] - wavenumber[0]
        ) - 1.0
        wavelength = 1.0e7 / wavenumber
        order = np.argsort(wavelength)
        wavelength = wavelength[order]
        basis = (np.ones_like(x), x, 0.5 * (3.0 * x**2 - 1.0))

        for moment, axis in enumerate(axes[row]):
            coefficient = f"P{moment}"
            state_name = f"{state_band}_surface_{coefficient}"
            truth_name = f"{truth_band}_{coefficient}"
            series = (
                ("Prior", prior[state_name] * basis[moment], "0.55", ":"),
                ("Truth", float(truth[truth_name]) * basis[moment], "black", "-"),
                ("Retrieved", retrieved[state_name] * basis[moment], "#d62728", "--"),
            )
            for label, values, color, linestyle in series:
                axis.plot(
                    wavelength, values[order], color=color, linestyle=linestyle,
                    linewidth=1.7, label=label,
                )
            if row == 0:
                axis.set_title(moment_labels[moment], fontsize=13)
            if moment == 0:
                axis.set_ylabel(f"{band_label}\nAlbedo contribution")
            if row == 2:
                axis.set_xlabel("Wavelength (nm)")
            axis.grid(alpha=0.22)
            axis.margins(x=0)
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))

    axes[0, 0].legend(frameon=False, loc="best")
    fig.suptitle(
        f"{retrieval_class.capitalize()} retrieval: state {state_index:03d}, "
        f"perturbation {perturbation_index:02d}\n"
        "Surface Legendre contributions; $x$ is normalized wavenumber",
        fontsize=15,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
