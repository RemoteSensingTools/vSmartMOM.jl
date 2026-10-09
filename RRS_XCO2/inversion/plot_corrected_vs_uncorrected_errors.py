#!/usr/bin/env python3
"""Compare corrected and uncorrected retrieval displacements for one state.

Every marker represents one matched noise perturbation.  By default, filled
markers show ``x_retrieved - x_truth`` and open markers show the terminal
linear prediction ``x_retrieved,0 - x_truth + G @ delta_y`` for the same
realization, where index 0 denotes the noiseless retrieval (perturbation 11).
The horizontal coordinate is corrected and the vertical coordinate is
uncorrected.  The gain overlay can be disabled, and prior-referenced retrieval
displacements remain available explicitly from the command line.

Layer-resolved CO2 is condensed to both the bottom-layer VMR and XCO2.  The
former is the primary coordinate for bottom-layer campaigns; for uniform
full-column campaigns it equals the truth VMR.  Logarithmic aerosol retrieval
coordinates are transformed to physical AOD and height.  Gain increments are
mapped with the tangent of that physical transformation at the terminal state;
they are not exponentiated as though ``G @ delta_y`` were a complete state.

The 3-by-7 panel placement is intentionally fixed to support direct comparison
between truth states.  The two unused cells in columns 1 and 2 are left blank.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from netCDF4 import Dataset

from retrieval_plot_schema import RetrievalPlotSchema

from plot_retrieval_state_convergence import (
    SCENE_COMPONENTS,
    TRUTH_TABLE,
    aerosol_aod760,
    co2_truth_label,
    physical_state,
    table_row,
    truth_values,
)


HERE = Path(__file__).resolve().parent

# Normalized dry-air columns for the shared 16-layer, 1000-hPa retrieval
# profile.  These are generated and documented in retrieval_setup/build_apriori.jl.
# Only the bottom-layer column changes when p_s is retrieved.
DRY_COLUMN_FRACTIONS_1000 = np.asarray([
    0.06271458366696699, 0.06267999123579230,
    0.06264827279449964, 0.06261831034665068,
    0.06258999065374835, 0.06255915670625688,
    0.06253127935109044, 0.06250621030757231,
    0.06248208275638748, 0.06245969764265345,
    0.06243817296080489, 0.06241161154141432,
    0.06238468284683572, 0.06235227383324737,
    0.06232590314728601, 0.06229778020879323,
])
REFERENCE_SURFACE_PRESSURE_HPA = 1000.0
BOTTOM_LAYER_TOP_PRESSURE_HPA = 937.50875
SIF_WAVENUMBER_TO_WAVELENGTH_760 = 1.0e7 / 760.0**2

# key, panel title, row, column, physical plotting unit
PANELS = (
    ("bottom_co2", "Bottom-layer CO$_2$", 0, 0, "ppm"),
    ("XCO2", "XCO$_2$", 1, 0, "ppm"),
    ("psurf", "$p_s$", 2, 0, "hPa"),
    ("SIF760", "SIF$_{760}$", 0, 1,
     "mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$"),
    ("mSIF", "$m_{SIF}$", 1, 1, "native"),
    ("sulfate_aod760", "Sulphate AOD$_{760}$", 0, 2, "1"),
    ("organic_carbon_aod760", "Organic AOD$_{760}$", 1, 2, "1"),
    ("utls_sulfate_aod760", "UTLS AOD$_{760}$", 2, 2, "1"),
    ("sulfate_z0", "Sulphate $z_0$", 0, 3, "km"),
    ("organic_carbon_z0", "Organic $z_0$", 1, 3, "km"),
    ("utls_sulfate_z0", "UTLS $z_0$", 2, 3, "km"),
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
        "--inversion-root", type=Path, default=HERE,
        help="Directory containing corrected/ and uncorrected/",
    )
    parser.add_argument(
        "--truth-table", type=Path, default=TRUTH_TABLE,
        help="Campaign true_states.dat (default: full-column truth map)",
    )
    parser.add_argument(
        "--scene-components", type=Path, default=SCENE_COMPONENTS,
        help="Campaign scene_components.dat (default: full-column truth map)",
    )
    parser.add_argument("--output", type=Path, help="Output PNG")
    parser.add_argument(
        "--table-output", type=Path,
        help="Whitespace table containing every plotted value",
    )
    parser.add_argument(
        "--reference", choices=("truth", "prior"), default="truth",
        help=(
            "Reference for filled markers (default: truth). Use prior for "
            "retrieval displacement from the a priori."
        ),
    )
    gain_group = parser.add_mutually_exclusive_group()
    gain_group.add_argument(
        "--gain-prediction", "--gain-noise", dest="gain_noise",
        action="store_true",
        help=(
            "Overlay open markers for x_retr,0 - x_truth + terminal "
            "G @ injected_noise (default)"
        ),
    )
    gain_group.add_argument(
        "--no-gain-prediction", "--no-gain-noise", dest="gain_noise",
        action="store_false",
        help="Plot only the filled terminal retrieval displacements",
    )
    parser.set_defaults(gain_noise=True)
    return parser.parse_args()


def injected_noise(dataset):
    """Read the exact additive noise, with fallback for archived outputs."""
    if "injected_measurement_noise" in dataset.variables:
        return np.asarray(
            dataset["injected_measurement_noise"][:], dtype=float
        )
    if (
        "measurement_perturbed" in dataset.variables
        and "measurement_noiseless" in dataset.variables
    ):
        return (
            np.asarray(dataset["measurement_perturbed"][:], dtype=float)
            - np.asarray(dataset["measurement_noiseless"][:], dtype=float)
        )
    if (
        "normalized_noise_draw" in dataset.variables
        and "noise_standard_deviation" in dataset.variables
    ):
        return (
            np.asarray(dataset["normalized_noise_draw"][:], dtype=float)
            * np.asarray(dataset["noise_standard_deviation"][:], dtype=float)
        )
    raise RuntimeError("retrieval does not store a recoverable noise realization")


def terminal_native_noise_response(dataset, nstate):
    """Return G @ delta_y in the saved active retrieval coordinates."""
    noise = injected_noise(dataset)
    gain = np.asarray(dataset["gain_matrix"][:], dtype=float)
    if gain.shape == (noise.size, nstate):
        # Julia/NCDatasets writes dimensions in column-major order; netCDF4
        # exposes this variable as (measurement, state).
        gain = gain.T
    elif gain.shape != (nstate, noise.size):
        raise RuntimeError(
            f"gain shape {gain.shape} is incompatible with state={nstate} "
            f"and measurement={noise.size}"
        )
    response = gain @ noise
    if response.shape != (nstate,) or not np.all(np.isfinite(response)):
        raise RuntimeError("G @ injected_noise produced an invalid state increment")
    return response


def xco2_value_and_tangent(names, native_state, native_delta, fixed_upper_ppm):
    """Map a native state increment to its first-order XCO2 increment."""
    native = dict(zip(names, native_state))
    delta = dict(zip(names, native_delta))
    psurf = native["psurf"]
    bottom_thickness = psurf - BOTTOM_LAYER_TOP_PRESSURE_HPA
    reference_thickness = (
        REFERENCE_SURFACE_PRESSURE_HPA - BOTTOM_LAYER_TOP_PRESSURE_HPA
    )
    if bottom_thickness <= 0.0:
        raise RuntimeError(
            f"terminal surface pressure {psurf} hPa is above the bottom-layer top"
        )

    weights = DRY_COLUMN_FRACTIONS_1000.copy()
    weights[-1] *= bottom_thickness / reference_thickness
    total = float(np.sum(weights))
    co2 = np.empty(16, dtype=float)
    co2[:4] = fixed_upper_ppm * 1.0e-6
    for layer in range(5, 17):
        co2[layer - 1] = native[f"co2_vmr_layer{layer:02d}"]
    mean_vmr = float(np.dot(co2, weights) / total)

    xco2_delta = 0.0
    for layer in range(5, 17):
        gradient = 1.0e6 * weights[layer - 1] / total
        xco2_delta += gradient * delta[f"co2_vmr_layer{layer:02d}"]

    bottom_weight_derivative = (
        DRY_COLUMN_FRACTIONS_1000[-1] / reference_thickness
    )
    psurf_gradient = (
        1.0e6 * bottom_weight_derivative * (co2[-1] - mean_vmr) / total
    )
    xco2_delta += psurf_gradient * delta["psurf"]
    return 1.0e6 * mean_vmr, xco2_delta


def physical_noise_response(names, native_state, native_delta, fixed_upper_ppm):
    """Transform a terminal native-coordinate increment to panel units."""
    native = dict(zip(names, native_state))
    delta = dict(zip(names, native_delta))
    _, xco2_delta = xco2_value_and_tangent(
        names, native_state, native_delta, fixed_upper_ppm
    )
    values = {
        "bottom_co2": 1.0e6 * delta["co2_vmr_layer16"],
        "XCO2": xco2_delta,
        "psurf": delta["psurf"],
    }
    for species in ("sulfate", "organic_carbon", "utls_sulfate"):
        aod_key = f"ln_{species}_aod760"
        height_key = f"ln_{species}_z0"
        values[f"{species}_aod760"] = np.exp(native[aod_key]) * delta[aod_key]
        values[f"{species}_z0"] = np.exp(native[height_key]) * delta[height_key]
    for band in ("o2a", "weak_co2", "strong_co2"):
        for order in range(3):
            key = f"{band}_surface_P{order}"
            values[key] = delta[key]
    values["SIF760"] = delta["SIF760"] * SIF_WAVENUMBER_TO_WAVELENGTH_760
    values["mSIF"] = delta["mSIF"]
    return values


def completed_retrievals(directory, state_index, expected_class, load_gain_noise):
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
                schema = RetrievalPlotSchema.from_dataset(dataset, str(path))
                active_final_state = np.asarray(
                    dataset["final_state"][:], dtype=float
                )
                active_prior_state = np.asarray(
                    dataset["a_priori_state"][:], dtype=float
                )
                names = schema.canonical_names
                final_state = schema.expand_absolute_state(active_final_state)
                prior_state = schema.expand_absolute_state(active_prior_state)
                final_xco2 = float(dataset["XCO2"][:])
                prior_xco2 = float(dataset["a_priori_XCO2"][:])
                gain_noise = None
                if load_gain_noise:
                    active_delta = terminal_native_noise_response(
                        dataset, len(schema.active_names)
                    )
                    native_delta = schema.expand_tangent(active_delta)
                    fixed_upper_ppm = float(
                        dataset.getncattr("fixed_upper_co2_ppm")
                    )
                    estimated_xco2, _ = xco2_value_and_tangent(
                        names, final_state, np.zeros_like(final_state),
                        fixed_upper_ppm,
                    )
                    if not np.isclose(
                        estimated_xco2, final_xco2, rtol=3.0e-7, atol=2.0e-5
                    ):
                        raise RuntimeError(
                            "the plotting XCO2 mapping does not reproduce the "
                            f"saved terminal diagnostic ({estimated_xco2} versus "
                            f"{final_xco2} ppm)"
                        )
                    gain_noise = physical_noise_response(
                        names, final_state, native_delta, fixed_upper_ppm
                    )
                metadata = {
                    "surface": str(dataset.getncattr("surface")),
                    "aerosol_case": str(dataset.getncattr("aerosol_case")),
                    "truth_xco2": float(dataset.getncattr("truth_xco2_ppm")),
                    "truth_bottom_co2": (
                        float(dataset.getncattr("truth_bottom_co2_ppm"))
                        if "truth_bottom_co2_ppm" in dataset.ncattrs() else None
                    ),
                    "plot_schema_identity": schema.identity,
                    "plot_schema_description": schema.description(),
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
            "gain_noise": gain_noise,
            "metadata": metadata,
        }
    return records


def matched_records(inversion_root, state_index, load_gain_noise):
    corrected = completed_retrievals(
        inversion_root / "corrected", state_index, "corrected", load_gain_noise
    )
    uncorrected = completed_retrievals(
        inversion_root / "uncorrected", state_index, "uncorrected",
        load_gain_noise,
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


def write_table(path, perturbations, truth, corrected, uncorrected,
                reference, include_gain_noise):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# Paired corrected and uncorrected terminal retrieval "
            f"displacements relative to {reference}.\n"
        )
        stream.write(
            "# State schema: %s\n" %
            corrected[perturbations[0]]["metadata"]["plot_schema_description"]
        )
        columns = (
            "# parameter units perturbation truth corrected_prior "
            "uncorrected_prior corrected_retrieved uncorrected_retrieved "
            f"corrected_minus_{reference} uncorrected_minus_{reference}"
        )
        if include_gain_noise:
            columns += (
                " corrected_G_delta_y uncorrected_G_delta_y"
                " corrected_xretr0_minus_truth_plus_G_delta_y"
                " uncorrected_xretr0_minus_truth_plus_G_delta_y"
            )
        stream.write(columns + "\n")
        for key, _, _, _, unit in PANELS:
            for perturbation in perturbations:
                corrected_value = corrected[perturbation]["final"][key]
                uncorrected_value = uncorrected[perturbation]["final"][key]
                corrected_prior = corrected[perturbation]["prior"][key]
                uncorrected_prior = uncorrected[perturbation]["prior"][key]
                corrected_reference = (
                    corrected_prior if reference == "prior" else truth[key]
                )
                uncorrected_reference = (
                    uncorrected_prior if reference == "prior" else truth[key]
                )
                line = (
                    f"{key} {unit.replace(' ', '_')} {perturbation:02d} "
                    f"{truth[key]:.12g} {corrected_prior:.12g} "
                    f"{uncorrected_prior:.12g} {corrected_value:.12g} "
                    f"{uncorrected_value:.12g} "
                    f"{corrected_value - corrected_reference:.12g} "
                    f"{uncorrected_value - uncorrected_reference:.12g}"
                )
                if include_gain_noise:
                    corrected_prediction = (
                        corrected[11]["final"][key] - truth[key]
                        + corrected[perturbation]["gain_noise"][key]
                    )
                    uncorrected_prediction = (
                        uncorrected[11]["final"][key] - truth[key]
                        + uncorrected[perturbation]["gain_noise"][key]
                    )
                    line += (
                        f" {corrected[perturbation]['gain_noise'][key]:.12g}"
                        f" {uncorrected[perturbation]['gain_noise'][key]:.12g}"
                        f" {corrected_prediction:.12g}"
                        f" {uncorrected_prediction:.12g}"
                    )
                stream.write(line + "\n")


def symmetric_limit(x_values, y_values):
    maximum = float(
        np.max(np.abs(np.concatenate([x_values, y_values, np.asarray([0.0])])))
    )
    if maximum == 0.0:
        maximum = 1.0
    return 1.16 * maximum


def main():
    args = parse_args()
    if args.gain_noise and args.reference != "truth":
        raise RuntimeError(
            "the linearized overlay is defined as x_retr,0 - x_truth + "
            "G @ delta_y; use --reference truth, or disable the overlay with "
            "--no-gain-prediction"
        )
    corrected, uncorrected, perturbations = matched_records(
        args.inversion_root, args.state, args.gain_noise
    )
    if args.gain_noise and 11 not in perturbations:
        raise RuntimeError(
            "the linearized overlay requires matched corrected and uncorrected "
            "perturbation-11 noiseless retrievals"
        )
    metadata = check_pair_metadata(corrected, uncorrected, perturbations)
    truth_row = table_row(args.truth_table, args.state)
    schema_description = metadata["plot_schema_description"]
    # Identity equality was already enforced by check_pair_metadata.  Load one
    # file again only to validate that the selected truth table belongs to the
    # same SIF convention as the retrieval campaign.
    with Dataset(corrected[perturbations[0]]["path"]) as dataset:
        schema = RetrievalPlotSchema.from_dataset(
            dataset, str(corrected[perturbations[0]]["path"])
        )
        schema.validate_truth_row(truth_row, str(args.truth_table))
    truth = truth_values(
        truth_row,
        aerosol_aod760(args.scene_components, metadata["aerosol_case"]),
        corrected[perturbations[0]]["prior"],
        schema=schema,
    )

    output = args.output or args.inversion_root / (
        f"state{args.state:03d}_corrected_vs_uncorrected_"
        f"{args.reference}_displacements.png"
    )
    table_output = args.table_output or output.with_suffix(".dat")
    write_table(
        table_output, perturbations, truth, corrected, uncorrected,
        args.reference, args.gain_noise,
    )

    tab10 = plt.get_cmap("tab10")
    colors = np.asarray([
        np.asarray((0.12, 0.12, 0.12, 1.0))
        if perturbation == 11 else tab10((perturbation - 1) % 10)
        for perturbation in perturbations
    ])
    figure, axes = plt.subplots(3, 7, figsize=(23.5, 10.5))
    for axis in axes.ravel():
        axis.set_visible(False)

    for key, title, row, column, unit in PANELS:
        axis = axes[row, column]
        axis.set_visible(True)
        x_values = np.asarray([
            corrected[perturbation]["final"][key]
            - (
                corrected[perturbation]["prior"][key]
                if args.reference == "prior" else truth[key]
            )
            for perturbation in perturbations
        ])
        y_values = np.asarray([
            uncorrected[perturbation]["final"][key]
            - (
                uncorrected[perturbation]["prior"][key]
                if args.reference == "prior" else truth[key]
            )
            for perturbation in perturbations
        ])
        if args.gain_noise:
            corrected_noiseless_error = corrected[11]["final"][key] - truth[key]
            uncorrected_noiseless_error = (
                uncorrected[11]["final"][key] - truth[key]
            )
            gain_x_values = np.asarray([
                corrected_noiseless_error
                + corrected[perturbation]["gain_noise"][key]
                for perturbation in perturbations
            ])
            gain_y_values = np.asarray([
                uncorrected_noiseless_error
                + uncorrected[perturbation]["gain_noise"][key]
                for perturbation in perturbations
            ])
            limit = symmetric_limit(
                np.concatenate([x_values, gain_x_values]),
                np.concatenate([y_values, gain_y_values]),
            )
        else:
            gain_x_values = gain_y_values = None
            limit = symmetric_limit(x_values, y_values)
        axis.plot(
            [-limit, limit], [-limit, limit], color="0.35", linestyle="--",
            linewidth=1.0, zorder=1,
        )
        axis.axhline(0.0, color="0.78", linewidth=0.8, zorder=0)
        axis.axvline(0.0, color="0.78", linewidth=0.8, zorder=0)
        noiseless = np.asarray(perturbations) == 11
        noisy = ~noiseless
        if np.any(noisy):
            axis.scatter(
                x_values[noisy], y_values[noisy], c=colors[noisy], marker="o",
                s=49, edgecolors="white", linewidths=0.7, zorder=3,
            )
        if np.any(noiseless):
            axis.scatter(
                x_values[noiseless], y_values[noiseless],
                c=colors[noiseless], marker="s", s=62, edgecolors="white",
                linewidths=0.8, zorder=5,
            )
        if args.gain_noise:
            if np.any(noisy):
                axis.scatter(
                    gain_x_values[noisy], gain_y_values[noisy], marker="o",
                    facecolors="none", edgecolors=colors[noisy], s=72,
                    linewidths=1.35, zorder=4,
                )
            if np.any(noiseless):
                axis.scatter(
                    gain_x_values[noiseless], gain_y_values[noiseless],
                    marker="s", facecolors="none", edgecolors=colors[noiseless],
                    s=86, linewidths=1.5, zorder=6,
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
            [0], [0], marker="s" if perturbation == 11 else "o",
            linestyle="none", markersize=7,
            markerfacecolor=color, markeredgecolor="white",
            label=(
                "Perturbation 11 (noiseless)" if perturbation == 11
                else f"Perturbation {perturbation:02d}"
            ),
        )
        for perturbation, color in zip(perturbations, colors)
    ]
    legend_handles.extend([
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=7,
            markerfacecolor="0.35", markeredgecolor="white",
            label=f"Filled: $x_{{retr}}-x_{{{args.reference}}}$",
        ),
    ])
    if args.gain_noise:
        legend_handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="none", markersize=7,
                markerfacecolor="none", markeredgecolor="0.25",
                markeredgewidth=1.35,
                label=(
                    "Open: $x_{retr,0}-x_{truth}+G\\,\\Delta y$"
                ),
            )
        )
    figure.legend(
        handles=legend_handles, loc="lower center",
        ncol=min(11, len(legend_handles)),
        frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.005),
    )
    coordinate_label = f"$x_{{retr}}-x_{{{args.reference}}}$"
    if args.gain_noise:
        coordinate_label += " or $x_{retr,0}-x_{truth}+G\\,\\Delta y$"
    figure.text(
        0.5, 0.055,
        f"Corrected: {coordinate_label} (panel units)",
        ha="center", va="center", fontsize=13,
    )
    figure.text(
        0.018, 0.49,
        f"Uncorrected: {coordinate_label} (panel units)",
        ha="center", va="center", rotation="vertical", fontsize=13,
    )
    figure.suptitle(
        f"State {args.state:03d}: corrected versus uncorrected terminal-state "
        f"{args.reference} displacements"
        f"{' and linearized noise predictions' if args.gain_noise else ''}\n"
        f"{metadata['surface']}, {metadata['aerosol_case']}, "
        f"{co2_truth_label(truth_row)}; "
        f"{len(perturbations)} matched perturbations\n"
        f"{schema_description}",
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
    print(
        "linearized gain/noise overlay:",
        "enabled" if args.gain_noise else "disabled",
    )
    print("matched perturbations:", " ".join(f"{p:02d}" for p in perturbations))


if __name__ == "__main__":
    main()
