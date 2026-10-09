#!/usr/bin/env python3
"""Create four trial views of corrected/uncorrected retrieval spectral fits.

Only states having complete, paired corrected and uncorrected retrievals for
all perturbations 01--11 are admitted to the ensemble figures.  Retrieval
files store ``final_residual = F(x) - y``; every residual shown here uses the
more intuitive sign ``y - F(x)`` and the script verifies that conversion.
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
from netCDF4 import Dataset

from co2_plot_metadata import (
    bottom_layer_mapping_from_rows,
    co2_case_label,
)
from retrieval_plot_schema import RetrievalPlotSchema


BANDS = ((1, "O$_2$ A band"), (2, "Weak CO$_2$ band"),
         (3, "Strong CO$_2$ band"))
CLASSES = ("corrected", "uncorrected")
CLASS_COLOR = {"corrected": "#159447", "uncorrected": "#d43f3a"}
SURFACE_COLOR = {
    "urban": "#4c78a8", "rural": "#54a24b",
    "desert": "#e19c24", "forest": "#8b5fbf",
}
FILE_RE = re.compile(r"retrieval_state(\d{3})_perturbation(\d{2})\.nc$")


def parse_args():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inversion-dir", type=Path, default=here)
    parser.add_argument("--truth-table", type=Path,
                        default=here.parent / "truth_map" / "true_states.dat")
    parser.add_argument("--state", type=int, default=33,
                        help="state used for the fit card and RRS budget")
    parser.add_argument("--output-dir", type=Path,
                        default=here / "spectral_fit_visualization_trial")
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def scalar_attr(dataset, name, default=None):
    return dataset.getncattr(name) if name in dataset.ncattrs() else default


def read_record(path, expected_class, state, perturbation):
    """Read only compact measurement-space diagnostics from one retrieval."""
    with Dataset(str(path), "r") as dataset:
        if int(scalar_attr(dataset, "retrieval_complete", 0)) != 1:
            raise RuntimeError("retrieval is not complete: %s" % path)
        actual_class = str(scalar_attr(dataset, "measurement_class", ""))
        if actual_class != expected_class:
            raise RuntimeError("class mismatch in %s" % path)
        if int(scalar_attr(dataset, "truth_state_index", -1)) != state:
            raise RuntimeError("truth-state mismatch in %s" % path)
        if int(scalar_attr(dataset, "perturbation_index", -1)) != perturbation:
            raise RuntimeError("perturbation mismatch in %s" % path)

        schema = RetrievalPlotSchema.from_dataset(dataset, str(path))
        record = {
            "path": path,
            "wavelength": np.asarray(dataset["wavelength"][:], dtype=float),
            "band": np.asarray(dataset["band_index"][:], dtype=int),
            "y": np.asarray(dataset["measurement_perturbed"][:], dtype=float),
            "y0": np.asarray(dataset["measurement_noiseless"][:], dtype=float),
            "sigma": np.asarray(dataset["noise_standard_deviation"][:],
                                dtype=float),
            "forward": np.asarray(dataset["final_forward_model"][:],
                                  dtype=float),
            "stored_residual": np.asarray(dataset["final_residual"][:],
                                          dtype=float),
            "chi2": np.asarray(
                dataset["final_band_reduced_chi_squared"][:], dtype=float),
            "surface": str(scalar_attr(dataset, "surface", "unknown")),
            "aerosol": str(scalar_attr(dataset, "aerosol_case", "unknown")),
            "xco2": float(scalar_attr(dataset, "truth_xco2_ppm", np.nan)),
            "bottom_co2": float(
                scalar_attr(dataset, "truth_bottom_co2_ppm", np.nan)
            ),
            "converged": int(scalar_attr(dataset, "converged", 0)) == 1,
            "fit_ok": int(scalar_attr(dataset, "fit_quality_ok", 0)) == 1,
            "schema": schema,
            "schema_identity": schema.identity,
            "schema_description": schema.description(),
        }

    shape = record["wavelength"].shape
    names = ("band", "y", "y0", "sigma", "forward", "stored_residual")
    if any(record[name].shape != shape for name in names):
        raise RuntimeError("measurement arrays do not share one shape: %s" % path)
    if record["chi2"].shape != (3,):
        raise RuntimeError("expected three band chi-squared values: %s" % path)
    numeric = [record[name] for name in
               ("wavelength", "y", "y0", "sigma", "forward",
                "stored_residual", "chi2")]
    if not all(np.all(np.isfinite(value)) for value in numeric):
        raise RuntimeError("non-finite retrieval diagnostics: %s" % path)
    if not np.all(record["sigma"] > 0.0):
        raise RuntimeError("non-positive noise standard deviation: %s" % path)
    expected = record["forward"] - record["y"]
    if not np.allclose(record["stored_residual"], expected,
                       rtol=2e-12, atol=2e-12):
        raise RuntimeError("stored residual is not F(x)-y: %s" % path)

    # All plots below deliberately use y-F(x), not the stored F(x)-y sign.
    record["residual"] = -record["stored_residual"]
    record["z"] = record["residual"] / record["sigma"]
    return record


def index_complete_files(inversion_dir):
    indexed = {retrieval_class: {} for retrieval_class in CLASSES}
    for retrieval_class in CLASSES:
        for path in sorted((inversion_dir / retrieval_class).glob(
                "retrieval_state*_perturbation*.nc")):
            match = FILE_RE.match(path.name)
            if match is None:
                continue
            state, perturbation = map(int, match.groups())
            try:
                with Dataset(str(path), "r") as dataset:
                    complete = int(scalar_attr(
                        dataset, "retrieval_complete", 0)) == 1
            except Exception:
                complete = False
            if complete:
                indexed[retrieval_class][(state, perturbation)] = path
    required = set(range(1, 12))
    states = []
    candidates = sorted(set(
        state for retrieval_class in CLASSES
        for state, _ in indexed[retrieval_class].keys()))
    for state in candidates:
        if all(set(p for s, p in indexed[retrieval_class].keys() if s == state)
               >= required for retrieval_class in CLASSES):
            states.append(state)
    return indexed, states


def read_truth_table(path):
    columns = None
    rows = {}
    with path.open("r") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith("# index "):
                columns = text[2:].split()
                continue
            if text.startswith("#"):
                continue
            if columns is None:
                raise RuntimeError("truth-table column header was not found")
            values = text.split()
            if len(values) != len(columns):
                raise RuntimeError("malformed truth-table row: %s" % text)
            row = dict(zip(columns, values))
            row["index"] = int(row["index"])
            row["xco2_ppm"] = float(row["xco2_ppm"])
            if "bottom_co2_ppm" in row:
                row["bottom_co2_ppm"] = float(row["bottom_co2_ppm"])
            rows[row["index"]] = row
    return rows


def load_ensemble(indexed, states, truth_rows):
    records = {}
    reference_wave = None
    reference_band = None
    reference_schema = None
    for state in states:
        if state not in truth_rows:
            raise RuntimeError("state %03d absent from truth table" % state)
        for perturbation in range(1, 12):
            pair = {}
            for retrieval_class in CLASSES:
                path = indexed[retrieval_class][(state, perturbation)]
                pair[retrieval_class] = read_record(
                    path, retrieval_class, state, perturbation)
            for retrieval_class, record in pair.items():
                if reference_schema is None:
                    reference_schema = record["schema_identity"]
                elif record["schema_identity"] != reference_schema:
                    raise RuntimeError(
                        "retrieval-state schema changed within ensemble: %s" %
                        record["path"]
                    )
                if reference_wave is None:
                    reference_wave = record["wavelength"].copy()
                    reference_band = record["band"].copy()
                if not np.array_equal(record["wavelength"], reference_wave):
                    raise RuntimeError("non-common wavelength grid: %s" %
                                       record["path"])
                if not np.array_equal(record["band"], reference_band):
                    raise RuntimeError("non-common band grid: %s" % record["path"])
            if not np.array_equal(pair["corrected"]["wavelength"],
                                  pair["uncorrected"]["wavelength"]):
                raise RuntimeError("corrected/uncorrected grid mismatch")
            truth = truth_rows[state]
            for retrieval_class in CLASSES:
                record = pair[retrieval_class]
                record["schema"].validate_truth_row(
                    truth, "truth state %03d" % state
                )
                if record["surface"] != truth["surface"]:
                    raise RuntimeError("surface metadata mismatch: %s" %
                                       record["path"])
                if record["aerosol"] != truth["aerosol_case"]:
                    raise RuntimeError("aerosol metadata mismatch: %s" %
                                       record["path"])
                if abs(record["xco2"] - truth["xco2_ppm"]) > 1e-10:
                    raise RuntimeError("XCO2 metadata mismatch: %s" %
                                       record["path"])
                if "bottom_co2_ppm" in truth and (
                        not np.isfinite(record["bottom_co2"]) or
                        abs(record["bottom_co2"] -
                            truth["bottom_co2_ppm"]) > 1e-10):
                    raise RuntimeError("bottom-layer CO2 metadata mismatch: %s" %
                                       record["path"])
            records[(state, perturbation)] = pair
    return records, reference_wave, reference_band


def percentile_stack(records, retrieval_class, key, state, perturbations,
                     selected):
    values = np.vstack([
        records[(state, perturbation)][retrieval_class][key][selected]
        for perturbation in perturbations
    ])
    return (np.percentile(values, 16, axis=0),
            np.median(values, axis=0),
            np.percentile(values, 84, axis=0))


def make_row_geometry(number_states):
    # Median noisy ensemble gets 78% of each state row; p11 gets a thin 22%.
    edges = []
    for index in range(number_states):
        edges.extend([index, index + 0.78])
    edges.append(number_states)
    return np.asarray(edges, dtype=float)


def co2_truth_description(row):
    return co2_case_label(row.get("bottom_co2_ppm"), row["xco2_ppm"])


def co2_color_coordinate(states, truth_rows):
    bottom_campaign = all(
        "bottom_co2_ppm" in truth_rows[state] for state in states
    )
    field = "bottom_co2_ppm" if bottom_campaign else "xco2_ppm"
    label = "bottom-layer CO$_2$" if bottom_campaign else "XCO$_2$"
    values = [float(truth_rows[state][field]) for state in states]
    return field, label, values


def draw_annotation_strips(axis, states, truth_rows, row_edges):
    surfaces = list(SURFACE_COLOR.keys())
    aerosol_color = {"none": "#f2f2f2", "aod760_0p28": "#4d4d4d"}
    sif_color = {
        "off": "#f2f2f2",
        "angular_integral760_0p5": "#f2c14e",
    }
    co2_field, _, co2_values = co2_color_coordinate(states, truth_rows)
    co2_norm = Normalize(min(co2_values), max(co2_values))
    co2_map = plt.get_cmap("viridis")
    for state_index, state in enumerate(states):
        row = truth_rows[state]
        colors = [SURFACE_COLOR.get(row["surface"], "0.6"),
                  aerosol_color.get(row["aerosol_case"], "0.6"),
                  sif_color.get(row["sif_case"], "0.6"),
                  co2_map(co2_norm(float(row[co2_field])))]
        for column, color in enumerate(colors):
            axis.add_patch(Rectangle((column, state_index), 1.0, 1.0,
                                     facecolor=color, edgecolor="white",
                                     linewidth=0.25))
        # Delineate the small p11 sub-row in the same metadata strip.
        axis.plot([0, 4], [row_edges[2 * state_index + 1]] * 2,
                  color="white", linewidth=0.25)
    axis.set_xlim(0, 4)
    axis.set_ylim(0, len(states))
    axis.invert_yaxis()
    axis.set_xticks(np.arange(4) + 0.5)
    axis.set_xticklabels(("surface", "aer.", "SIF", "CO$_2$"),
                         rotation=55, ha="right", fontsize=7)
    axis.set_yticks(np.arange(len(states)) + 0.39)
    axis.set_yticklabels(["%03d" % state for state in states], fontsize=7)
    axis.set_ylabel("truth state  (main=median p01--10; thin=p11)")
    axis.tick_params(length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)


def make_atlas(records, states, truth_rows, wavelength, band_index,
               output_dir, dpi):
    paths = []
    cmap = plt.get_cmap("RdBu_r")
    row_edges = make_row_geometry(len(states))
    for iband, band_label in BANDS:
        selected = band_index == iband
        x = wavelength[selected]
        panels = []
        flags = []
        for retrieval_class in CLASSES:
            rows = []
            class_flags = []
            for state in states:
                noisy = np.vstack([
                    records[(state, p)][retrieval_class]["z"][selected]
                    for p in range(1, 11)])
                rows.extend([np.median(noisy, axis=0),
                             records[(state, 11)][retrieval_class]["z"][selected]])
                class_flags.extend([
                    not all(records[(state, p)][retrieval_class]["converged"]
                            for p in range(1, 11)),
                    not records[(state, 11)][retrieval_class]["converged"],
                ])
            panels.append(np.vstack(rows))
            flags.append(class_flags)
        panels.append(panels[1] - panels[0])
        flags.append([a or b for a, b in zip(flags[0], flags[1])])

        fig = plt.figure(figsize=(18.5, max(9.5, 0.43 * len(states))))
        grid = fig.add_gridspec(1, 4, width_ratios=(1.0, 4.2, 4.2, 4.2),
                                left=0.055, right=0.965, bottom=0.09, top=0.82,
                                wspace=0.08)
        annotation_axis = fig.add_subplot(grid[0, 0])
        draw_annotation_strips(annotation_axis, states, truth_rows, row_edges)
        axes = [fig.add_subplot(grid[0, index]) for index in range(1, 4)]
        titles = ("Corrected", "Uncorrected", "Uncorrected - corrected")
        limits = (3.0, 3.0, 5.0)
        meshes = []
        for axis, values, title, limit, bad_rows in zip(
                axes, panels, titles, limits, flags):
            mesh = axis.pcolormesh(x, row_edges, values, cmap=cmap,
                                   norm=Normalize(-limit, limit),
                                   shading="auto", rasterized=True)
            meshes.append(mesh)
            axis.set_title("%s  (limit $\\pm$%g)" % (title, limit))
            axis.set_ylim(0, len(states))
            axis.invert_yaxis()
            axis.set_yticks([])
            axis.set_xlabel("Wavelength (nm)")
            axis.margins(x=0)
            dx = x[-1] - x[0]
            for row_number, failed in enumerate(bad_rows):
                if failed:
                    y = 0.5 * (row_edges[row_number] + row_edges[row_number + 1])
                    axis.text(x[0] + 0.008 * dx, y, "×", color="black",
                              fontsize=8, ha="left", va="center", weight="bold")
        colorbar1 = fig.colorbar(meshes[0], ax=axes[:2], orientation="horizontal",
                                fraction=0.025, pad=0.07, aspect=55)
        colorbar1.set_label("median normalized residual $(y-F)/\\sigma$")
        colorbar2 = fig.colorbar(meshes[2], ax=axes[2], orientation="horizontal",
                                fraction=0.025, pad=0.07, aspect=28)
        colorbar2.set_label("difference in normalized residual")
        fig.suptitle(
            "%s: all-state spectral residual atlas (%d fully paired states)\n"
            "× marks a nonconverged retrieval; color limits clip outliers" %
            (band_label, len(states)), fontsize=14)
        _, co2_label, co2_values = co2_color_coordinate(states, truth_rows)
        mapping_note = bottom_layer_mapping_from_rows(truth_rows.values())
        strip_note = (
            "strip key — surface: blue urban, green rural, gold desert, purple "
            "forest; aerosol: dark=present/light=none; SIF: gold=on/light=off; "
            "%s: purple %g $\\rightarrow$ yellow %g ppm" % (
                co2_label, min(co2_values), max(co2_values)
            )
        )
        if mapping_note:
            strip_note += "\n" + mapping_note
        fig.text(
            0.055, 0.865,
            strip_note, fontsize=7.5, ha="left", va="center",
            linespacing=1.3)
        path = output_dir / ("01_residual_atlas_band%d.png" % iband)
        fig.savefig(str(path), dpi=dpi)
        plt.close(fig)
        paths.append(path)
    return paths


def make_fit_card(records, state, truth_rows, wavelength, band_index,
                  output_dir, dpi):
    fig, axes = plt.subplots(3, 3, figsize=(19, 12.3))
    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.15, top=0.89,
                        hspace=0.28, wspace=0.18)
    styles = {"corrected": "-", "uncorrected": "--"}
    for row_index, (iband, label) in enumerate(BANDS):
        selected = band_index == iband
        x = wavelength[selected]
        fit_axis, residual_axis, contrast_axis = axes[row_index]
        for retrieval_class in CLASSES:
            color = CLASS_COLOR[retrieval_class]
            low, median, high = percentile_stack(
                records, retrieval_class, "forward", state, range(1, 11), selected)
            fit_axis.fill_between(x, low, high, color=color, alpha=0.20,
                                  linewidth=0)
            fit_axis.plot(x, median, color=color, linewidth=1.25,
                          label="%s median $F$, p01--10" % retrieval_class.title())
            p11 = records[(state, 11)][retrieval_class]
            fit_axis.plot(x, p11["y"][selected], color="0.55",
                          linestyle=styles[retrieval_class], linewidth=0.75,
                          label="%s noiseless $y$" % retrieval_class.title())
            fit_axis.plot(x, p11["forward"][selected], color="black",
                          linestyle=styles[retrieval_class], linewidth=1.65,
                          label="%s p11 $F$" % retrieval_class.title())

            low, median, high = percentile_stack(
                records, retrieval_class, "z", state, range(1, 11), selected)
            residual_axis.fill_between(x, low, high, color=color, alpha=0.20,
                                       linewidth=0)
            residual_axis.plot(x, median, color=color, linewidth=1.0)
            residual_axis.plot(x, p11["z"][selected], color="black",
                               linestyle=styles[retrieval_class], linewidth=1.55)

        noisy_contrasts = np.vstack([
            records[(state, p)]["uncorrected"]["z"][selected] -
            records[(state, p)]["corrected"]["z"][selected]
            for p in range(1, 11)])
        low, median, high = np.percentile(noisy_contrasts, (16, 50, 84), axis=0)
        p11_contrast = (records[(state, 11)]["uncorrected"]["z"][selected] -
                        records[(state, 11)]["corrected"]["z"][selected])
        contrast_axis.fill_between(x, low, high, color="#7b4ab5", alpha=0.24,
                                   linewidth=0)
        contrast_axis.plot(x, median, color="#7b4ab5", linewidth=1.1,
                           label="median p01--10")
        contrast_axis.plot(x, p11_contrast, color="black", linewidth=1.65,
                           label="p11 (noiseless)")

        fit_axis.set_title("%s: measurement and terminal fit" % label, pad=8)
        residual_axis.set_title("Normalized residual $(y-F)/\\sigma$", pad=8)
        contrast_axis.set_title(
            "Uncorrected - corrected normalized residual", pad=8)
        fit_axis.set_ylabel("Radiance\n(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)")
        residual_axis.set_ylabel("$z$")
        contrast_axis.set_ylabel("$\\Delta z$")
        residual_axis.axhspan(-1, 1, color="0.85", alpha=0.42, zorder=0)
        residual_axis.axhline(0, color="0.4", linewidth=0.6)
        contrast_axis.axhline(0, color="0.4", linewidth=0.6)
        for axis in axes[row_index]:
            axis.grid(alpha=0.16)
            axis.margins(x=0)
            if row_index == 2:
                axis.set_xlabel("Wavelength (nm)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    contrast_handles, contrast_labels = axes[0, 2].get_legend_handles_labels()
    handles.extend(contrast_handles)
    labels.extend(["Contrast " + label for label in contrast_labels])
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.015),
               ncol=4, frameon=False, fontsize=8, columnspacing=1.5)
    fig.suptitle(
        "State %03d (%s) spectral fit card: black = noiseless p11; color = "
        "noisy median with 16--84%% ribbon\n%s" % (
            state, co2_truth_description(truth_rows[state]),
            records[(state, 11)]["corrected"]["schema_description"],
        ), fontsize=14, y=0.975)
    path = output_dir / ("02_state%03d_spectral_fit_card.png" % state)
    fig.savefig(str(path), dpi=dpi)
    plt.close(fig)
    return path


def make_rrs_budget(records, state, truth_rows, wavelength, band_index,
                    output_dir, dpi):
    corrected = records[(state, 11)]["corrected"]
    uncorrected = records[(state, 11)]["uncorrected"]
    fig, axes = plt.subplots(3, 1, figsize=(14.5, 10), constrained_layout=True)
    maximum_closure = 0.0
    for axis, (iband, label) in zip(axes, BANDS):
        selected = band_index == iband
        x = wavelength[selected]
        # One pooled sigma puts both measurement classes on a common scale.
        sigma = np.sqrt(0.5 * (corrected["sigma"][selected] ** 2 +
                              uncorrected["sigma"][selected] ** 2))
        truth_discrepancy = ((uncorrected["y0"][selected] -
                              corrected["y0"][selected]) / sigma)
        state_motion = ((uncorrected["forward"][selected] -
                         corrected["forward"][selected]) / sigma)
        terminal_discrepancy = (
            uncorrected["residual"][selected] -
            corrected["residual"][selected]) / sigma
        closure = truth_discrepancy - state_motion - terminal_discrepancy
        maximum_closure = max(maximum_closure, float(np.max(np.abs(closure))))
        axis.plot(x, truth_discrepancy, color="black", linewidth=1.45,
                  label="Truth discrepancy $(y_u-y_c)/\\sigma_{common}$")
        axis.plot(x, state_motion, color="#7246a5", linewidth=1.15,
                  label="Absorbed by state motion $(F_u-F_c)/\\sigma_{common}$")
        axis.plot(x, terminal_discrepancy, color="#e1812c", linewidth=1.15,
                  label="Terminal discrepancy $(r_u-r_c)/\\sigma_{common}$")
        axis.axhline(0, color="0.45", linewidth=0.6)
        axis.axhspan(-1, 1, color="0.86", alpha=0.35, zorder=0)
        axis.set_title(label)
        axis.set_ylabel("Common-noise units")
        axis.set_xlabel("Wavelength (nm)")
        axis.grid(alpha=0.16)
        axis.margins(x=0)
    axes[0].legend(frameon=False, fontsize=8, ncol=3)
    fig.suptitle(
        "State %03d (%s) noiseless RRS discrepancy budget\n"
        "$r_u-r_c=(y_u-y_c)-(F_u-F_c)$; $\\sigma_{common}="
        "\\sqrt{(\\sigma_u^2+\\sigma_c^2)/2}$; max closure %.2e\n%s" %
        (state, co2_truth_description(truth_rows[state]), maximum_closure,
         corrected["schema_description"]),
        fontsize=13)
    path = output_dir / ("03_state%03d_rrs_absorption_budget.png" % state)
    fig.savefig(str(path), dpi=dpi)
    plt.close(fig)
    return path, maximum_closure


def make_quality_summary(records, states, truth_rows, output_dir, dpi):
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 7.0))
    mapping_note = bottom_layer_mapping_from_rows(truth_rows.values())
    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.22,
                        top=0.80 if mapping_note else 0.84,
                        wspace=0.19)
    all_values = []
    for state in states:
        for perturbation in range(1, 12):
            for retrieval_class in CLASSES:
                all_values.extend(np.sqrt(np.maximum(
                    records[(state, perturbation)][retrieval_class]["chi2"], 0)))
    positive_values = [value for value in all_values if value > 0.0]
    lower = max(0.01, float(np.min(positive_values)) * 0.78)
    upper = max(1.25, float(np.max(all_values)) * 1.08)
    for band_offset, (axis, (_, label)) in enumerate(zip(axes, BANDS)):
        for state in states:
            truth = truth_rows[state]
            color = SURFACE_COLOR.get(truth["surface"], "0.4")
            aerosol = truth["aerosol_case"] != "none"
            for perturbation in range(1, 12):
                corr = records[(state, perturbation)]["corrected"]
                uncorr = records[(state, perturbation)]["uncorrected"]
                x = np.sqrt(max(corr["chi2"][band_offset], 0.0))
                y = np.sqrt(max(uncorr["chi2"][band_offset], 0.0))
                marker = "s" if perturbation == 11 else "o"
                face = color if aerosol else "none"
                axis.scatter(x, y, s=31 if perturbation == 11 else 18,
                             marker=marker, facecolors=face, edgecolors=color,
                             linewidths=0.75, alpha=0.84, zorder=3)
                if not (corr["converged"] and uncorr["converged"]):
                    axis.scatter(x, y, s=31, marker="x", color="black",
                                 linewidths=0.75, zorder=4)
        axis.plot([lower, upper], [lower, upper], color="0.35", linewidth=0.8,
                  linestyle="--", zorder=1)
        threshold = np.sqrt(1.4)
        axis.axvline(threshold, color="0.65", linewidth=0.7, linestyle=":")
        axis.axhline(threshold, color="0.65", linewidth=0.7, linestyle=":")
        axis.axvline(1.0, color="0.80", linewidth=0.6)
        axis.axhline(1.0, color="0.80", linewidth=0.6)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlim(lower, upper)
        axis.set_ylim(lower, upper)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(label)
        axis.set_xlabel("Corrected $\\sqrt{\\chi^2_{red}}$")
        axis.set_ylabel("Uncorrected $\\sqrt{\\chi^2_{red}}$")
        axis.grid(alpha=0.13)
    legend = [Line2D([], [], marker="o", linestyle="none", markerfacecolor="0.4",
                     markeredgecolor="0.4", label="p01--10 (noisy)"),
              Line2D([], [], marker="s", linestyle="none", markerfacecolor="0.4",
                     markeredgecolor="0.4", label="p11 (noiseless)"),
              Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
                     markeredgecolor="0.3", label="no aerosol"),
              Line2D([], [], marker="o", linestyle="none", markerfacecolor="0.3",
                     markeredgecolor="0.3", label="with aerosol"),
              Line2D([], [], marker="x", linestyle="none", color="black",
                     label="one/both nonconverged")]
    legend.extend([Patch(facecolor=color, label=surface.title())
                   for surface, color in SURFACE_COLOR.items()])
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.015),
               ncol=5, frameon=False, fontsize=8, columnspacing=1.6)
    title = (
        "Compact terminal fit quality (%d fully paired states); dotted = "
        "$\\sqrt{1.4}$ fit-quality threshold" % len(states)
    )
    if mapping_note:
        title += "\n" + mapping_note
    fig.suptitle(title, fontsize=14)
    path = output_dir / "04_compact_fit_quality_summary.png"
    fig.savefig(str(path), dpi=dpi)
    plt.close(fig)
    return path


def write_manifest(path, states, state, outputs, maximum_closure,
                   schema_description):
    with path.open("w") as handle:
        handle.write("Spectral-fit visualization trial\n")
        handle.write("================================\n")
        handle.write("state-specific example: %03d\n" % state)
        handle.write("fully paired states (%d): %s\n" %
                     (len(states), " ".join("%03d" % x for x in states)))
        handle.write("perturbations required: 01--11 in both classes\n")
        handle.write("displayed residual convention: y - F(x)\n")
        handle.write("stored residual convention verified: F(x) - y\n")
        handle.write("retrieval state schema: %s\n" % schema_description)
        handle.write("RRS budget maximum algebraic closure: %.12e\n" %
                     maximum_closure)
        handle.write("outputs:\n")
        for output in outputs:
            handle.write("  %s\n" % output.name)


def main():
    args = parse_args()
    indexed, states = index_complete_files(args.inversion_dir)
    if not states:
        raise RuntimeError("no fully paired state has perturbations 01--11")
    if args.state not in states:
        raise RuntimeError("state %03d is not fully paired and complete" % args.state)
    truth_rows = read_truth_table(args.truth_table)
    records, wavelength, band_index = load_ensemble(indexed, states, truth_rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    outputs.extend(make_atlas(records, states, truth_rows, wavelength, band_index,
                              args.output_dir, args.dpi))
    outputs.append(make_fit_card(
        records, args.state, truth_rows, wavelength, band_index,
        args.output_dir, args.dpi
    ))
    budget_path, closure = make_rrs_budget(
        records, args.state, truth_rows, wavelength, band_index,
        args.output_dir, args.dpi
    )
    outputs.append(budget_path)
    outputs.append(make_quality_summary(records, states, truth_rows,
                                        args.output_dir, args.dpi))
    manifest = args.output_dir / "trial_manifest.txt"
    write_manifest(
        manifest, states, args.state, outputs, closure,
        records[(args.state, 11)]["corrected"]["schema_description"],
    )
    outputs.append(manifest)
    print("Validated %d fully paired states: %s" %
          (len(states), ", ".join("%03d" % state for state in states)))
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
