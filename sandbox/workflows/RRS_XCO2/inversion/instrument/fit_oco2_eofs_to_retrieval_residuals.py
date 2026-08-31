#!/usr/bin/env python3
"""Fit the three OCO-2 EOF components to mean synthetic retrieval residuals.

The canonical land EOF waveforms are stored in ``l2_oco_eof.h5`` (three EOFs
per band and eight footprints per EOF).  The four nadir L1B orbit files used
by ``OCORaman/src/OCOPlots/oco_gain.jl`` provide the detector-sample wavelength
dispersion.  This script constructs one representative waveform per EOF by:

1. mapping each of the 32 orbit/footprint waveforms onto the synthetic grid;
2. aligning the arbitrary footprint EOF signs to footprint 1; and
3. averaging within each orbit and then equally over the four orbits.

For each truth state having matched corrected and uncorrected retrievals, the
terminal ``F(x)-y`` residual is averaged over the matched noise perturbations.
The three EOF columns are then fitted jointly, with inverse-noise-variance
weighting, independently in each band and retrieval class.  Reported EOF
coefficients multiply the dimensionless stored EOF waveforms and therefore
have residual-radiance units.  They are the amplitudes that would be *added to
the forward model*; the post-fit residual is therefore ``r + E*c``.
"""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


HERE = Path(__file__).resolve().parent
WORKFLOW_ROOT = HERE.parents[1]
DATA_ROOT = Path(os.environ.get("RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
INVERSION_ROOT = DATA_ROOT / "inversion"
DEFAULT_OCO_ROOT = Path("/home/sanghavi/data/OCO2_jacobians")
OUTPUT_ROOT = INVERSION_ROOT / "eof_residual_fits"

L1B_FILES = (
    "oco2_L1bScND_37214a_210630_B11006r_220326081005.h5",
    "oco2_L1bScND_42610a_220706_B11008r_220831171738.h5",
    "oco2_L1bScND_41856a_220515_B11008r_220804213809.h5",
    "oco2_L1bScND_45442a_230116_B11008r_230209231005.h5",
)
EOF_FILE = "l2_oco_eof.h5"

MIN_ORBIT_COVERAGE = 1.0
CLASSES = ("corrected", "uncorrected")


@dataclass(frozen=True)
class Band:
    index: int
    key: str
    title: str
    lower_um: float
    upper_um: float


BANDS = (
    Band(1, "o2a", "O$_2$ A band", 0.74, 0.80),
    Band(2, "weak_co2", "Weak CO$_2$ band", 1.55, 1.65),
    Band(3, "strong_co2", "Strong CO$_2$ band", 2.00, 2.10),
)


def as_float(variable, key=slice(None)):
    values = variable[key]
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    return np.asarray(values, dtype=np.float64)


def completed_retrievals():
    """Return complete retrieval paths indexed by class, state, perturbation."""
    result = {retrieval_class: defaultdict(dict) for retrieval_class in CLASSES}
    pattern = re.compile(r"retrieval_state(\d{3})_perturbation(\d{2})\.nc$")
    for retrieval_class in CLASSES:
        for path in sorted((INVERSION_ROOT / retrieval_class).glob(
            "retrieval_state*_perturbation*.nc"
        )):
            match = pattern.match(path.name)
            if match is None:
                continue
            try:
                with Dataset(path) as dataset:
                    if int(dataset.getncattr("retrieval_complete")) != 1:
                        continue
                    state = int(dataset.getncattr("truth_state_index"))
                    perturbation = int(dataset.getncattr("perturbation_index"))
            except (OSError, RuntimeError, AttributeError):
                continue
            result[retrieval_class][state][perturbation] = path
    return result


def matched_state_paths(retrievals):
    states = sorted(
        set(retrievals["corrected"]) & set(retrievals["uncorrected"])
    )
    result = {}
    for state in states:
        perturbations = sorted(
            set(retrievals["corrected"][state])
            & set(retrievals["uncorrected"][state])
        )
        if perturbations:
            result[state] = {
                retrieval_class: [
                    retrievals[retrieval_class][state][perturbation]
                    for perturbation in perturbations
                ]
                for retrieval_class in CLASSES
            }
            result[state]["perturbations"] = perturbations
    return result


def target_grid(first_retrieval):
    with Dataset(first_retrieval) as dataset:
        wavelength = as_float(dataset["wavelength"])
        band_index = np.asarray(dataset["band_index"][:], dtype=np.int8)
    grids = {}
    for band in BANDS:
        selected = band_index == band.index
        grid = wavelength[selected]
        if not np.all(np.diff(grid) > 0):
            raise RuntimeError(f"synthetic {band.key} grid is not increasing")
        grids[band.key] = grid
    return wavelength, band_index, grids


def canonical_waveforms(oco_root):
    path = oco_root / EOF_FILE
    waveforms = {}
    with Dataset(path) as dataset:
        for band in BANDS:
            components = []
            for order in range(1, 4):
                variable = dataset[
                    "/Instrument/EmpiricalOrthogonalFunction/Land/"
                    f"EOF_{order}_waveform_{band.index}"
                ]
                values = as_float(variable)
                if values.shape != (8, 1016):
                    raise RuntimeError(
                        f"unexpected {band.key} EOF shape {values.shape}"
                    )
                components.append(values)
            waveforms[band.key] = np.stack(components, axis=-1)
    return waveforms


def footprint_signs(waveforms):
    """Align footprint EOF signs to footprint 1 without changing magnitudes."""
    result = {}
    for band in BANDS:
        values = waveforms[band.key]
        signs = np.ones((8, 3), dtype=np.float64)
        for footprint in range(8):
            for order in range(3):
                dot = np.dot(values[0, :, order], values[footprint, :, order])
                if dot < 0:
                    signs[footprint, order] = -1.0
        result[band.key] = signs
    return result


def representative_eof_basis(oco_root, grids):
    """Map the stored EOF waveforms to wavelength and average 4x8 sources."""
    waveforms = canonical_waveforms(oco_root)
    signs = footprint_signs(waveforms)
    orbit_means = {
        band.key: np.full((len(L1B_FILES), grids[band.key].size, 3), np.nan)
        for band in BANDS
    }
    orbit_coverage = {
        band.key: np.zeros((len(L1B_FILES), grids[band.key].size, 3))
        for band in BANDS
    }
    samples = np.arange(1.0, 1017.0)

    for orbit_index, filename in enumerate(L1B_FILES):
        path = oco_root / filename
        print(f"EOF dispersion orbit {orbit_index + 1}/{len(L1B_FILES)}: {path}", flush=True)
        with Dataset(path) as dataset:
            dispersion = as_float(
                dataset["/InstrumentHeader/dispersion_coef_samp"]
            )
            if dispersion.shape != (3, 8, 6):
                raise RuntimeError(
                    f"unexpected dispersion shape {dispersion.shape} in {path}"
                )
            for band in BANDS:
                target = grids[band.key]
                sums = np.zeros((target.size, 3))
                counts = np.zeros((target.size, 3), dtype=np.int16)
                for footprint in range(8):
                    coefficients = dispersion[band.index - 1, footprint, :]
                    wavelength_nm = 1000.0 * np.polynomial.polynomial.polyval(
                        samples, coefficients
                    )
                    if not np.all(np.diff(wavelength_nm) > 0):
                        raise RuntimeError(
                            f"non-increasing dispersion for {path}, "
                            f"{band.key}, footprint {footprint + 1}"
                        )
                    for component in range(3):
                        waveform = (
                            waveforms[band.key][footprint, :, component]
                            * signs[band.key][footprint, component]
                        )
                        interpolated = np.interp(
                            target,
                            wavelength_nm,
                            waveform,
                            left=np.nan,
                            right=np.nan,
                        )
                        finite = np.isfinite(interpolated)
                        sums[finite, component] += interpolated[finite]
                        counts[finite, component] += 1
                positive = counts > 0
                orbit_means[band.key][orbit_index][positive] = (
                    sums[positive] / counts[positive]
                )
                orbit_coverage[band.key][orbit_index] = counts / 8.0

    representative = {}
    coverage = {}
    orbit_spread = {}
    for band in BANDS:
        values = orbit_means[band.key]
        finite_count = np.sum(np.isfinite(values), axis=0)
        mean = np.full(values.shape[1:], np.nan)
        valid_mean = finite_count > 0
        mean[valid_mean] = (
            np.nansum(values, axis=0)[valid_mean] / finite_count[valid_mean]
        )
        spread = np.full(values.shape[1:], np.nan)
        valid_spread = finite_count > 1
        squared_difference = (values - mean[None, :, :]) ** 2
        spread[valid_spread] = np.sqrt(
            np.nansum(squared_difference, axis=0)[valid_spread]
            / (finite_count[valid_spread] - 1)
        )
        representative[band.key] = mean
        orbit_spread[band.key] = spread
        coverage[band.key] = np.nanmin(orbit_coverage[band.key], axis=(0, 2))
    return {
        "basis": representative,
        "orbit_mean": orbit_means,
        "orbit_spread": orbit_spread,
        "coverage": coverage,
        "footprint_signs": signs,
    }


def plot_representative_basis(path, grids, source):
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    colors = ("#1f77b4", "#2ca02c", "#9467bd")
    for axis, band in zip(axes, BANDS):
        wavelength = grids[band.key]
        basis = source["basis"][band.key]
        coverage = source["coverage"][band.key] >= MIN_ORBIT_COVERAGE
        for component in range(3):
            axis.plot(
                wavelength[coverage],
                basis[coverage, component],
                color=colors[component],
                linewidth=1.2,
                label=f"EOF {component + 1}",
            )
        axis.axhline(0.0, color="0.35", linewidth=0.7)
        axis.set_title(band.title)
        axis.set_ylabel("Dimensionless waveform")
        axis.grid(alpha=0.2)
        axis.margins(x=0)
    axes[-1].set_xlabel("Wavelength (nm)")
    axes[0].legend(frameon=False, ncol=3)
    fig.suptitle(
        "Representative OCO-2 land EOF waveforms\n"
        "l2_oco_eof.h5; four-orbit, eight-footprint dispersion mapping",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_representative_basis(path, grids, source):
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, "w") as dataset:
        dataset.createDimension("orbit", len(L1B_FILES))
        dataset.createDimension("eof", 3)
        for band in BANDS:
            dimension = f"{band.key}_sample"
            dataset.createDimension(dimension, grids[band.key].size)
            wavelength = dataset.createVariable(
                f"{band.key}_wavelength", "f8", (dimension,)
            )
            wavelength[:] = grids[band.key]
            wavelength.units = "nm"
            basis = dataset.createVariable(
                f"{band.key}_eof_waveform", "f8", (dimension, "eof")
            )
            basis[:] = source["basis"][band.key]
            basis.units = "1"
            orbit_mean = dataset.createVariable(
                f"{band.key}_orbit_mean_eof_waveform",
                "f8",
                ("orbit", dimension, "eof"),
            )
            orbit_mean[:] = source["orbit_mean"][band.key]
            orbit_mean.units = basis.units
            spread = dataset.createVariable(
                f"{band.key}_orbit_standard_deviation",
                "f8",
                (dimension, "eof"),
            )
            spread[:] = source["orbit_spread"][band.key]
            spread.units = basis.units
            coverage = dataset.createVariable(
                f"{band.key}_minimum_orbit_coverage_fraction", "f8", (dimension,)
            )
            coverage[:] = source["coverage"][band.key]
        dataset.source_eof_file = str(DEFAULT_OCO_ROOT / EOF_FILE)
        dataset.source_l1b_dispersion_files = " ".join(L1B_FILES)
        dataset.eof_components = "1 2 3"
        dataset.eof_region = "Land"
        dataset.source_count = 32
        dataset.averaging = "eight-footprint mean within each orbit, followed by equal four-orbit mean"
        dataset.sign_alignment = "canonical footprint waveforms aligned to footprint 1"


def mean_retrieval_residual(paths):
    residuals = []
    noises = []
    reference = None
    metadata = None
    for path in paths:
        with Dataset(path) as dataset:
            wavelength = as_float(dataset["wavelength"])
            band_index = np.asarray(dataset["band_index"][:], dtype=np.int8)
            residual = as_float(dataset["final_residual"])
            noise = as_float(dataset["noise_standard_deviation"])
            if reference is None:
                reference = (wavelength, band_index)
                metadata = {
                    "surface": str(dataset.getncattr("surface")),
                    "aerosol_case": str(dataset.getncattr("aerosol_case")),
                    "truth_xco2_ppm": float(dataset.getncattr("truth_xco2_ppm")),
                }
            elif not (
                np.array_equal(band_index, reference[1])
                and np.allclose(wavelength, reference[0], rtol=0, atol=1e-12)
            ):
                raise RuntimeError(f"retrieval grid changed at {path}")
            residuals.append(residual)
            noises.append(noise)
    return {
        "wavelength": reference[0],
        "band_index": reference[1],
        "residual": np.mean(residuals, axis=0),
        "noise": np.mean(noises, axis=0),
        "residual_ensemble": np.asarray(residuals),
        "noise_ensemble": np.asarray(noises),
        "noise_spread_max": float(np.max(np.ptp(noises, axis=0))),
        "metadata": metadata,
        "count": len(paths),
    }


def fit_band(residual, noise, basis, fit_mask):
    valid = (
        fit_mask
        & np.isfinite(residual)
        & np.isfinite(noise)
        & (noise > 0)
        & np.all(np.isfinite(basis), axis=1)
    )
    if np.count_nonzero(valid) < 4:
        raise RuntimeError("too few samples for three-component EOF fit")
    weighted_basis = basis[valid] / noise[valid, None]
    weighted_target = -residual[valid] / noise[valid]
    coefficients, _, rank, singular = np.linalg.lstsq(
        weighted_basis, weighted_target, rcond=None
    )
    if rank != 3:
        raise RuntimeError(f"EOF basis is rank deficient: rank={rank}")
    correction = basis @ coefficients
    postfit = residual + correction
    prefit_rms = float(np.sqrt(np.mean(residual[valid] ** 2)))
    postfit_rms = float(np.sqrt(np.mean(postfit[valid] ** 2)))
    prefit_noise_rms = float(np.sqrt(np.mean((residual[valid] / noise[valid]) ** 2)))
    postfit_noise_rms = float(np.sqrt(np.mean((postfit[valid] / noise[valid]) ** 2)))
    prefit_power = float(np.sum((residual[valid] / noise[valid]) ** 2))
    postfit_power = float(np.sum((postfit[valid] / noise[valid]) ** 2))
    return {
        "valid": valid,
        "coefficients": coefficients,
        "residual_fit_coefficients": -coefficients,
        "correction": correction,
        "explained": -correction,
        "postfit": postfit,
        "prefit_rms": prefit_rms,
        "postfit_rms": postfit_rms,
        "prefit_noise_rms": prefit_noise_rms,
        "postfit_noise_rms": postfit_noise_rms,
        "weighted_fraction_removed": 1.0 - postfit_power / prefit_power,
        "condition_number": float(singular[0] / singular[-1]),
    }


def save_state_result(path, state, perturbations, class_data, fits, source):
    measurement_count = class_data["corrected"]["wavelength"].size
    with Dataset(path, "w") as dataset:
        dataset.createDimension("retrieval_class", 2)
        dataset.createDimension("band", 3)
        dataset.createDimension("eof", 3)
        dataset.createDimension("measurement", measurement_count)
        class_name = dataset.createVariable("retrieval_class_name", str, ("retrieval_class",))
        class_name[:] = np.asarray(CLASSES, dtype=object)
        wavelength = dataset.createVariable("wavelength", "f8", ("measurement",))
        wavelength[:] = class_data["corrected"]["wavelength"]
        wavelength.units = "nm"
        band_index = dataset.createVariable("band_index", "i1", ("measurement",))
        band_index[:] = class_data["corrected"]["band_index"]
        average_residual = dataset.createVariable(
            "average_terminal_residual", "f8", ("retrieval_class", "measurement")
        )
        noise = dataset.createVariable(
            "noise_standard_deviation", "f8", ("retrieval_class", "measurement")
        )
        eof_basis = dataset.createVariable(
            "representative_eof_waveform", "f8", ("measurement", "eof"), fill_value=np.nan
        )
        fit_mask = dataset.createVariable("eof_fit_mask", "i1", ("measurement",))
        coefficient = dataset.createVariable(
            "eof_scale_coefficient", "f8", ("retrieval_class", "band", "eof")
        )
        correction = dataset.createVariable(
            "eof_forward_model_correction", "f8", ("retrieval_class", "measurement")
        )
        postfit = dataset.createVariable(
            "postfit_residual", "f8", ("retrieval_class", "measurement")
        )
        metric_keys = {
            "prefit_rms": "prefit_rms",
            "postfit_rms": "postfit_rms",
            "prefit_noise_normalized_rms": "prefit_noise_rms",
            "postfit_noise_normalized_rms": "postfit_noise_rms",
            "weighted_fraction_removed": "weighted_fraction_removed",
            "condition_number": "condition_number",
        }
        metric_variables = {
            name: dataset.createVariable(name, "f8", ("retrieval_class", "band"))
            for name in metric_keys
        }
        basis_full = np.full((measurement_count, 3), np.nan)
        mask_full = np.zeros(measurement_count, dtype=np.int8)
        for band in BANDS:
            selected = class_data["corrected"]["band_index"] == band.index
            basis_full[selected] = source["basis"][band.key]
            mask_full[selected] = fits["corrected"][band.key]["valid"].astype(np.int8)
        eof_basis[:] = basis_full
        fit_mask[:] = mask_full
        for class_index, retrieval_class in enumerate(CLASSES):
            average_residual[class_index] = class_data[retrieval_class]["residual"]
            noise[class_index] = class_data[retrieval_class]["noise"]
            full_correction = np.zeros(measurement_count)
            full_postfit = np.array(class_data[retrieval_class]["residual"], copy=True)
            for band_index_value, band in enumerate(BANDS):
                selected = class_data[retrieval_class]["band_index"] == band.index
                fit = fits[retrieval_class][band.key]
                coefficient[class_index, band_index_value] = fit["coefficients"]
                full_correction[selected] = fit["correction"]
                full_postfit[selected] = fit["postfit"]
                for name, variable in metric_variables.items():
                    variable[class_index, band_index_value] = fit[metric_keys[name]]
            correction[class_index] = full_correction
            postfit[class_index] = full_postfit
        for variable in (average_residual, noise, correction, postfit):
            variable.units = "mW m-2 sr-1 nm-1"
        eof_basis.units = "1"
        coefficient.units = "mW m-2 sr-1 nm-1 per dimensionless EOF waveform"
        dataset.truth_state_index = state
        dataset.matched_perturbation_indices = " ".join(map(str, perturbations))
        dataset.matched_perturbation_count = len(perturbations)
        dataset.residual_sign = "F(x) - y"
        dataset.fit_definition = "minimize ||(r + E*c)/sigma||_2 independently per band"
        dataset.source_representative_eof = str(
            OUTPUT_ROOT / "representative_oco2_eofs.nc"
        )
        for key, value in class_data["corrected"]["metadata"].items():
            setattr(dataset, key, value)


def fit_residual_ensembles(class_data, source):
    """Fit every perturbation separately so coefficient scatter is retained."""
    result = {retrieval_class: {} for retrieval_class in CLASSES}
    for retrieval_class in CLASSES:
        data = class_data[retrieval_class]
        for band in BANDS:
            selected = data["band_index"] == band.index
            coverage_mask = source["coverage"][band.key] >= MIN_ORBIT_COVERAGE
            result[retrieval_class][band.key] = [
                fit_band(
                    residual[selected],
                    noise[selected],
                    source["basis"][band.key],
                    coverage_mask,
                )
                for residual, noise in zip(
                    data["residual_ensemble"], data["noise_ensemble"]
                )
            ]
    return result


def plot_state_band(path, state, band, perturbations, class_data,
                    fits, ensemble_fits, source):
    """Requested two-column, residual-over-difference EOF comparison."""
    figure, axes = plt.subplots(
        2, 2, figsize=(15, 9), sharex="col",
        gridspec_kw={"height_ratios": (1.35, 1.0)},
    )
    class_order = ("uncorrected", "corrected")
    selected = class_data["corrected"]["band_index"] == band.index
    wavelength = class_data["corrected"]["wavelength"][selected]
    basis = source["basis"][band.key]

    for column, retrieval_class in enumerate(class_order):
        data = class_data[retrieval_class]
        mean_fit = fits[retrieval_class][band.key]
        valid = mean_fit["valid"]
        residual = data["residual"][selected]
        noise = data["noise"][selected]
        individual = ensemble_fits[retrieval_class][band.key]
        coefficient_samples = np.stack(
            [fit["residual_fit_coefficients"] for fit in individual]
        )
        coefficient_mean = coefficient_samples.mean(axis=0)
        coefficient_sd = coefficient_samples.std(axis=0, ddof=1)
        eof_fit = basis @ coefficient_mean
        difference = residual - eof_fit

        top = axes[0, column]
        bottom = axes[1, column]
        top.plot(
            wavelength[valid], residual[valid], color="black", linewidth=1.25,
            label="Mean terminal residual",
        )
        top.plot(
            wavelength[valid], eof_fit[valid], color="#ff7f0e", linewidth=1.6,
            label="Joint three-EOF fit",
        )
        top.axhline(0.0, color="0.4", linewidth=0.7)
        top.set_title(retrieval_class.capitalize(), fontsize=13)
        top.set_ylabel(
            "Residual radiance\n(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)"
        )
        top.grid(alpha=0.20)
        top.margins(x=0)

        bottom.fill_between(
            wavelength[valid], -noise[valid], noise[valid], color="0.82",
            alpha=0.55, linewidth=0, label="$\\pm1\\sigma$ noise",
        )
        bottom.plot(
            wavelength[valid], difference[valid], color="#d62728",
            linewidth=1.15, label="Mean residual $-$ EOF fit",
        )
        bottom.axhline(0.0, color="0.35", linewidth=0.7)
        bottom.set_ylabel(
            "Fit residual\n(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)"
        )
        bottom.set_xlabel("Wavelength (nm)")
        bottom.grid(alpha=0.20)
        bottom.margins(x=0)

        valid_difference = difference[valid]
        radiance_rms = np.sqrt(np.mean(valid_difference ** 2))
        noise_rms = np.sqrt(np.mean((valid_difference / noise[valid]) ** 2))
        coefficient_text = (
            f"EOF coefficient mean $\\pm$ sample SD "
            f"({len(individual)} perturbations)\n"
            f"EOF 1: {coefficient_mean[0]:+.6g} $\\pm$ {coefficient_sd[0]:.3g}    "
            f"EOF 2: {coefficient_mean[1]:+.6g} $\\pm$ {coefficient_sd[1]:.3g}    "
            f"EOF 3: {coefficient_mean[2]:+.6g} $\\pm$ {coefficient_sd[2]:.3g}\n"
            f"post-fit RMS: {radiance_rms:.6g}; "
            f"noise-normalized RMS: {noise_rms:.4g}"
        )
        figure.text(
            0.265 if column == 0 else 0.755,
            0.035,
            coefficient_text,
            ha="center",
            va="bottom",
            fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white",
                  "edgecolor": "0.75"},
        )

    handles_top, labels_top = axes[0, 0].get_legend_handles_labels()
    handles_bottom, labels_bottom = axes[1, 0].get_legend_handles_labels()
    figure.legend(
        handles_top + handles_bottom,
        labels_top + labels_bottom,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.135),
    )
    metadata = class_data["corrected"]["metadata"]
    figure.suptitle(
        f"State {state:03d}, {band.title}: mean terminal residual and "
        "joint OCO-2 EOF fit\n"
        f"{metadata['surface']}, {metadata['aerosol_case']}, "
        f"truth XCO$_2$={metadata['truth_xco2_ppm']:.0f} ppm",
        fontsize=15,
        y=0.985,
    )
    figure.subplots_adjust(
        left=0.075, right=0.985, bottom=0.205, top=0.885,
        hspace=0.12, wspace=0.17,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_state(path, state, perturbations, class_data, fits, source):
    fig, axes = plt.subplots(3, 2, figsize=(16, 12), sharex=False)
    component_colors = ("#1f77b4", "#2ca02c", "#9467bd")
    for row, band in enumerate(BANDS):
        selected = class_data["corrected"]["band_index"] == band.index
        wavelength = class_data["corrected"]["wavelength"][selected]
        basis = source["basis"][band.key]
        for column, retrieval_class in enumerate(CLASSES):
            axis = axes[row, column]
            residual = class_data[retrieval_class]["residual"][selected]
            noise = class_data[retrieval_class]["noise"][selected]
            fit = fits[retrieval_class][band.key]
            valid = fit["valid"]
            axis.fill_between(
                wavelength,
                -noise,
                noise,
                color="0.82",
                alpha=0.55,
                linewidth=0,
                label="$\\pm1\\sigma$ noise",
            )
            axis.plot(wavelength, residual, color="black", linewidth=1.25,
                      label="Mean terminal residual")
            for component in range(3):
                contribution = -basis[:, component] * fit["coefficients"][component]
                axis.plot(
                    wavelength[valid],
                    contribution[valid],
                    color=component_colors[component],
                    linestyle="--",
                    linewidth=0.85,
                    alpha=0.8,
                    label=f"EOF {component + 1} contribution",
                )
            axis.plot(
                wavelength[valid], fit["explained"][valid], color="#ff7f0e",
                linewidth=1.5, label="Joint 3-EOF fit",
            )
            axis.plot(
                wavelength[valid], fit["postfit"][valid], color="#d62728",
                linewidth=1.0, label="Post-fit residual",
            )
            axis.axhline(0.0, color="0.35", linewidth=0.7)
            axis.set_title(
                f"{band.title}, {retrieval_class}\n"
                f"noise-weighted power removed = "
                f"{100 * fit['weighted_fraction_removed']:.1f}%",
                fontsize=11,
            )
            axis.set_ylabel("Residual radiance\n(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)")
            axis.set_xlabel("Wavelength (nm)")
            axis.grid(alpha=0.20)
            axis.margins(x=0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False)
    metadata = class_data["corrected"]["metadata"]
    fig.suptitle(
        f"State {state:03d}: mean residual versus OCO-2 three-EOF fit\n"
        f"{metadata['surface']}, {metadata['aerosol_case']}, "
        f"truth XCO$_2$={metadata['truth_xco2_ppm']:.0f} ppm; "
        f"matched perturbations={len(perturbations)}",
        fontsize=15,
        y=0.985,
    )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.105, top=0.90,
                        hspace=0.42, wspace=0.20)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, state_results):
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# Joint three-component OCO-2 EOF fits to average terminal residuals.\n"
        )
        stream.write("# Residual sign is F(x)-y; coefficients c minimize ||(r+E*c)/sigma||.\n")
        stream.write(
            "# Stored EOF waveforms are dimensionless; coefficients have "
            "mW m-2 sr-1 nm-1 units.\n"
        )
        stream.write(
            "# state class band n_pert n_fit wavelength_min_nm wavelength_max_nm "
            "eof1_coefficient eof2_coefficient eof3_coefficient prefit_rms "
            "postfit_rms prefit_noise_rms postfit_noise_rms "
            "weighted_fraction_removed condition_number\n"
        )
        for state, result in sorted(state_results.items()):
            for retrieval_class in CLASSES:
                for band in BANDS:
                    fit = result["fits"][retrieval_class][band.key]
                    wavelength = result["class_data"][retrieval_class]["wavelength"][
                        result["class_data"][retrieval_class]["band_index"] == band.index
                    ]
                    valid_wavelength = wavelength[fit["valid"]]
                    coefficients = fit["coefficients"]
                    stream.write(
                        f"{state:03d} {retrieval_class} {band.key} "
                        f"{len(result['perturbations'])} {valid_wavelength.size} "
                        f"{valid_wavelength[0]:.6f} {valid_wavelength[-1]:.6f} "
                        f"{coefficients[0]:.12e} {coefficients[1]:.12e} "
                        f"{coefficients[2]:.12e} {fit['prefit_rms']:.12e} "
                        f"{fit['postfit_rms']:.12e} "
                        f"{fit['prefit_noise_rms']:.12e} "
                        f"{fit['postfit_noise_rms']:.12e} "
                        f"{fit['weighted_fraction_removed']:.12e} "
                        f"{fit['condition_number']:.12e}\n"
                    )


def write_coefficient_ensemble_summary(path, state_results):
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# EOF residual-fit coefficients from individual perturbations.\n"
        )
        stream.write(
            "# Coefficients a minimize ||(r-E*a)/sigma||; units are "
            "mW m-2 sr-1 nm-1.\n"
        )
        stream.write(
            "# state class band n_pert eof1_mean eof1_sd eof2_mean eof2_sd "
            "eof3_mean eof3_sd mean_fit_residual_rms "
            "mean_fit_residual_noise_rms\n"
        )
        for state, result in sorted(state_results.items()):
            for retrieval_class in ("uncorrected", "corrected"):
                data = result["class_data"][retrieval_class]
                for band in BANDS:
                    selected = data["band_index"] == band.index
                    residual = data["residual"][selected]
                    noise = data["noise"][selected]
                    basis = result["source"]["basis"][band.key]
                    fits = result["ensemble_fits"][retrieval_class][band.key]
                    samples = np.stack(
                        [fit["residual_fit_coefficients"] for fit in fits]
                    )
                    mean = samples.mean(axis=0)
                    sd = samples.std(axis=0, ddof=1)
                    valid = result["fits"][retrieval_class][band.key]["valid"]
                    difference = residual - basis @ mean
                    rms = np.sqrt(np.mean(difference[valid] ** 2))
                    noise_rms = np.sqrt(
                        np.mean((difference[valid] / noise[valid]) ** 2)
                    )
                    stream.write(
                        f"{state:03d} {retrieval_class} {band.key} {len(fits)} "
                        f"{mean[0]:.12e} {sd[0]:.12e} "
                        f"{mean[1]:.12e} {sd[1]:.12e} "
                        f"{mean[2]:.12e} {sd[2]:.12e} "
                        f"{rms:.12e} {noise_rms:.12e}\n"
                    )


def main():
    oco_root = Path(
        __import__("os").environ.get("OCO_JACOBIAN_ROOT", DEFAULT_OCO_ROOT)
    )
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    retrievals = completed_retrievals()
    states = matched_state_paths(retrievals)
    if not states:
        raise RuntimeError("no states have matched complete corrected/uncorrected retrievals")
    first_state = min(states)
    first_retrieval = states[first_state]["corrected"][0]
    _, _, grids = target_grid(first_retrieval)

    print(f"matched states: {list(states)}", flush=True)
    source = representative_eof_basis(oco_root, grids)
    representative_path = OUTPUT_ROOT / "representative_oco2_eofs.nc"
    save_representative_basis(representative_path, grids, source)
    plot_representative_basis(
        OUTPUT_ROOT / "representative_oco2_eofs.png", grids, source
    )

    state_results = {}
    for state, paths in states.items():
        perturbations = paths["perturbations"]
        print(
            f"fitting state {state:03d}, matched perturbations {perturbations}",
            flush=True,
        )
        class_data = {
            retrieval_class: mean_retrieval_residual(paths[retrieval_class])
            for retrieval_class in CLASSES
        }
        fits = {retrieval_class: {} for retrieval_class in CLASSES}
        for retrieval_class in CLASSES:
            for band in BANDS:
                selected = class_data[retrieval_class]["band_index"] == band.index
                coverage_mask = source["coverage"][band.key] >= MIN_ORBIT_COVERAGE
                fits[retrieval_class][band.key] = fit_band(
                    class_data[retrieval_class]["residual"][selected],
                    class_data[retrieval_class]["noise"][selected],
                    source["basis"][band.key],
                    coverage_mask,
                )
        ensemble_fits = fit_residual_ensembles(class_data, source)
        state_results[state] = {
            "perturbations": perturbations,
            "class_data": class_data,
            "fits": fits,
            "ensemble_fits": ensemble_fits,
            "source": source,
        }
        save_state_result(
            OUTPUT_ROOT / f"state{state:03d}_eof_residual_fit.nc",
            state,
            perturbations,
            class_data,
            fits,
            source,
        )
        plot_state(
            OUTPUT_ROOT / f"state{state:03d}_eof_residual_fit.png",
            state,
            perturbations,
            class_data,
            fits,
            source,
        )
        for band in BANDS:
            plot_state_band(
                OUTPUT_ROOT / f"state{state:03d}_{band.key}_eof_residual_fit.png",
                state,
                band,
                perturbations,
                class_data,
                fits,
                ensemble_fits,
                source,
            )

    summary_path = OUTPUT_ROOT / "eof_fit_summary.dat"
    write_summary(summary_path, state_results)
    coefficient_summary_path = (
        OUTPUT_ROOT / "eof_coefficient_ensemble_summary.dat"
    )
    write_coefficient_ensemble_summary(coefficient_summary_path, state_results)
    print(representative_path)
    print(summary_path)
    print(coefficient_summary_path)


if __name__ == "__main__":
    main()
