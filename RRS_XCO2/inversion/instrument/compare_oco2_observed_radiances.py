#!/usr/bin/env python3
"""Compare synthetic OCO-like radiances with measured OCO-2 L1B radiances.

The L1B radiances are photon spectral radiances per micrometre.  They are
converted sample by sample to energy spectral radiance using ``h*c/lambda``.
Numerically, W m-2 sr-1 um-1 is identical to mW m-2 sr-1 nm-1, which is the
unit stored by ``process_truth_map.jl``.

The comparison deliberately uses the raw synthetic analyzer response.  It
does not apply the division by M11=0.5 used only in diagnostic plots.
"""

import argparse
from collections import defaultdict, namedtuple
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


PLANCK = 6.62607015e-34  # J s, exact SI definition
LIGHT_SPEED = 299792458.0  # m s-1, exact SI definition
SAMPLES = np.arange(1.0, 1017.0)

SCRIPT_DIR = Path(__file__).resolve().parent
RRS_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_SYNTHETIC_ROOT = RRS_ROOT / "truth_map" / "OCO_radiances"
DEFAULT_OUTPUT_ROOT = DEFAULT_SYNTHETIC_ROOT / "validation_against_OCO2"
DEFAULT_L1B_ROOT = Path("/home/sanghavi/data/OCO2_jacobians")

L1B_FILES = (
    ("TallVeg", "oco2_L1bScND_37214a_210630_B11006r_220326081005.h5"),
    ("Himalayas", "oco2_L1bScND_42610a_220706_B11008r_220831171738.h5"),
    ("Sahara", "oco2_L1bScND_41856a_220515_B11008r_220804213809.h5"),
    ("Andes", "oco2_L1bScND_45442a_230116_B11008r_230209231005.h5"),
)


Band = namedtuple("Band", (
    "name", "label", "index", "observed_variable", "synthetic_variable",
    "wavelength_variable"))


BANDS = (
    Band("o2a", "O$_2$ A-band", 0, "radiance_o2",
         "I_OCO_uncorrected_o2a", "o2a_wavelength"),
    Band("weak_co2", "Weak CO$_2$", 1, "radiance_weak_co2",
         "I_OCO_uncorrected_weak_co2", "weak_co2_wavelength"),
    Band("strong_co2", "Strong CO$_2$", 2, "radiance_strong_co2",
         "I_OCO_uncorrected_strong_co2", "strong_co2_wavelength"),
)

METRIC_NAMES = ("spectral_p05", "spectral_median", "spectral_p95",
                "spectral_mean", "spectral_max")
SUMMARY_PERCENTILES = (0.5, 5.0, 50.0, 95.0, 99.5)


def as_float(variable, key):
    """Read a NetCDF/HDF5 slice and replace masks with NaN."""
    values = variable[key]
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    return np.asarray(values, dtype=np.float64)


def dispersion_wavelength_um(coefficients):
    """Evaluate the OCO dispersion polynomial at 1-based detector samples."""
    return np.polynomial.polynomial.polyval(SAMPLES, coefficients)


def spectral_metrics(radiance):
    """Return p05, median, p95, mean, and maximum along the spectral axis."""
    quantiles = np.nanpercentile(radiance, (5.0, 50.0, 95.0), axis=1).T
    return np.column_stack((
        quantiles,
        np.nanmean(radiance, axis=1),
        np.nanmax(radiance, axis=1),
    ))


def contiguous_chunks(indices, chunk_size):
    for start in range(0, len(indices), chunk_size):
        yield indices[start:start + chunk_size]


def synthetic_wavelength_ranges(synthetic_root):
    files = sorted(synthetic_root.glob("OCO2sims_[0-9][0-9][0-9].nc"))
    if len(files) != 64:
        raise RuntimeError(f"expected 64 synthetic files in {synthetic_root}; found {len(files)}")

    ranges = {}
    with Dataset(files[0]) as dataset:
        for band in BANDS:
            wavelength = as_float(dataset.variables[band.wavelength_variable], slice(None))
            ranges[band.name] = (float(wavelength[0]), float(wavelength[-1]))
    return ranges


def common_wavelength_ranges(l1b_root, native_ranges):
    """Find the wavelength intersection shared by every OCO footprint."""
    observed_lower = defaultdict(list)
    observed_upper = defaultdict(list)
    for _, filename in L1B_FILES:
        path = l1b_root / filename
        if not path.is_file():
            raise FileNotFoundError(path)
        with Dataset(path) as dataset:
            coefficient = as_float(
                dataset["/InstrumentHeader/dispersion_coef_samp"], slice(None))
            for band in BANDS:
                for footprint in range(8):
                    wavelength_nm = 1000.0 * dispersion_wavelength_um(
                        coefficient[band.index, footprint, :])
                    observed_lower[band.name].append(float(np.min(wavelength_nm)))
                    observed_upper[band.name].append(float(np.max(wavelength_nm)))

    ranges = {}
    for band in BANDS:
        native_lower, native_upper = native_ranges[band.name]
        lower = max(native_lower, max(observed_lower[band.name]))
        upper = min(native_upper, min(observed_upper[band.name]))
        if lower >= upper:
            raise RuntimeError(f"no common wavelength coverage for {band.name}")
        ranges[band.name] = (lower, upper)
    return ranges


def read_synthetic_scenes(synthetic_root, wavelength_ranges):
    files = sorted(synthetic_root.glob("OCO2sims_[0-9][0-9][0-9].nc"))
    if len(files) != 64:
        raise RuntimeError(f"expected 64 synthetic files in {synthetic_root}; found {len(files)}")

    scenes = []
    for path in files:
        with Dataset(path) as dataset:
            state = int(dataset.state_index)
            if state != int(path.stem.rsplit("_", 1)[1]):
                raise RuntimeError(f"state metadata disagrees with filename: {path}")
            metadata = {
                "state": state,
                "surface": str(dataset.surface),
                "aerosol": str(dataset.aerosol_case),
                "sif": str(dataset.sif_case),
                "xco2_ppm": int(dataset.xco2_ppm),
            }
            for band in BANDS:
                wavelength = as_float(dataset.variables[band.wavelength_variable], slice(None))
                values = as_float(dataset.variables[band.synthetic_variable], slice(None))
                if wavelength.shape != values.shape or not np.all(np.isfinite(values)):
                    raise RuntimeError(f"invalid {band.name} data in {path}")
                if np.any(values <= 0.0):
                    raise RuntimeError(f"non-positive {band.name} radiance in {path}")
                lower, upper = wavelength_ranges[band.name]
                selection = (wavelength >= lower) & (wavelength <= upper)
                if selection.sum() < 100:
                    raise RuntimeError(f"too little common {band.name} coverage in {path}")
                metric = spectral_metrics(values[None, selection])[0]
                scene = metadata.copy()
                scene["band"] = band.name
                scene.update(dict(zip(METRIC_NAMES, metric)))
                scenes.append(scene)
    return scenes


def read_observed_metrics(l1b_root: Path, wavelength_ranges, chunk_size: int):
    """Read valid-land soundings and return per-sounding spectral metrics.

    Two subsets are retained:
      * ``valid_land``: quality flag 0 and land fraction >=90%;
      * ``geometry_match``: valid_land plus SZA=30+/-5 deg and VZA<=10 deg.
    """
    pooled = defaultdict(list)
    geometry_counts = []

    for context, filename in L1B_FILES:
        path = l1b_root / filename
        if not path.is_file():
            raise FileNotFoundError(path)
        with Dataset(path) as dataset:
            coefficient = as_float(
                dataset["/InstrumentHeader/dispersion_coef_samp"], slice(None))
            quality = as_float(
                dataset["/SoundingGeometry/sounding_qual_flag"], slice(None))
            land = as_float(
                dataset["/SoundingGeometry/sounding_land_fraction"], slice(None))
            sza = as_float(
                dataset["/SoundingGeometry/sounding_solar_zenith"], slice(None))
            vza = as_float(
                dataset["/SoundingGeometry/sounding_zenith"], slice(None))

            valid_land = ((quality == 0) & (land >= 90.0) &
                          np.isfinite(sza) & np.isfinite(vza) &
                          (sza >= 0.0) & (sza <= 85.0) &
                          (vza >= 0.0) & (vza <= 90.0))
            geometry_match = (valid_land & (np.abs(sza - 30.0) <= 5.0) &
                              (vza <= 10.0))
            geometry_counts.append((context, filename, int(valid_land.sum()),
                                    int(geometry_match.sum())))

            for band in BANDS:
                observed = dataset[f"/SoundingMeasurements/{band.observed_variable}"]
                units = str(getattr(observed, "Units", getattr(observed, "units", "")))
                if "Ph sec" not in units or "um" not in units:
                    raise RuntimeError(f"unexpected units for {path}:{band.observed_variable}: {units}")
                lower_nm, upper_nm = wavelength_ranges[band.name]

                for footprint in range(8):
                    wavelength_um = dispersion_wavelength_um(
                        coefficient[band.index, footprint, :])
                    spectral_selection = ((wavelength_um * 1000.0 >= lower_nm) &
                                          (wavelength_um * 1000.0 <= upper_nm))
                    if spectral_selection.sum() < 100:
                        raise RuntimeError(
                            f"too little common coverage for {context} {band.name} footprint {footprint + 1}")
                    energy_per_photon = PLANCK * LIGHT_SPEED / (
                        wavelength_um[spectral_selection] * 1.0e-6)
                    frames = np.flatnonzero(valid_land[:, footprint])

                    for frame_chunk in contiguous_chunks(frames, chunk_size):
                        photons = as_float(observed, (frame_chunk, footprint, slice(None)))
                        photons = photons[:, spectral_selection]
                        # photons s-1 m-2 sr-1 um-1 * J photon-1 =
                        # W m-2 sr-1 um-1, numerically identical to mW ... nm-1.
                        energy_radiance = photons * energy_per_photon[None, :]
                        energy_radiance[(energy_radiance <= 0.0) |
                                        (~np.isfinite(energy_radiance))] = np.nan
                        finite_fraction = np.mean(np.isfinite(energy_radiance), axis=1)
                        keep = finite_fraction >= 0.90
                        if not np.any(keep):
                            continue
                        metrics = spectral_metrics(energy_radiance[keep, :])
                        kept_frames = frame_chunk[keep]
                        matched = geometry_match[kept_frames, footprint]

                        pooled[("valid_land", context, band.name)].append(metrics)
                        if np.any(matched):
                            pooled[("geometry_match", context, band.name)].append(
                                metrics[matched, :])

    arrays = {}
    for key, pieces in pooled.items():
        arrays[key] = np.vstack(pieces)
    for subset in ("valid_land", "geometry_match"):
        for band in BANDS:
            pieces = [arrays[(subset, context, band.name)]
                      for context, _ in L1B_FILES
                      if (subset, context, band.name) in arrays]
            arrays[(subset, "ALL", band.name)] = np.vstack(pieces)
    return arrays, geometry_counts


def write_observed_summary(path, observed):
    header = (
        "# subset context band metric n q0p5 q05 q50 q95 q99p5 minimum maximum\n"
        "# Radiance units: mW m-2 sr-1 nm-1. L1B photons were multiplied by hc/lambda.\n"
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for (subset, context, band), values in sorted(observed.items()):
            for column, metric in enumerate(METRIC_NAMES):
                vector = values[:, column]
                quantiles = np.percentile(vector, SUMMARY_PERCENTILES)
                fields = (subset, context, band, metric, len(vector), *quantiles,
                          np.min(vector), np.max(vector))
                handle.write("%-14s %-10s %-10s %-16s %7d " % fields[:5])
                handle.write(" ".join(f"{value:.10g}" for value in fields[5:]))
                handle.write("\n")


def comparison_class(value: float, observed_vector: np.ndarray):
    q005, q50, q995 = np.percentile(observed_vector, (0.5, 50.0, 99.5))
    ratio = value / q50
    range_flag = "low" if value < q005 else "high" if value > q995 else "within"
    factor2_flag = "yes" if ratio < 0.5 or ratio > 2.0 else "no"
    unit_scale_flag = "yes" if ratio < 0.1 or ratio > 10.0 else "no"
    return q005, q50, q995, ratio, range_flag, factor2_flag, unit_scale_flag


def write_synthetic_metrics(path, scenes, observed):
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# state surface aerosol sif xco2_ppm band spectral_p05 spectral_median "
            "spectral_p95 spectral_mean spectral_max obs_median_q0p5 obs_median_q50 "
            "obs_median_q99p5 median_ratio median_range factor2_median "
            "obs_p95_q0p5 obs_p95_q50 obs_p95_q99p5 p95_ratio p95_range "
            "factor2_p95 unit_scale_flag\n"
        )
        for scene in scenes:
            pool = observed[("geometry_match", "ALL", scene["band"])]
            median_check = comparison_class(
                scene["spectral_median"], pool[:, METRIC_NAMES.index("spectral_median")])
            p95_check = comparison_class(
                scene["spectral_p95"], pool[:, METRIC_NAMES.index("spectral_p95")])
            unit_flag = "yes" if (median_check[-1] == "yes" or p95_check[-1] == "yes") else "no"
            handle.write(
                f"{scene['state']:03d} {scene['surface']} {scene['aerosol']} "
                f"{scene['sif']} {scene['xco2_ppm']} {scene['band']} "
            )
            handle.write(" ".join(f"{scene[name]:.10g}" for name in METRIC_NAMES))
            handle.write(
                f" {median_check[0]:.10g} {median_check[1]:.10g} {median_check[2]:.10g}"
                f" {median_check[3]:.10g} {median_check[4]} {median_check[5]}"
                f" {p95_check[0]:.10g} {p95_check[1]:.10g} {p95_check[2]:.10g}"
                f" {p95_check[3]:.10g} {p95_check[4]} {p95_check[5]} {unit_flag}\n"
            )


def write_geometry_counts(path, counts):
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# context l1b_filename valid_land geometry_match\n")
        for context, filename, valid_count, match_count in counts:
            handle.write(f"{context} {filename} {valid_count} {match_count}\n")


def write_wavelength_ranges(path, native_ranges, comparison_ranges):
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# band synthetic_start_nm synthetic_end_nm comparison_start_nm comparison_end_nm\n")
        for band in BANDS:
            native = native_ranges[band.name]
            common = comparison_ranges[band.name]
            handle.write(
                f"{band.name} {native[0]:.12g} {native[1]:.12g} "
                f"{common[0]:.12g} {common[1]:.12g}\n")


def plot_comparison(path, scenes, observed):
    surfaces = ("urban", "rural", "desert", "forest")
    colors = dict(zip(surfaces, ("#6a51a3", "#31a354", "#e6550d", "#006d2c")))
    markers = {"none": "o", "aod760_0p28": "^"}
    figures, axes = plt.subplots(3, 2, figsize=(12.5, 11.0), constrained_layout=True)

    for row, band in enumerate(BANDS):
        band_scenes = [scene for scene in scenes if scene["band"] == band.name]
        pool = observed[("geometry_match", "ALL", band.name)]
        for column, (metric_name, title) in enumerate((
                ("spectral_median", "Spectral median"),
                ("spectral_p95", "95th spectral percentile (continuum proxy)"))):
            axis = axes[row, column]
            metric_column = METRIC_NAMES.index(metric_name)
            observed_values = pool[:, metric_column]
            q005, q05, q50, q95, q995 = np.percentile(
                observed_values, SUMMARY_PERCENTILES)
            axis.axhspan(q005, q995, color="0.88", zorder=0,
                         label="OCO-2 0.5–99.5%" if row == 0 and column == 0 else None)
            axis.axhspan(q05, q95, color="0.72", zorder=0,
                         label="OCO-2 5–95%" if row == 0 and column == 0 else None)
            axis.axhline(q50, color="black", linewidth=1.5, zorder=1,
                         label="OCO-2 median" if row == 0 and column == 0 else None)

            for surface_index, surface in enumerate(surfaces):
                selected = [scene for scene in band_scenes if scene["surface"] == surface]
                for aerosol_name, marker in markers.items():
                    values = [scene[metric_name] for scene in selected
                              if scene["aerosol"] == aerosol_name]
                    x = surface_index + (-0.10 if aerosol_name == "none" else 0.10)
                    axis.scatter(np.full(len(values), x), values, marker=marker,
                                 s=28, facecolor=colors[surface], edgecolor="white",
                                 linewidth=0.35, alpha=0.82, zorder=3,
                                 label=("No aerosol" if aerosol_name == "none" else "AOD760=0.28")
                                 if row == 0 and column == 1 and surface_index == 0 else None)

            axis.set_xticks(range(len(surfaces)))
            axis.set_xticklabels([name.title() for name in surfaces])
            axis.grid(axis="y", color="0.92", linewidth=0.7)
            axis.set_title(f"{band.label}: {title}")
            axis.set_ylabel("mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
            axis.set_ylim(bottom=0.0)

    axes[0, 0].legend(loc="upper left", fontsize=8)
    axes[0, 1].legend(loc="upper left", fontsize=8)
    figures.suptitle(
        "Synthetic analyzer radiances vs OCO-2 valid land soundings\n"
        "OCO-2 subset: SZA 30±5°, VZA≤10°; synthetic values are not divided by M11",
        fontsize=14,
    )
    figures.savefig(path, dpi=180)
    plt.close(figures)


def print_summary(scenes, observed):
    print("Comparison against pooled geometry-matched OCO-2 soundings")
    print("All values are raw analyzer radiances in mW m-2 sr-1 nm-1")
    for band in BANDS:
        pool = observed[("geometry_match", "ALL", band.name)]
        band_scenes = [scene for scene in scenes if scene["band"] == band.name]
        print(f"\n{band.name}  observed soundings={len(pool)}")
        for metric_name in ("spectral_median", "spectral_p95"):
            column = METRIC_NAMES.index(metric_name)
            oq = np.percentile(pool[:, column], SUMMARY_PERCENTILES)
            synthetic_values = np.asarray([scene[metric_name] for scene in band_scenes])
            ratios = synthetic_values / oq[2]
            outside = np.count_nonzero((synthetic_values < oq[0]) |
                                       (synthetic_values > oq[-1]))
            factor2 = np.count_nonzero((ratios < 0.5) | (ratios > 2.0))
            unit_scale = np.count_nonzero((ratios < 0.1) | (ratios > 10.0))
            print(
                f"  {metric_name:16s} OCO q0.5/q50/q99.5="
                f"{oq[0]:.4g}/{oq[2]:.4g}/{oq[-1]:.4g}; "
                f"synthetic={synthetic_values.min():.4g}..{synthetic_values.max():.4g}; "
                f"outside={outside}/64 factor2={factor2}/64 unit-scale={unit_scale}/64"
            )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l1b-root", type=Path, default=DEFAULT_L1B_ROOT)
    parser.add_argument("--synthetic-root", type=Path, default=DEFAULT_SYNTHETIC_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--chunk-size", type=int, default=512)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be positive")
    args.output_root.mkdir(parents=True, exist_ok=True)

    native_ranges = synthetic_wavelength_ranges(args.synthetic_root)
    wavelength_ranges = common_wavelength_ranges(args.l1b_root, native_ranges)
    scenes = read_synthetic_scenes(args.synthetic_root, wavelength_ranges)
    observed, geometry_counts = read_observed_metrics(
        args.l1b_root, wavelength_ranges, args.chunk_size)

    write_observed_summary(args.output_root / "observed_oco2_summary.dat", observed)
    write_synthetic_metrics(args.output_root / "synthetic_scene_metrics.dat", scenes, observed)
    write_geometry_counts(args.output_root / "observed_sounding_counts.dat", geometry_counts)
    write_wavelength_ranges(args.output_root / "comparison_wavelength_ranges.dat",
                            native_ranges, wavelength_ranges)
    plot_comparison(args.output_root / "synthetic_vs_oco2_radiance_ranges.png",
                    scenes, observed)
    print_summary(scenes, observed)
    print(f"\nWrote validation products to {args.output_root}")


if __name__ == "__main__":
    main()
