#!/usr/bin/env python3
"""Visualize corrected and uncorrected retrievals in physical curve space.

The figure set is intentionally organized around physical products rather than
individual state-vector coordinates:

* aerosol AOD is combined with its retrieved median height to reconstruct
  ``d tau_760 / dz`` for each of the three fixed-width altitude-lognormal
  profiles;
* the three retrieved Legendre coefficients in each band are combined to
  reconstruct the Lambertian surface-reflectance spectrum on the canonical
  coefficient-definition grid;
* surface pressure and the twelve retrieved layer VMRs are combined with the
  four fixed upper layers and the canonical materialized p/T/q atmosphere to
  reconstruct dry-air and CO2 molecular columns, layer thicknesses, and the
  full 16-layer CO2 number-concentration profile.  For visualization, a
  shape-preserving cubic curve is fit to log concentration versus geometric
  layer-center altitude; it passes through every positive layer value without
  ordinary cubic-spline overshoot.  Altitude is displayed logarithmically to
  resolve the near-surface layers, and horizontal layer-center bars show one
  sample standard deviation across the paired noisy retrievals;
* SIF coordinates are expanded to a common wavenumber-linear source view,
  converted to spectral radiance per nm, and shown beneath each CO2 profile
  across the O2 A-band.  Legacy states retrieve ``SIF760`` and ``mSIF``;
  round-4/round-5 states instead display the known 759-nm anchor and reconstruct
  the fixed/derived 760-nm coordinate from the saved state schema. Round 5 also
  fixes the spectral slope, so its SIF curves have no retrieval spread.

Perturbations 01--10 are paired by index between corrected and uncorrected
retrievals.  By default, statistics are formed only from pairs for which both
files are complete, converged, and pass the saved fit-quality criterion.  The
``--include-fit-failures`` diagnostic mode retains complete, state-step-
converged retrievals that fail that spectral-fit gate and marks every affected
card explicitly.  This mode is useful for visualizing compensating retrieval
behavior; it must not be interpreted as declaring those fits acceptable.
Curves are derived for every terminal state before pointwise means and 16--84
percentile envelopes are calculated; this preserves AOD-height covariance and
the nonlinear log-coordinate transformation.  Perturbation 11 is never
included in ensemble statistics and can be displayed separately with
``--show-noiseless``.

The default invocation creates four locally zoomed physical-reconstruction
views: aerosol/surface and combined CO2-profile/SIF figures for the
aerosol-free and aerosol-on selected-CO2, SIF-off truth categories.  Pass
``--sif-case angular_integral760_0p5`` to select the corresponding SIF-on
truth categories;
their output names are distinct from the historical SIF-off names.  CO2 is
shown as layer-mean molecular number concentration against altitude.  The
dry-air column and geometric layer grid are recomputed for every retrieved
surface pressure before forming the ensemble.  Paired corrected-minus-
uncorrected plots are available only as an explicit diagnostic option.  Missing
states retain their position in the 2x2 surface layout as explicit
placeholders.  A partial ensemble is plotted but prominently marked
provisional.
"""

import argparse
import math
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import PchipInterpolator

from co2_plot_metadata import co2_case_label
from retrieval_plot_schema import RetrievalPlotSchema


HERE = Path(__file__).resolve().parent
RRS_ROOT = HERE.parent
TRUTH_TABLE = RRS_ROOT / "truth_map" / "true_states.dat"
SCENE_COMPONENTS = RRS_ROOT / "truth_map" / "scene_components.dat"
VERTICAL_PROFILE_TABLE = (
    RRS_ROOT / "truth_map_aerosols" / "aerosol_vertical_profiles.dat"
)
DEFAULT_OUTPUT_DIR = HERE / "physical_ensemble_visualizations"
DEFAULT_ATMOSPHERIC_PROFILE = (
    HERE / "retrieval_setup" / "retrieval_atmosphere_16layer.nc"
)

SURFACE_ORDER = ("urban", "rural", "desert", "forest")
CATEGORY_LABELS = {
    "no_aerosol": "No aerosol",
    "with_aerosol": "Aerosol AOD$_{760}=0.28$",
}
CATEGORY_AEROSOL_CASES = {
    "no_aerosol": "none",
    "with_aerosol": "aod760_0p28",
}
CATEGORY_SELECTIONS = {
    "both": tuple(CATEGORY_AEROSOL_CASES),
    "no-aerosol": ("no_aerosol",),
    "with-aerosol": ("with_aerosol",),
}
SIF_CASE_ON = "angular_integral760_0p5"
SIF_CASE_METADATA = {
    "off": {
        "slug": "nosif",
        "short": "no SIF",
        "truth": "no-SIF truth",
    },
    SIF_CASE_ON: {
        "slug": "sif_angular_integral760_0p5",
        "short": "SIF on ($2\\pi L_\\lambda(760\\,\\mathrm{nm})=0.5$)",
        "truth": (
            "SIF-on truth "
            "($L_\\lambda(760\\,\\mathrm{nm})=0.5/(2\\pi)$ per sr)"
        ),
    },
}

TRUTH_COLOR = "#111111"
UNCORRECTED_COLOR = "#D55E00"  # colorblind-safe vermilion
CORRECTED_COLOR = "#009E73"    # colorblind-safe bluish green
# The CO2 curves often coincide almost exactly, so use the maximally separated
# Okabe--Ito orange/blue pair in addition to redundant line/marker encodings.
CO2_UNCORRECTED_COLOR = "#E69F00"
CO2_CORRECTED_COLOR = "#0072B2"
# CO2 profiles are often almost coincident.  Give each visual layer its own
# opacity so the reference remains visible without desaturating the retrievals.
CO2_TRUTH_LINE_ALPHA = 0.32
CO2_TRUTH_MARKER_ALPHA = 0.45
CO2_MEAN_ALPHA = 0.95
CO2_MEMBER_ALPHA = 0.035
CO2_ENVELOPE_ALPHA = 0.12
CO2_NOISELESS_ALPHA = 0.55
CO2_ERRORBAR_ALPHA = 0.82
CO2_ALTITUDE_MIN_KM = 0.05
CO2_ALTITUDE_MAX_KM = 16.5
CO2_ALTITUDE_TICKS_KM = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)
EFFECT_COLOR = "#7A5195"
PENDING_COLOR = "#777777"

# Explicit typography hierarchy for the dense multi-panel products.  Keep
# these local rather than changing Matplotlib's global rcParams: the figures
# mix full-width profiles, compact spectral panels, and metadata,
# each of which needs a deliberately different scale.
FONT_CARD_TITLE = 12.0
FONT_AXIS_LABEL = 10.0
FONT_AXIS_TICK = 9.0
FONT_COMPACT_TITLE = 9.8
FONT_COMPACT_LABEL = 9.0
FONT_COMPACT_TICK = 8.5
FONT_METADATA = 8.2
FONT_STATUS = 8.2
FONT_CAMPAIGN = 10.2
FONT_SELECTION = 12.0
FONT_HEADER = 10.5
FONT_AEROSOL_LEGEND = 9.8
FONT_CO2_LEGEND = 9.2
FONT_SNAPSHOT = 8.0

WAVENUMBER_CONVERSION = 1.0e7
SIF_REFERENCE_WAVELENGTH_NM = 760.0
SIF_DEFINITION_VERSION = 2
SIF_UPWELLING_SOLID_ANGLE_SR = 2.0 * math.pi
SIF_ANGULAR_INTEGRAL_760 = 0.5
SIF_RADIANCE_760 = (
    SIF_ANGULAR_INTEGRAL_760 / SIF_UPWELLING_SOLID_ANGLE_SR
)
SIF_DEFINITION = (
    "isotropic BOA radiance normalized by 2pi*L_lambda(760 nm)=0.5"
)
SIF_PROVENANCE_ATTRIBUTES = (
    "sif_definition_version",
    "sif_definition",
    "sif_case_on_label",
    "sif_reference_wavelength_nm",
    "sif_upwelling_solid_angle_sr",
    "sif_angular_integral_760_mW_m-2_nm-1",
    "sif_radiance_760_mW_m-2_sr-1_nm-1",
    "sif_cosine_weighted_irradiance_760_mW_m-2_nm-1",
    "sif_SIF760_mW_m-2_sr-1_per_cm-1",
    "sif_mSIF_mW_m-2_sr-1_per_cm-2",
    "sif_template_wavelength_integral_mW_m-2_sr-1",
)


class SIFProvenanceError(ValueError):
    """A truth/retrieval file belongs to an incompatible SIF campaign."""

SPECIES = (
    ("sulfate", "Sulphate", 0.49, "-", 2.6),
    ("organic_carbon", "Organic carbon", 0.40, "--", 2.1),
    ("utls_sulfate", "UTLS sulphate", 0.10, ":", 1.8),
)

# These reproduce RRSXCO2Common.surface_basis_grids(Float32): the canonical
# coefficient grids, before the strong-band convolution-support shoulder.
# Julia's Float32 StepRange does not land exactly on the nominal short edge,
# hence the explicit sample counts.
BANDS = (
    ("o2a", "o2a", "O$_2$ A", 773.0, 2735),
    ("weak_co2", "weak", "Weak CO$_2$", 1622.0, 1281),
    ("strong_co2", "strong", "Strong CO$_2$", 2084.0, 987),
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inversion-root", type=Path, default=HERE,
        help="Directory containing corrected/ and uncorrected/ retrieval files",
    )
    parser.add_argument("--truth-table", type=Path, default=TRUTH_TABLE)
    parser.add_argument(
        "--scene-components", type=Path, default=SCENE_COMPONENTS,
    )
    parser.add_argument(
        "--vertical-profile-table", type=Path,
        default=VERTICAL_PROFILE_TABLE,
    )
    parser.add_argument(
        "--atmospheric-profile", type=Path,
        default=DEFAULT_ATMOSPHERIC_PROFILE,
        help=(
            "Canonical materialized p/T/q atmosphere written by "
            "retrieval_setup/export_retrieval_atmosphere.jl"
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--campaign-label", default="",
        help=(
            "Optional label placed in every figure title, for example "
            "'Round 3 bottom-layer'"
        ),
    )
    parser.add_argument(
        "--xco2-ppm", type=float, default=400.0,
        help=(
            "Truth column-XCO2 category to plot (default: 400 ppm); ignored "
            "when --bottom-co2-ppm is supplied"
        ),
    )
    parser.add_argument(
        "--bottom-co2-ppm", type=float,
        help=(
            "Select a bottom-layer campaign by its injected layer-16 CO2 "
            "VMR (360, 380, 400, 420, or 440 ppm)"
        ),
    )
    parser.add_argument(
        "--sif-case", choices=tuple(SIF_CASE_METADATA), default="off",
        help=(
            "Truth SIF category to plot: off (default) or "
            "angular_integral760_0p5. The latter selects cases with "
            "2*pi*L_lambda(760 nm)=0.5."
        ),
    )
    parser.add_argument(
        "--aerosol-category", choices=tuple(CATEGORY_SELECTIONS),
        default="both",
        help=(
            "Generate both aerosol categories (default), only no-aerosol, "
            "or only with-aerosol"
        ),
    )
    parser.add_argument(
        "--product", choices=("all", "aerosol-surface", "co2-sif"),
        default="all",
        help=(
            "Generate both physical products (default), only the aerosol/"
            "surface figure, or only the CO2/SIF figure"
        ),
    )
    parser.add_argument(
        "--co2-only", action="store_true",
        help="Generate only the two combined CO2-profile/SIF figures",
    )
    parser.add_argument(
        "--show-noiseless", action="store_true",
        help="Overlay paired perturbation-11 retrievals as thin dotted curves",
    )
    parser.add_argument(
        "--include-fit-failures", action="store_true",
        help=(
            "Include complete, state-step-converged retrievals even when "
            "their saved spectral fit-quality flag fails; affected cards "
            "are explicitly marked"
        ),
    )
    parser.add_argument(
        "--hide-individual", action="store_true",
        help="Suppress the faint individual paired perturbation curves",
    )
    parser.add_argument(
        "--scale-mode", choices=("local", "shared"), default="local",
        help=(
            "Use locally zoomed card scales (default) or per-band scales "
            "shared across surfaces"
        ),
    )
    parser.add_argument(
        "--include-paired-correction", action="store_true",
        help="Also generate the paired corrected-minus-uncorrected diagnostics",
    )
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def read_truth_table(path):
    """Return truth-table rows keyed by integer state index."""
    names = None
    rows = {}
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("# index "):
                names = stripped[2:].split()
            elif stripped and not stripped.startswith("#"):
                if names is None:
                    raise RuntimeError("truth header was not found in %s" % path)
                values = stripped.split()
                row = dict(zip(names, values))
                rows[int(row["index"])] = row
    return rows


def _sif_close(value, expected, label, source, atol=1.0e-13):
    """Return a finite float after a tight campaign-identity check."""
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise SIFProvenanceError(
            "%s has non-numeric %s" % (source, label)
        ) from error
    if not np.isfinite(numeric) or not math.isclose(
            numeric, float(expected), rel_tol=1.0e-12, abs_tol=atol):
        raise SIFProvenanceError(
            "%s has %s=%.16g; expected %.16g" %
            (source, label, numeric, float(expected))
        )
    return numeric


def validate_truth_sif(row, source="truth-table row"):
    """Validate and return the SIF truth coordinates for one selected row.

    State indices were retained when the campaign normalization changed, so a
    legacy ``total_0p5`` row is not interchangeable with the version-2 case.
    Validating the physical coordinates here prevents a stale table from
    producing a plausible-looking but incorrectly normalized figure.
    """
    required = ("sif_case", "SIF760", "mSIF")
    missing = [name for name in required if name not in row]
    if missing:
        raise SIFProvenanceError(
            "%s is missing %s" % (source, ", ".join(missing))
        )

    case = str(row["sif_case"])
    angular_key = "sif_angular_integral760"
    if angular_key not in row:
        # Zero is invariant under the normalization change, so legacy no-SIF
        # rows remain usable.  The same exception is intentionally forbidden
        # for SIF-on rows, where ``sif_total`` meant a wavelength integral.
        if case == "off" and "sif_total" in row:
            angular_key = "sif_total"
        else:
            legacy_note = (
                " (legacy sif_total schema detected)"
                if "sif_total" in row else ""
            )
            raise SIFProvenanceError(
                "%s is missing sif_angular_integral760%s" %
                (source, legacy_note)
            )
    try:
        angular_integral = float(row[angular_key])
        native_sif760 = float(row["SIF760"])
        native_msif = float(row["mSIF"])
    except (TypeError, ValueError) as error:
        raise SIFProvenanceError(
            "%s has non-numeric SIF truth coordinates" % source
        ) from error
    values = (angular_integral, native_sif760, native_msif)
    if not all(np.isfinite(value) for value in values):
        raise SIFProvenanceError(
            "%s has non-finite SIF truth coordinates" % source
        )

    if case == "off":
        if not all(math.isclose(value, 0.0, rel_tol=0.0, abs_tol=1.0e-16)
                   for value in values):
            raise SIFProvenanceError(
                "%s labels SIF off but has nonzero SIF coordinates" % source
            )
    elif case == SIF_CASE_ON:
        _sif_close(
            angular_integral, SIF_ANGULAR_INTEGRAL_760,
            "sif_angular_integral760", source,
        )
        expected_native_sif760 = (
            SIF_RADIANCE_760 * SIF_REFERENCE_WAVELENGTH_NM**2 /
            WAVENUMBER_CONVERSION
        )
        _sif_close(
            native_sif760, expected_native_sif760, "SIF760", source,
            atol=5.0e-15,
        )
    else:
        raise SIFProvenanceError(
            "%s uses stale or unsupported SIF case '%s'; expected '%s'" %
            (source, case, SIF_CASE_ON)
        )

    return {
        "case": case,
        "angular_integral760": angular_integral,
        "SIF760": native_sif760,
        "mSIF": native_msif,
    }


def validate_retrieval_sif(dataset, expected_sif, source):
    """Reject a retrieval whose SIF convention differs from its truth row."""
    if "sif_case" not in dataset.ncattrs():
        # The original full-column no-SIF campaign predates explicit SIF
        # provenance metadata.  Its zero source is invariant under the later
        # normalization correction, so those files remain unambiguous when
        # (and only when) they are paired with a no-SIF truth row.  SIF-on
        # files never receive this compatibility exception because their old
        # and corrected normalizations are physically different.
        state_model = (
            str(dataset.getncattr("retrieval_state_model"))
            if "retrieval_state_model" in dataset.ncattrs() else "legacy"
        )
        if expected_sif["case"] == "off" and state_model == "legacy":
            return
        raise SIFProvenanceError("%s is missing sif_case" % source)
    actual_case = str(dataset.getncattr("sif_case"))
    expected_case = expected_sif["case"]
    if actual_case != expected_case:
        raise SIFProvenanceError(
            "%s has sif_case='%s'; truth requires '%s'" %
            (source, actual_case, expected_case)
        )
    # Historical no-SIF results remain physically valid and need no versioned
    # SIF record.  SIF-on results must carry every version-2 attribute.
    if expected_case == "off":
        return

    missing = [
        name for name in SIF_PROVENANCE_ATTRIBUTES
        if name not in dataset.ncattrs()
    ]
    if missing:
        raise SIFProvenanceError(
            "%s is missing corrected-SIF provenance: %s" %
            (source, ", ".join(missing))
        )
    attributes = {
        name: dataset.getncattr(name) for name in SIF_PROVENANCE_ATTRIBUTES
    }
    if int(attributes["sif_definition_version"]) != SIF_DEFINITION_VERSION:
        raise SIFProvenanceError(
            "%s has sif_definition_version=%s; expected %d" %
            (source, attributes["sif_definition_version"],
             SIF_DEFINITION_VERSION)
        )
    if str(attributes["sif_definition"]) != SIF_DEFINITION:
        raise SIFProvenanceError(
            "%s has an unexpected version-2 sif_definition" % source
        )
    if str(attributes["sif_case_on_label"]) != SIF_CASE_ON:
        raise SIFProvenanceError(
            "%s has an inconsistent sif_case_on_label" % source
        )

    reference = _sif_close(
        attributes["sif_reference_wavelength_nm"],
        SIF_REFERENCE_WAVELENGTH_NM,
        "sif_reference_wavelength_nm", source,
    )
    solid_angle = _sif_close(
        attributes["sif_upwelling_solid_angle_sr"],
        SIF_UPWELLING_SOLID_ANGLE_SR,
        "sif_upwelling_solid_angle_sr", source,
    )
    angular_integral = _sif_close(
        attributes["sif_angular_integral_760_mW_m-2_nm-1"],
        expected_sif["angular_integral760"],
        "sif_angular_integral_760_mW_m-2_nm-1", source,
    )
    radiance = _sif_close(
        attributes["sif_radiance_760_mW_m-2_sr-1_nm-1"],
        SIF_RADIANCE_760,
        "sif_radiance_760_mW_m-2_sr-1_nm-1", source,
    )
    cosine_irradiance = _sif_close(
        attributes["sif_cosine_weighted_irradiance_760_mW_m-2_nm-1"],
        math.pi * radiance,
        "sif_cosine_weighted_irradiance_760_mW_m-2_nm-1", source,
    )
    native_sif760 = _sif_close(
        attributes["sif_SIF760_mW_m-2_sr-1_per_cm-1"],
        expected_sif["SIF760"],
        "sif_SIF760_mW_m-2_sr-1_per_cm-1", source,
        atol=5.0e-15,
    )
    _sif_close(
        attributes["sif_mSIF_mW_m-2_sr-1_per_cm-2"],
        expected_sif["mSIF"],
        "sif_mSIF_mW_m-2_sr-1_per_cm-2", source,
        atol=5.0e-16,
    )
    try:
        template_integral = float(
            attributes["sif_template_wavelength_integral_mW_m-2_sr-1"]
        )
    except (TypeError, ValueError) as error:
        raise SIFProvenanceError(
            "%s has a non-numeric SIF template integral" % source
        ) from error
    if not np.isfinite(template_integral):
        raise SIFProvenanceError(
            "%s has a non-finite SIF template integral" % source
        )

    # Redundant identities make unit/factor-of-pi regressions fail loudly.
    _sif_close(
        solid_angle * radiance, angular_integral,
        "2pi-times-radiance identity", source,
    )
    _sif_close(
        native_sif760 * WAVENUMBER_CONVERSION / reference**2,
        radiance, "wavenumber-to-wavelength SIF760 identity", source,
    )
    _sif_close(
        cosine_irradiance, math.pi * radiance,
        "cosine-weighted irradiance identity", source,
    )


def read_component_sections(path):
    """Read the simple bracketed tables in scene_components.dat."""
    sections = {}
    section = None
    names = None
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("[") and stripped.endswith("]"):
                section = stripped[1:-1]
                sections[section] = []
                names = None
            elif section is not None and stripped.startswith("# "):
                candidate = stripped[2:].split()
                if candidate and candidate[0] in ("case", "species"):
                    names = candidate
            elif section is not None and not stripped.startswith("#") and names:
                values = stripped.split()
                if len(values) >= len(names):
                    sections[section].append(dict(zip(names, values)))
    return sections


def read_vertical_truth(path):
    """Read exact truth medians, widths, and the 16-layer boundaries."""
    names = None
    species_truth = {}
    tops = []
    bottoms = []
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("# nlayer "):
                names = stripped[2:].split()
            elif stripped and not stripped.startswith("#"):
                if names is None:
                    raise RuntimeError(
                        "vertical-profile header was not found in %s" % path
                    )
                row = dict(zip(names, stripped.split()))
                if int(row["nlayer"]) != 16:
                    continue
                species = row["species"]
                species_truth[species] = {
                    "aod760": float(row["tau760"]),
                    "z0": float(row["z_median_km"]),
                    "sigma": float(row["sigma_ln_z"]),
                }
                if species == "sulfate":
                    tops.append(float(row["z_top_km"]))
                    bottoms.append(float(row["z_bottom_km"]))
    if not tops:
        raise RuntimeError("no 16-layer records were found in %s" % path)
    boundaries = sorted(set(tops + bottoms), reverse=True)
    return species_truth, np.asarray(boundaries, dtype=float)


AVOGADRO = 6.02214179e23
GAS_CONSTANT = 8.3144598
DRY_AIR_MOLAR_MASS = 28.9644e-3
WATER_MOLAR_MASS = 18.01534e-3
GRAVITY = 9.8032465


def compute_model_atmosphere(template, psurf):
    """Recompute hydrostatic columns and geometry for one surface pressure.

    This is a unit-for-unit translation of
    ``CoreRT.compute_atmos_profile_fields``.  Pressure is in hPa, layer
    thickness in m, and dry-air VCD in molecules cm-2.
    """
    pressure = template["pressure_interface"].copy()
    if not np.isfinite(psurf) or psurf <= pressure[-2]:
        raise ValueError(
            "surface pressure %.6g hPa is not below the bottom-layer top %.6g hPa"
            % (psurf, pressure[-2])
        )
    pressure[-1] = psurf
    temperature = template["temperature"]
    specific_humidity = template["specific_humidity"]
    water_dry_ratio = (
        specific_humidity / (1.0 - specific_humidity) *
        (DRY_AIR_MOLAR_MASS / WATER_MOLAR_MASS)
    )
    dry_fraction = 1.0 / (1.0 + water_dry_ratio)
    water_fraction = water_dry_ratio * dry_fraction
    mixture_molar_mass = (
        dry_fraction * DRY_AIR_MOLAR_MASS +
        water_fraction * WATER_MOLAR_MASS
    )
    pressure_thickness = np.diff(pressure)
    dry_air_column = (
        dry_fraction * AVOGADRO * pressure_thickness /
        (mixture_molar_mass * GRAVITY * 100.0)
    )
    layer_thickness = (
        np.log(pressure[1:] / pressure[:-1]) *
        GAS_CONSTANT * temperature /
        (GRAVITY * mixture_molar_mass)
    )
    altitude_interface = np.zeros(pressure.size, dtype=float)
    for layer in range(layer_thickness.size - 1, -1, -1):
        altitude_interface[layer] = (
            altitude_interface[layer + 1] +
            layer_thickness[layer] / 1000.0
        )
    return {
        "pressure_interface": pressure,
        "dry_air_column": dry_air_column,
        "layer_thickness": layer_thickness,
        "altitude_interface": altitude_interface,
    }


def read_model_atmosphere(path):
    """Load and regression-check the canonical materialized p/T/q template."""
    if not path.is_file():
        raise RuntimeError(
            "canonical retrieval atmosphere is missing: %s\n"
            "Generate it with retrieval_setup/export_retrieval_atmosphere.jl"
            % path
        )
    required = (
        "pressure_interface", "temperature", "specific_humidity",
        "dry_air_vertical_column", "layer_thickness", "altitude_interface",
    )
    with Dataset(str(path)) as dataset:
        missing = [name for name in required if name not in dataset.variables]
        if missing:
            raise RuntimeError(
                "canonical atmosphere %s is missing %s" %
                (path, ", ".join(missing))
            )
        template = {
            name: np.asarray(dataset[name][:], dtype=float)
            for name in required
        }
        template["reference_psurf"] = float(
            dataset.getncattr("reference_surface_pressure_hpa")
        )
    if (
            template["pressure_interface"].shape != (17,) or
            template["temperature"].shape != (16,) or
            template["specific_humidity"].shape != (16,)):
        raise RuntimeError("canonical atmosphere must contain exactly 16 layers")
    reference = compute_model_atmosphere(
        template, template["reference_psurf"]
    )
    checks = (
        ("dry-air column", reference["dry_air_column"],
         template["dry_air_vertical_column"]),
        ("layer thickness", reference["layer_thickness"],
         template["layer_thickness"]),
        ("altitude interfaces", reference["altitude_interface"],
         template["altitude_interface"]),
    )
    for label, calculated, saved in checks:
        # The artifact records Float32 model arithmetic, whereas this plotting
        # reconstruction intentionally evaluates the same equations in
        # Float64.  The resulting geometry differs by only a few millimetres.
        if not np.allclose(calculated, saved, rtol=1.0e-5, atol=5.0e-6):
            raise RuntimeError(
                "Python hydrostatic reconstruction failed the saved %s check" %
                label
            )
    return template


def canonical_band_grid(long_nm, count):
    """Reproduce a canonical Float32 0.1-cm^-1 surface-basis grid."""
    first_nu = np.float32(1.0e7 / long_nm)
    index = np.arange(count, dtype=np.float32)
    nu = first_nu + np.float32(0.1) * index
    x = 2.0 * (
        (nu.astype(float) - float(nu[0])) /
        (float(nu[-1]) - float(nu[0]))
    ) - 1.0
    wavelength = 1.0e7 / nu.astype(float)
    order = np.argsort(wavelength)
    basis = np.vstack((
        np.ones_like(x),
        x,
        0.5 * (3.0 * x * x - 1.0),
    ))
    return wavelength[order], basis[:, order]


BAND_GRIDS = {
    state_band: canonical_band_grid(long_nm, count)
    for state_band, _, _, long_nm, count in BANDS
}


def native_to_physical(schema, state, fixed_upper_co2_ppm,
                       saved_xco2_ppm=None):
    """Transform one saved terminal state into plotting parameters."""
    native = schema.native_dict(np.asarray(state, dtype=float))
    co2_ppm = np.full(16, float(fixed_upper_co2_ppm), dtype=float)
    for layer in range(5, 17):
        name = "co2_vmr_layer%02d" % layer
        if name not in native:
            raise ValueError("terminal state is missing %s" % name)
        co2_ppm[layer - 1] = 1.0e6 * native[name]
    record = {
        "aod": {},
        "z0": {},
        "surface": {},
        "sif": {
            "SIF759": native["SIF759"],
            "SIF760": native["SIF760"],
            "mSIF": native["mSIF"],
            "known_wavelength_nm": (
                schema.known_wavelength_nm if schema.has_fixed_sif else None
            ),
            "known_Lnu": schema.known_lnu if schema.has_fixed_sif else None,
            "sif760_status": schema.sif760_status,
            "msif_status": schema.msif_status,
            "state_model": schema.state_model,
        },
        "psurf": native["psurf"],
        "co2_ppm": co2_ppm,
    }
    if saved_xco2_ppm is not None:
        record["saved_xco2_ppm"] = float(saved_xco2_ppm)
    for species, _, _, _, _ in SPECIES:
        record["aod"][species] = math.exp(
            native["ln_%s_aod760" % species]
        )
        record["z0"][species] = math.exp(native["ln_%s_z0" % species])
    for state_band, _, _, _, _ in BANDS:
        record["surface"][state_band] = np.asarray([
            native["%s_surface_P%d" % (state_band, order)]
            for order in range(3)
        ])
    return record


def read_retrieval(path, expected_state, expected_perturbation, expected_class,
                   expected_sif, include_fit_failures=False):
    """Return a validated physical terminal state, or (None, reason)."""
    if not path.is_file():
        return None, "missing"
    try:
        with Dataset(str(path)) as dataset:
            if int(dataset.getncattr("retrieval_complete")) != 1:
                return None, "incomplete"
            if int(dataset.getncattr("truth_state_index")) != expected_state:
                return None, "wrong state"
            if int(dataset.getncattr("perturbation_index")) != expected_perturbation:
                return None, "wrong perturbation"
            if str(dataset.getncattr("measurement_class")) != expected_class:
                return None, "wrong class"
            validate_retrieval_sif(dataset, expected_sif, str(path))
            try:
                schema = RetrievalPlotSchema.from_dataset(dataset, str(path))
                schema.validate_truth_row(
                    {
                        "sif_case": expected_sif["case"],
                        "sif_angular_integral760": (
                            expected_sif["angular_integral760"]
                        ),
                        "SIF760": expected_sif["SIF760"],
                        "mSIF": expected_sif["mSIF"],
                    },
                    "validated truth state for %s" % path,
                )
            except ValueError as error:
                raise SIFProvenanceError(
                    "retrieval-state schema mismatch: %s" % error
                ) from error
            if "converged" in dataset.ncattrs():
                if int(dataset.getncattr("converged")) != 1:
                    return None, "not converged"
            fit_quality_ok = True
            if "fit_quality_ok" in dataset.ncattrs():
                fit_quality_ok = int(dataset.getncattr("fit_quality_ok")) == 1
                if not fit_quality_ok and not include_fit_failures:
                    return None, "fit quality failed"
            state = np.asarray(dataset["final_state"][:], dtype=float)
            if state.ndim != 1 or len(schema.active_names) != state.size:
                return None, "state schema mismatch"
            if not np.all(np.isfinite(state)):
                return None, "non-finite state"
            if "fixed_upper_co2_ppm" in dataset.ncattrs():
                fixed_upper_co2_ppm = float(
                    dataset.getncattr("fixed_upper_co2_ppm")
                )
            elif "truth_xco2_ppm" in dataset.ncattrs():
                fixed_upper_co2_ppm = float(
                    dataset.getncattr("truth_xco2_ppm")
                )
            else:
                return None, "missing fixed-upper-CO2 metadata"
            saved_xco2_ppm = (
                float(np.asarray(dataset["XCO2"][:]))
                if "XCO2" in dataset.variables else None
            )
            record = native_to_physical(
                schema, state, fixed_upper_co2_ppm, saved_xco2_ppm
            )
            # Plot-only provenance.  Leading underscores distinguish these
            # diagnostics from physical state coordinates consumed below.
            record["_plot_schema"] = schema
            record["_source"] = str(path)
            record["_fit_quality_ok"] = fit_quality_ok
            if "final_band_reduced_chi_squared" in dataset.variables:
                record["_final_band_reduced_chi_squared"] = np.asarray(
                    dataset["final_band_reduced_chi_squared"][:], dtype=float
                )
            reason = "ok" if fit_quality_ok else "fit quality failed (included)"
            return record, reason
    except SIFProvenanceError as error:
        return None, "SIF provenance mismatch: %s" % error
    except Exception as error:
        return None, "unreadable: %s" % error


def retrieval_path(root, retrieval_class, state_index, perturbation_index):
    return root / retrieval_class / (
        "retrieval_state%03d_perturbation%02d.nc" %
        (state_index, perturbation_index)
    )


def load_paired_ensemble(root, state_index, expected_sif,
                         include_fit_failures=False):
    """Load matched corrected/uncorrected noisy draws and optional p11."""
    paired = {"corrected": [], "uncorrected": [], "indices": []}
    valid_counts = {"corrected": 0, "uncorrected": 0}
    fit_failure_indices = {"corrected": [], "uncorrected": []}
    reasons = {"corrected": {}, "uncorrected": {}}
    for perturbation in range(1, 11):
        records = {}
        for retrieval_class in ("corrected", "uncorrected"):
            path = retrieval_path(
                root, retrieval_class, state_index, perturbation
            )
            record, reason = read_retrieval(
                path, state_index, perturbation, retrieval_class,
                expected_sif,
                include_fit_failures=include_fit_failures,
            )
            if reason.startswith("SIF provenance mismatch: "):
                raise SIFProvenanceError(reason.split(": ", 1)[1])
            records[retrieval_class] = record
            reasons[retrieval_class][perturbation] = reason
            if record is not None:
                valid_counts[retrieval_class] += 1
                if not record.get("_fit_quality_ok", True):
                    fit_failure_indices[retrieval_class].append(perturbation)
        if records["corrected"] is not None and records["uncorrected"] is not None:
            paired["indices"].append(perturbation)
            paired["corrected"].append(records["corrected"])
            paired["uncorrected"].append(records["uncorrected"])

    noiseless = {}
    for retrieval_class in ("corrected", "uncorrected"):
        path = retrieval_path(root, retrieval_class, state_index, 11)
        record, reason = read_retrieval(
            path, state_index, 11, retrieval_class,
            expected_sif,
            include_fit_failures=include_fit_failures,
        )
        if reason.startswith("SIF provenance mismatch: "):
            raise SIFProvenanceError(reason.split(": ", 1)[1])
        noiseless[retrieval_class] = record
        reasons[retrieval_class][11] = reason
        if record is not None and not record.get("_fit_quality_ok", True):
            fit_failure_indices[retrieval_class].append(11)

    paired["valid_counts"] = valid_counts
    paired["missing_indices"] = [
        index for index in range(1, 11) if index not in paired["indices"]
    ]
    paired["noiseless"] = noiseless
    paired["noiseless_paired"] = (
        noiseless["corrected"] is not None and
        noiseless["uncorrected"] is not None
    )
    paired["reasons"] = reasons
    paired["fit_failure_indices"] = fit_failure_indices
    schemas = [
        record["_plot_schema"]
        for retrieval_class in ("corrected", "uncorrected")
        for record in (
            paired[retrieval_class] +
            ([noiseless[retrieval_class]]
             if noiseless[retrieval_class] is not None else [])
        )
    ]
    if schemas:
        reference = schemas[0]
        mismatched = [
            schema.source for schema in schemas
            if schema.identity != reference.identity
        ]
        if mismatched:
            raise SIFProvenanceError(
                "state %03d mixes incompatible retrieval-state schemas; "
                "first mismatch: %s" % (state_index, mismatched[0])
            )
        paired["plot_schema"] = reference
    else:
        paired["plot_schema"] = None
    return paired


def truth_record(row, aerosol_cases, vertical_truth, validated_sif=None,
                 plot_schema=None):
    """Construct the physical truth record for one truth-table row."""
    aerosol_case = row["aerosol_case"]
    case = aerosol_cases[aerosol_case]
    truth_xco2_ppm = float(row["xco2_ppm"])
    if "bottom_co2_ppm" in row:
        background_ppm = float(row["background_co2_ppm"])
        bottom_layer = int(row["bottom_layer_index"])
        if not 1 <= bottom_layer <= 16:
            raise RuntimeError("invalid bottom_layer_index in truth table")
        co2_ppm = np.full(16, background_ppm, dtype=float)
        co2_ppm[bottom_layer - 1] = float(row["bottom_co2_ppm"])
    else:
        co2_ppm = np.full(16, truth_xco2_ppm, dtype=float)
    sif_truth = (
        validate_truth_sif(row) if validated_sif is None else validated_sif
    )
    if plot_schema is not None:
        try:
            plot_schema.validate_truth_row(row)
        except ValueError as error:
            raise SIFProvenanceError(str(error)) from error
    if plot_schema is not None and plot_schema.has_fixed_sif:
        try:
            truth_sif = plot_schema.truth_sif_coordinates(row)
        except ValueError as error:
            raise SIFProvenanceError(str(error)) from error
        truth_msif = truth_sif["mSIF"]
        truth_sif759 = truth_sif["SIF759"]
        truth_sif760 = truth_sif["SIF760"]
        known_wavelength_nm = plot_schema.known_wavelength_nm
        known_lnu = plot_schema.known_lnu
        sif760_status = plot_schema.sif760_status
        msif_status = plot_schema.msif_status
        state_model = plot_schema.state_model
    else:
        truth_msif = sif_truth["mSIF"]
        truth_sif760 = sif_truth["SIF760"]
        delta_nu = (
            WAVENUMBER_CONVERSION / SIF_REFERENCE_WAVELENGTH_NM -
            WAVENUMBER_CONVERSION / 759.0
        )
        truth_sif759 = truth_sif760 - delta_nu * truth_msif
        known_wavelength_nm = None
        known_lnu = None
        sif760_status = "retrieved"
        msif_status = "retrieved"
        state_model = "legacy"
    record = {
        "aod": {},
        "z0": {},
        "surface": {},
        "sif": {
            "angular_integral760": sif_truth["angular_integral760"],
            "SIF759": truth_sif759,
            "SIF760": truth_sif760,
            "mSIF": truth_msif,
            "known_wavelength_nm": known_wavelength_nm,
            "known_Lnu": known_lnu,
            "sif760_status": sif760_status,
            "msif_status": msif_status,
            "state_model": state_model,
        },
        "psurf": float(row["psurf_hpa"]),
        "co2_ppm": co2_ppm,
        "saved_xco2_ppm": truth_xco2_ppm,
    }
    aod_columns = {
        "sulfate": "sulfate_AOD760",
        "organic_carbon": "organic_AOD760",
        "utls_sulfate": "utls_sulfate_AOD760",
    }
    for species, _, _, _, _ in SPECIES:
        record["aod"][species] = float(case[aod_columns[species]])
        record["z0"][species] = vertical_truth[species]["z0"]
    for state_band, truth_band, _, _, _ in BANDS:
        record["surface"][state_band] = np.asarray([
            float(row["%s_P%d" % (truth_band, order)])
            for order in range(3)
        ])
    return record


def lognormal_column_fraction(model_top, median, sigma):
    argument = (math.log(model_top) - math.log(median)) / (
        sigma * math.sqrt(2.0)
    )
    return 0.5 * (1.0 + math.erf(argument))


def aerosol_profile(record, species, sigma, altitude, model_top):
    """Return finite-column-normalized d tau_760 / dz in km^-1."""
    aod = record["aod"][species]
    median = record["z0"][species]
    if aod == 0.0:
        return np.zeros_like(altitude)
    if aod < 0.0 or median <= 0.0 or sigma <= 0.0:
        raise ValueError("invalid aerosol profile parameters")
    z = np.maximum(altitude, np.finfo(float).tiny)
    exponent = -0.5 * (np.log(z / median) / sigma) ** 2
    density = np.exp(exponent) / (
        z * sigma * math.sqrt(2.0 * math.pi)
    )
    density[altitude <= 0.0] = 0.0
    normalization = lognormal_column_fraction(model_top, median, sigma)
    return aod * density / normalization


def surface_curve(record, state_band):
    _, basis = BAND_GRIDS[state_band]
    return np.dot(record["surface"][state_band], basis)


def sif_curve(record):
    """Return the two-parameter O2 A-band SIF state as radiance per nm.

    The retrieval state is affine in wavenumber spectral density.  The final
    factor is the absolute spectral-coordinate Jacobian |d nu / d lambda|;
    no factor of pi belongs here because this is radiance, not the internal
    hemispheric-irradiance representation used by the surface kernel.
    """
    wavelength, _ = BAND_GRIDS["o2a"]
    wavenumber = WAVENUMBER_CONVERSION / wavelength
    reference_wavenumber = (
        WAVENUMBER_CONVERSION / SIF_REFERENCE_WAVELENGTH_NM
    )
    radiance_per_cm1 = (
        record["sif"]["SIF760"] +
        record["sif"]["mSIF"] * (wavenumber - reference_wavenumber)
    )
    return radiance_per_cm1 * WAVENUMBER_CONVERSION / wavelength**2


def sif_reference_radiance(record):
    """Return exact ``L_lambda(760 nm)`` from the native SIF760 state."""
    return (
        float(record["sif"]["SIF760"]) * WAVENUMBER_CONVERSION /
        SIF_REFERENCE_WAVELENGTH_NM**2
    )


def sif_reference_angular_integral(record):
    """Return the experiment's unweighted isotropic upward angular integral."""
    return SIF_UPWELLING_SOLID_ANGLE_SR * sif_reference_radiance(record)


def sif_known_anchor_radiance(record):
    """Return ``(wavelength_nm, L_lambda)`` for a fixed SIF anchor, if any."""
    wavelength = record["sif"].get("known_wavelength_nm")
    native_lnu = record["sif"].get("known_Lnu")
    if wavelength is None or native_lnu is None:
        return None
    return (
        float(wavelength),
        float(native_lnu) * WAVENUMBER_CONVERSION / float(wavelength)**2,
    )


def sif_stack(records):
    return np.asarray([sif_curve(record) for record in records])


def curve_statistics(curves):
    array = np.asarray(curves, dtype=float)
    return (
        np.mean(array, axis=0),
        np.percentile(array, 16.0, axis=0),
        np.percentile(array, 84.0, axis=0),
    )


def status_text(ensemble):
    count = len(ensemble["indices"])
    noisy_fit_failures = {
        retrieval_class: [
            index for index in ensemble.get("fit_failure_indices", {}).get(
                retrieval_class, []
            ) if 1 <= index <= 10
        ]
        for retrieval_class in ("corrected", "uncorrected")
    }
    failure_count = sum(map(len, noisy_fit_failures.values()))
    noiseless_fit_failures = [
        retrieval_class
        for retrieval_class in ("corrected", "uncorrected")
        if 11 in ensemble.get("fit_failure_indices", {}).get(
            retrieval_class, []
        )
    ]
    failure_count += len(noiseless_fit_failures)

    def failure_label():
        label = "noisy corr %d, unc %d" % (
            len(noisy_fit_failures["corrected"]),
            len(noisy_fit_failures["uncorrected"]),
        )
        if noiseless_fit_failures:
            short = {"corrected": "corr", "uncorrected": "unc"}
            label += "; p11 " + "/".join(
                short[value] for value in noiseless_fit_failures
            )
        return label

    if count == 10:
        text, color = "10/10 complete", "#E7F5EC"
        if failure_count:
            text += "\nfit failures included\n" + failure_label()
            color = "#FCE8D5"
        return text, color
    if count >= 2:
        missing = ", ".join("%02d" % value for value in ensemble["missing_indices"])
        text = "%d/10 provisional\nmissing %s" % (count, missing)
        if failure_count:
            text += "\nfit failures included: " + failure_label()
        return text, "#FFF3CD"
    return "%d/10 pending" % count, "#EEEEEE"


def annotate_status(axis, ensemble):
    """Place the small completeness badge in a low-conflict plot corner."""
    text, color = status_text(ensemble)
    axis.text(
        0.985, 0.965, text, transform=axis.transAxes,
        ha="right", va="top", fontsize=FONT_STATUS,
        bbox=dict(boxstyle="round,pad=0.24", facecolor=color,
                  edgecolor="0.72", alpha=0.94),
        zorder=20,
    )


def prepare_categories(args):
    if args.bottom_co2_ppm is None:
        selection_field = "xco2_ppm"
        selection_value = args.xco2_ppm
        selection_label = "XCO$_2$"
        selection_slug = "xco2"
    else:
        selection_field = "bottom_co2_ppm"
        selection_value = args.bottom_co2_ppm
        selection_label = "bottom-layer CO$_2$"
        selection_slug = "bottom_co2"
    if not np.isfinite(selection_value) or selection_value <= 0.0:
        raise ValueError("the selected CO2 value must be finite and positive")
    sif_metadata = SIF_CASE_METADATA[args.sif_case]
    truth_rows = read_truth_table(args.truth_table)
    sections = read_component_sections(args.scene_components)
    aerosol_cases = {
        row["case"]: row for row in sections.get("AEROSOL_CASES", [])
    }
    vertical_truth, _archived_boundaries = read_vertical_truth(
        args.vertical_profile_table
    )
    atmosphere = read_model_atmosphere(args.atmospheric_profile)
    # The historical aerosol table predates the current shared profile
    # materialization.  Species medians/widths remain authoritative there, but
    # all layer lines and finite-column normalization must use current geometry.
    boundaries = atmosphere["altitude_interface"].copy()
    category_names = CATEGORY_SELECTIONS[args.aerosol_category]
    categories = {}
    for category in category_names:
        aerosol_case = CATEGORY_AEROSOL_CASES[category]
        categories[category] = {}
        for surface in SURFACE_ORDER:
            candidates = [
                state_index for state_index, candidate in truth_rows.items()
                if candidate["surface"] == surface and
                candidate["sif_case"] == args.sif_case and
                candidate["aerosol_case"] == aerosol_case and
                selection_field in candidate and math.isclose(
                    float(candidate[selection_field]), selection_value,
                    rel_tol=0.0, abs_tol=1.0e-9
                )
            ]
            if len(candidates) != 1:
                if not candidates and args.sif_case == SIF_CASE_ON:
                    legacy_candidates = [
                        state_index
                        for state_index, candidate in truth_rows.items()
                        if candidate["surface"] == surface and
                        candidate["sif_case"] == "total_0p5" and
                        candidate["aerosol_case"] == aerosol_case and
                        selection_field in candidate and math.isclose(
                            float(candidate[selection_field]), selection_value,
                            rel_tol=0.0, abs_tol=1.0e-9,
                        )
                    ]
                    if legacy_candidates:
                        raise SIFProvenanceError(
                            "%s contains only legacy total_0p5 SIF state(s) "
                            "%s for %s/%s/%s=%g; regenerate the version-2 "
                            "truth table before plotting" % (
                                args.truth_table, legacy_candidates, surface,
                                aerosol_case, selection_field, selection_value,
                            )
                        )
                raise RuntimeError(
                    "expected exactly one %s/%s/SIF=%s/%s=%g state; found %s" %
                    (surface, aerosol_case, args.sif_case, selection_field,
                     selection_value, candidates)
                )
            state_index = candidates[0]
            row = truth_rows[state_index]
            validated_sif = validate_truth_sif(
                row, "%s state %03d" % (args.truth_table, state_index)
            )
            ensemble = load_paired_ensemble(
                args.inversion_root, state_index, validated_sif,
                include_fit_failures=args.include_fit_failures,
            )
            plot_schema = ensemble.get("plot_schema")
            categories[category][surface] = {
                "state_index": state_index,
                "row": row,
                "truth": truth_record(
                    row, aerosol_cases, vertical_truth, validated_sif,
                    plot_schema=plot_schema,
                ),
                "ensemble": ensemble,
            }
    selected_rows = [
        card["row"]
        for category_cards in categories.values()
        for card in category_cards.values()
    ]
    selected_xco2 = sorted(set(
        round(float(row["xco2_ppm"]), 10) for row in selected_rows
    ))
    if len(selected_xco2) != 1:
        raise RuntimeError(
            "selected physical-ensemble states do not share one XCO2"
        )
    schemas = [
        card["ensemble"].get("plot_schema")
        for category_cards in categories.values()
        for card in category_cards.values()
        if card["ensemble"].get("plot_schema") is not None
    ]
    if schemas:
        reference_schema = schemas[0]
        if any(schema.identity != reference_schema.identity
               for schema in schemas[1:]):
            raise SIFProvenanceError(
                "selected cards mix incompatible retrieval-state schemas"
            )
        state_description = (
            reference_schema.description() if reference_schema.has_fixed_sif
            else None
        )
    else:
        state_description = None
    selection = {
        "value": float(selection_value),
        "label": selection_label,
        "slug": selection_slug,
        "sif_case": args.sif_case,
        "sif_slug": sif_metadata["slug"],
        "sif_short": sif_metadata["short"],
        "sif_truth": sif_metadata["truth"],
        "state_description": state_description,
        "campaign_label": str(args.campaign_label).strip(),
        "category_names": category_names,
        "xco2_ppm": selected_xco2[0],
        "display": co2_case_label(
            selection_value if args.bottom_co2_ppm is not None else None,
            selected_xco2[0],
            truth=True,
        ),
    }
    return categories, boundaries, atmosphere, selection


def profile_stack(records, species, sigma, altitude, model_top):
    return np.asarray([
        aerosol_profile(record, species, sigma, altitude, model_top)
        for record in records
    ])


def surface_stack(records, state_band):
    return np.asarray([
        surface_curve(record, state_band) for record in records
    ])


def co2_molecular_profile(record, atmosphere):
    """Return model-consistent molecular columns and mean concentrations.

    CO2 VMR is relative to dry air.  Therefore the layer CO2 column is VMR
    times the dry-air VCD, not the total moist-air column.  Dividing by the
    geometric layer thickness gives the layer-mean molecular concentration.
    """
    model = compute_model_atmosphere(atmosphere, record["psurf"])
    co2_vmr = np.asarray(record["co2_ppm"], dtype=float) * 1.0e-6
    co2_column = co2_vmr * model["dry_air_column"]
    number_concentration = (
        co2_column / (100.0 * model["layer_thickness"])
    )
    if (
            np.any(co2_column <= 0.0) or
            np.any(number_concentration <= 0.0) or
            not np.all(np.isfinite(number_concentration))):
        raise ValueError("invalid CO2 molecular profile")
    reconstructed_column = np.sum(
        number_concentration * 100.0 * model["layer_thickness"]
    )
    if not np.isclose(
            reconstructed_column, np.sum(co2_column), rtol=2.0e-13):
        raise RuntimeError("CO2 concentration failed layer-column conservation")
    xco2_ppm = (
        1.0e6 * np.sum(co2_column) / np.sum(model["dry_air_column"])
    )
    if "saved_xco2_ppm" in record and not np.isclose(
            xco2_ppm, record["saved_xco2_ppm"], rtol=3.0e-6,
            atol=2.0e-4):
        raise RuntimeError(
            "reconstructed XCO2 %.8f ppm disagrees with saved %.8f ppm" %
            (xco2_ppm, record["saved_xco2_ppm"])
        )
    return {
        **model,
        "co2_column": co2_column,
        "number_concentration": number_concentration,
        "altitude_center": 0.5 * (
            model["altitude_interface"][:-1] +
            model["altitude_interface"][1:]
        ),
        "xco2_ppm": xco2_ppm,
        "total_co2_column": float(np.sum(co2_column)),
    }


def interpolate_co2_profile(profile, altitude):
    """PCHIP interpolation of log concentration through layer centers.

    PCHIP is shape preserving and the logarithmic coordinate guarantees a
    positive molecular concentration.  Below the lowest layer center, use the
    endpoint tangent rather than extrapolating the first cubic polynomial.
    """
    centers = profile["altitude_center"][::-1]
    values = profile["number_concentration"][::-1]
    interpolator = PchipInterpolator(centers, np.log(values), extrapolate=False)
    clipped = np.clip(altitude, centers[0], centers[-1])
    log_curve = np.asarray(interpolator(clipped), dtype=float)
    below = altitude < centers[0]
    above = altitude > centers[-1]
    derivative = interpolator.derivative()
    if np.any(below):
        log_curve[below] = (
            math.log(values[0]) + float(derivative(centers[0])) *
            (altitude[below] - centers[0])
        )
    if np.any(above):
        log_curve[above] = (
            math.log(values[-1]) + float(derivative(centers[-1])) *
            (altitude[above] - centers[-1])
        )
    reconstructed = np.exp(interpolator(centers))
    if not np.allclose(reconstructed, values, rtol=2.0e-13):
        raise RuntimeError("CO2 interpolation does not pass its layer values")
    return np.exp(log_curve)


def co2_concentration_stack(records, atmosphere, altitude):
    return np.asarray([
        interpolate_co2_profile(
            co2_molecular_profile(record, atmosphere), altitude
        )
        for record in records
    ])


def co2_profile_limits(cards, atmosphere, altitude, include_noiseless=False):
    """Positive log-scale limits from the concentration curves displayed."""
    values = []
    for card in cards.values():
        ensemble = card["ensemble"]
        if len(ensemble["indices"]) < 2:
            continue
        truth = co2_molecular_profile(card["truth"], atmosphere)
        values.append(interpolate_co2_profile(truth, altitude))
        for retrieval_class in ("corrected", "uncorrected"):
            values.extend(co2_concentration_stack(
                ensemble[retrieval_class], atmosphere, altitude
            ))
            noiseless = ensemble["noiseless"][retrieval_class]
            if include_noiseless and noiseless is not None:
                values.append(interpolate_co2_profile(
                    co2_molecular_profile(noiseless, atmosphere), altitude
                ))
    if not values:
        return (5.0e14, 2.0e16)
    lower = min(float(np.min(value)) for value in values)
    upper = max(float(np.max(value)) for value in values)
    return (lower / 1.10, upper * 1.10)


def sif_limits(cards, include_noiseless=False):
    """Linear SIF-radiance limits spanning every curve that will be drawn."""
    values = [np.asarray([0.0])]
    for card in cards.values():
        ensemble = card["ensemble"]
        if len(ensemble["indices"]) < 2:
            continue
        values.append(sif_curve(card["truth"]))
        for retrieval_class in ("corrected", "uncorrected"):
            values.extend(sif_stack(ensemble[retrieval_class]))
            noiseless = ensemble["noiseless"][retrieval_class]
            if include_noiseless and noiseless is not None:
                values.append(sif_curve(noiseless))
    lower = min(float(np.min(value)) for value in values)
    upper = max(float(np.max(value)) for value in values)
    span = max(upper - lower, 1.0e-3)
    padding = 0.10 * span
    return lower - padding, upper + padding


def main_profile_limits(cards, altitude, model_top, include_noiseless=False):
    """Shared category limit based on displayed envelopes, not remote outliers."""
    maxima = [0.0]
    for card in cards.values():
        if len(card["ensemble"]["indices"]) < 2:
            continue
        for species, _, sigma, _, _ in SPECIES:
            maxima.append(np.max(aerosol_profile(
                card["truth"], species, sigma, altitude, model_top
            )))
            for retrieval_class in ("corrected", "uncorrected"):
                curves = profile_stack(
                    card["ensemble"][retrieval_class], species, sigma,
                    altitude, model_top,
                )
                mean, _, upper = curve_statistics(curves)
                maxima.extend([np.max(mean), np.max(upper)])
                noiseless = card["ensemble"]["noiseless"][retrieval_class]
                if include_noiseless and noiseless is not None:
                    maxima.append(np.max(aerosol_profile(
                        noiseless, species, sigma, altitude, model_top
                    )))
    upper = max(1.0e-5, 1.08 * max(maxima))
    return (0.0, upper)


def effect_profile_limits(cards, altitude, model_top, include_noiseless=False):
    maxima = [0.0]
    for card in cards.values():
        ensemble = card["ensemble"]
        if len(ensemble["indices"]) < 2:
            continue
        for species, _, sigma, _, _ in SPECIES:
            corrected = profile_stack(
                ensemble["corrected"], species, sigma, altitude, model_top
            )
            uncorrected = profile_stack(
                ensemble["uncorrected"], species, sigma, altitude, model_top
            )
            differences = corrected - uncorrected
            mean, lower, upper = curve_statistics(differences)
            maxima.extend([
                np.max(np.abs(mean)), np.max(np.abs(lower)),
                np.max(np.abs(upper)),
            ])
            if include_noiseless and ensemble["noiseless_paired"]:
                noiseless = (
                    aerosol_profile(
                        ensemble["noiseless"]["corrected"], species, sigma,
                        altitude, model_top,
                    ) -
                    aerosol_profile(
                        ensemble["noiseless"]["uncorrected"], species, sigma,
                        altitude, model_top,
                    )
                )
                maxima.append(np.max(np.abs(noiseless)))
    limit = max(1.0e-5, 1.10 * max(maxima))
    return (-limit, limit)


def main_surface_limits(cards, state_band, include_noiseless=False):
    """Shared absolute-reflectance range for one band across all surfaces."""
    values = []
    for card in cards.values():
        ensemble = card["ensemble"]
        if len(ensemble["indices"]) < 2:
            continue
        values.append(surface_curve(card["truth"], state_band))
        for retrieval_class in ("corrected", "uncorrected"):
            curves = surface_stack(ensemble[retrieval_class], state_band)
            mean, lower, upper = curve_statistics(curves)
            values.extend([mean, lower, upper])
            noiseless = ensemble["noiseless"][retrieval_class]
            if include_noiseless and noiseless is not None:
                values.append(surface_curve(noiseless, state_band))
    if not values:
        return (0.0, 1.0)
    lower = min(np.min(value) for value in values)
    upper = max(np.max(value) for value in values)
    span = max(upper - lower, 5.0e-4)
    padding = 0.12 * span
    return (lower - padding, upper + padding)


def effect_surface_limits(cards, state_band, include_noiseless=False):
    maximum = 0.0
    for card in cards.values():
        ensemble = card["ensemble"]
        if len(ensemble["indices"]) < 2:
            continue
        corrected = surface_stack(ensemble["corrected"], state_band)
        uncorrected = surface_stack(ensemble["uncorrected"], state_band)
        differences = corrected - uncorrected
        _, lower, upper = curve_statistics(differences)
        maximum = max(
            maximum, np.max(np.abs(lower)), np.max(np.abs(upper))
        )
        if include_noiseless and ensemble["noiseless_paired"]:
            noiseless = (
                surface_curve(
                    ensemble["noiseless"]["corrected"], state_band
                ) -
                surface_curve(
                    ensemble["noiseless"]["uncorrected"], state_band
                )
            )
            maximum = max(maximum, np.max(np.abs(noiseless)))
    limit = max(1.0e-6, 1.12 * maximum)
    return (-limit, limit)


def draw_layer_boundaries(axis, boundaries, maximum_altitude):
    for boundary in boundaries:
        if 0.0 < boundary < maximum_altitude:
            axis.axhline(
                boundary, color="0.60", linewidth=0.45, alpha=0.25, zorder=0
            )


def placeholder_message(state_index, ensemble=None):
    if ensemble is not None:
        sif_mismatches = [
            reason
            for retrieval_class in ("corrected", "uncorrected")
            for reason in ensemble["reasons"][retrieval_class].values()
            if reason.startswith("SIF provenance mismatch:")
        ]
        if sif_mismatches:
            return (
                "Stale/incompatible SIF retrievals excluded\nstate %03d\n"
                "version-2 truth and retrieval provenance must match" %
                state_index
            )
        failed_fit = sorted({
            index
            for retrieval_class in ("corrected", "uncorrected")
            for index, reason in ensemble["reasons"][retrieval_class].items()
            if 1 <= index <= 10 and reason == "fit quality failed"
        })
        if failed_fit:
            if failed_fit == list(range(failed_fit[0], failed_fit[-1] + 1)):
                label = (
                    "p%02d--%02d" % (failed_fit[0], failed_fit[-1])
                    if len(failed_fit) > 1 else "p%02d" % failed_fit[0]
                )
            else:
                label = ", ".join("p%02d" % value for value in failed_fit)
            return (
                "Excluded from ensemble\nstate %03d\n%s fail the saved "
                "fit-quality criterion" % (state_index, label)
            )
        nonmissing = sum(
            reason != "missing"
            for retrieval_class in ("corrected", "uncorrected")
            for index, reason in ensemble["reasons"][retrieval_class].items()
            if 1 <= index <= 10
        )
        if nonmissing:
            return (
                "Insufficient valid paired retrievals\nstate %03d\n"
                "fewer than two pairs pass all filters" % state_index
            )
    return (
        "Pending\nstate %03d\nno paired noisy retrievals yet" % state_index
    )


def draw_placeholder(profile_axis, surface_axes, state_index, ensemble=None):
    all_axes = [profile_axis] + list(surface_axes)
    for axis in all_axes:
        axis.set_facecolor("#F2F2F2")
        axis.patch.set_hatch("//")
        axis.patch.set_edgecolor("#E0E0E0")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_color("#CCCCCC")
    for axis in surface_axes:
        axis.set_axis_off()
    profile_axis.text(
        0.5, 0.49,
        placeholder_message(state_index, ensemble),
        transform=profile_axis.transAxes, ha="center", va="center",
        fontsize=FONT_CARD_TITLE, color=PENDING_COLOR, fontweight="semibold",
        bbox=dict(boxstyle="round,pad=0.55", facecolor="white",
                  edgecolor="#BBBBBB", alpha=0.93),
    )


def draw_main_card(profile_axis, surface_axes, card, altitude, boundaries,
                   profile_limits, surface_limits_by_band, show_individual,
                   show_noiseless):
    ensemble = card["ensemble"]
    count = len(ensemble["indices"])
    state_index = card["state_index"]
    if count < 2:
        draw_placeholder(profile_axis, surface_axes, state_index, ensemble)
        profile_axis.set_title(
            "%s — state %03d" % (card["row"]["surface"].title(), state_index),
            fontsize=FONT_CARD_TITLE, pad=8,
        )
        return

    model_top = float(boundaries[0])
    if profile_limits is None:
        profile_limits = main_profile_limits(
            {"local": card}, altitude, model_top,
            include_noiseless=show_noiseless,
        )
    draw_layer_boundaries(profile_axis, boundaries, altitude[-1])
    truth_has_aerosol = sum(card["truth"]["aod"].values()) > 0.0
    if truth_has_aerosol:
        for species, _, sigma, linestyle, linewidth in SPECIES:
            truth = aerosol_profile(
                card["truth"], species, sigma, altitude, model_top
            )
            profile_axis.plot(
                truth, altitude, color=TRUTH_COLOR, linestyle=linestyle,
                linewidth=linewidth + 0.9, zorder=5,
            )
    else:
        profile_axis.axvline(
            0.0, color=TRUTH_COLOR, linewidth=3.1, zorder=5
        )

    for retrieval_class, color, width_scale, mean_zorder in (
            ("uncorrected", UNCORRECTED_COLOR, 0.86, 6),
            ("corrected", CORRECTED_COLOR, 0.52, 7)):
        for species, _, sigma, linestyle, linewidth in SPECIES:
            curves = profile_stack(
                ensemble[retrieval_class], species, sigma, altitude, model_top
            )
            if show_individual:
                for curve in curves:
                    profile_axis.plot(
                        curve, altitude, color=color, linestyle=linestyle,
                        linewidth=0.65, alpha=0.075, zorder=2,
                    )
            mean, lower, upper = curve_statistics(curves)
            if count >= 3:
                profile_axis.fill_betweenx(
                    altitude, lower, upper, color=color, alpha=0.14,
                    linewidth=0.0, zorder=3,
                )
            profile_axis.plot(
                mean, altitude, color=color, linestyle=linestyle,
                linewidth=max(1.05, width_scale * linewidth),
                zorder=mean_zorder,
            )
            if show_noiseless and ensemble["noiseless_paired"]:
                noiseless = aerosol_profile(
                    ensemble["noiseless"][retrieval_class], species, sigma,
                    altitude, model_top,
                )
                profile_axis.plot(
                    noiseless, altitude, color=color, linestyle=linestyle,
                    linewidth=0.85, alpha=0.72, marker="s", markevery=75,
                    markersize=2.0, markerfacecolor="white", zorder=5,
                )

    profile_axis.set_xlim(profile_limits)
    profile_axis.set_ylim(0.0, altitude[-1])
    profile_axis.set_xlabel(
        "$d\\tau_{760}/dz$ (km$^{-1}$)", fontsize=FONT_AXIS_LABEL
    )
    profile_axis.set_ylabel("Altitude (km)", fontsize=FONT_AXIS_LABEL)
    profile_axis.grid(axis="x", alpha=0.18)
    profile_axis.tick_params(labelsize=FONT_AXIS_TICK)
    annotate_status(profile_axis, ensemble)
    profile_axis.set_title(
        "%s — state %03d" % (card["row"]["surface"].title(), state_index),
        fontsize=FONT_CARD_TITLE, pad=8,
    )

    for band_axis, (state_band, _, band_label, _, _) in zip(
            surface_axes, BANDS):
        if surface_limits_by_band is None:
            surface_limits = main_surface_limits(
                {"local": card}, state_band,
                include_noiseless=show_noiseless,
            )
        else:
            surface_limits = surface_limits_by_band[state_band]
        wavelength, _ = BAND_GRIDS[state_band]
        truth = surface_curve(card["truth"], state_band)
        band_axis.plot(
            wavelength, truth, color=TRUTH_COLOR, linewidth=3.2, zorder=5
        )
        for retrieval_class, color, mean_width, mean_zorder in (
                ("uncorrected", UNCORRECTED_COLOR, 2.05, 6),
                ("corrected", CORRECTED_COLOR, 1.05, 7)):
            curves = surface_stack(ensemble[retrieval_class], state_band)
            if show_individual:
                for curve in curves:
                    band_axis.plot(
                        wavelength, curve, color=color, linewidth=0.55,
                        alpha=0.075, zorder=2,
                    )
            mean, lower, upper = curve_statistics(curves)
            if count >= 3:
                band_axis.fill_between(
                    wavelength, lower, upper, color=color, alpha=0.23,
                    linewidth=0.0, zorder=3,
                )
            band_axis.plot(
                wavelength, mean, color=color, linewidth=mean_width,
                zorder=mean_zorder,
            )
            if show_noiseless and ensemble["noiseless_paired"]:
                noiseless = surface_curve(
                    ensemble["noiseless"][retrieval_class], state_band
                )
                band_axis.plot(
                    wavelength, noiseless, color=color, linestyle=":",
                    linewidth=0.95, alpha=0.82, zorder=5,
                )
        band_axis.set_title(
            band_label, fontsize=FONT_COMPACT_TITLE, pad=3,
        )
        band_axis.set_ylim(surface_limits)
        band_axis.set_xlim(wavelength[0], wavelength[-1])
        band_axis.set_xlabel(
            "Wavelength (nm)", fontsize=FONT_COMPACT_LABEL,
        )
        band_axis.tick_params(labelsize=FONT_COMPACT_TICK)
        band_axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        band_axis.ticklabel_format(
            axis="y", style="plain", useOffset=False
        )
        band_axis.grid(alpha=0.16)
    surface_axes[0].set_ylabel("Reflectance", fontsize=FONT_COMPACT_LABEL)


def draw_co2_placeholder(axis, state_index, surface, ensemble=None):
    axis.set_facecolor("#F2F2F2")
    axis.patch.set_hatch("//")
    axis.patch.set_edgecolor("#E0E0E0")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_color("#CCCCCC")
    axis.text(
        0.5, 0.49,
        placeholder_message(state_index, ensemble),
        transform=axis.transAxes, ha="center", va="center",
        fontsize=FONT_CARD_TITLE, color=PENDING_COLOR, fontweight="semibold",
        bbox=dict(boxstyle="round,pad=0.55", facecolor="white",
                  edgecolor="#BBBBBB", alpha=0.93),
    )
    axis.set_title(
        "%s — state %03d" % (surface.title(), state_index),
        fontsize=FONT_CARD_TITLE, pad=8,
    )


def draw_sif_placeholder(axis, state_index=None, ensemble=None):
    axis.set_facecolor("#F2F2F2")
    axis.patch.set_hatch("//")
    axis.patch.set_edgecolor("#E0E0E0")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_color("#CCCCCC")
    message = (
        placeholder_message(state_index, ensemble)
        if state_index is not None else
        "SIF state unavailable with current filters"
    )
    axis.text(
        0.5, 0.5, message,
        transform=axis.transAxes, ha="center", va="center",
        fontsize=FONT_METADATA, color=PENDING_COLOR,
    )


def draw_sif_card(axis, annotation_axis, card, ylimits, show_individual,
                  show_noiseless):
    """Draw the complete two-parameter SIF state beneath one CO2 profile."""
    # SIF summaries are deliberately kept off the data axes.  Some retrieval
    # ensembles span zero, so there is no consistently empty in-panel corner.
    annotation_axis.set_axis_off()
    ensemble = card["ensemble"]
    plot_schema = ensemble.get("plot_schema")
    count = len(ensemble["indices"])
    if count < 2:
        draw_sif_placeholder(axis, card["state_index"], ensemble)
        return

    if ylimits is None:
        ylimits = sif_limits(
            {"local": card}, include_noiseless=show_noiseless
        )
    wavelength, _ = BAND_GRIDS["o2a"]
    reference_index = int(np.argmin(
        np.abs(wavelength - SIF_REFERENCE_WAVELENGTH_NM)
    ))
    axis.axhline(0.0, color="0.45", linewidth=0.7, alpha=0.55, zorder=0)
    truth_curve = sif_curve(card["truth"])
    truth_at_reference = sif_reference_radiance(card["truth"])
    known_anchor = sif_known_anchor_radiance(card["truth"])
    if known_anchor is not None:
        axis.axvline(
            known_anchor[0], color="0.30", linestyle="--",
            linewidth=0.75, alpha=0.58, zorder=1,
        )
    axis.plot(
        wavelength, truth_curve, color=TRUTH_COLOR,
        linewidth=2.4, alpha=0.72, zorder=5,
    )
    axis.plot(
        SIF_REFERENCE_WAVELENGTH_NM, truth_at_reference,
        color=TRUTH_COLOR, marker="D", markersize=4.0,
        markeredgecolor="white", markeredgewidth=0.55,
        linestyle="none", zorder=9,
    )
    if known_anchor is not None:
        axis.plot(
            known_anchor[0], known_anchor[1], color=TRUTH_COLOR,
            marker="*", markersize=6.5, markeredgecolor="white",
            markeredgewidth=0.55, linestyle="none", zorder=10,
        )

    reference_summary = {}
    for retrieval_class, color, mean_width, mean_zorder in (
            ("uncorrected", UNCORRECTED_COLOR, 2.05, 6),
            ("corrected", CORRECTED_COLOR, 1.35, 7)):
        curves = sif_stack(ensemble[retrieval_class])
        if show_individual:
            for curve in curves:
                axis.plot(
                    wavelength, curve, color=color, linewidth=0.55,
                    alpha=0.09, zorder=2,
                )
        mean, lower, upper = curve_statistics(curves)
        reference_values = np.asarray([
            sif_reference_radiance(record)
            for record in ensemble[retrieval_class]
        ])
        reference_summary[retrieval_class] = (
            float(np.mean(reference_values)),
            float(np.std(reference_values, ddof=1)),
        )
        if count >= 3:
            axis.fill_between(
                wavelength, lower, upper, color=color, alpha=0.18,
                linewidth=0.0, zorder=3,
            )
        axis.plot(
            wavelength, mean, color=color, linewidth=mean_width,
            zorder=mean_zorder,
        )
        if show_noiseless and ensemble["noiseless_paired"]:
            noiseless = sif_curve(ensemble["noiseless"][retrieval_class])
            axis.plot(
                wavelength, noiseless, color=color, linestyle=":",
                linewidth=1.0, alpha=0.82, marker="s",
                markevery=[reference_index], markersize=3.2,
                markerfacecolor="white", markeredgewidth=0.8, zorder=8,
            )

    axis.set_xlim(wavelength[0], wavelength[-1])
    axis.set_ylim(ylimits)
    axis.set_xlabel("Wavelength (nm)", fontsize=FONT_COMPACT_LABEL)
    axis.set_ylabel(
        "Surface SIF $L_\\lambda$\n"
        "(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)",
        fontsize=FONT_COMPACT_LABEL,
    )
    title = "O$_2$ A-band BOA SIF-radiance state"
    if plot_schema is not None and plot_schema.has_fixed_sif:
        if plot_schema.sif_mode == "off":
            title += " — $L_{759}$ and slope fixed to zero"
        elif plot_schema.is_fixed_sif_model:
            title += " — fixed $L_{759}$ and slope"
        else:
            title += " — known $L_{759}$; retrieved slope"
    axis.set_title(title, fontsize=FONT_COMPACT_TITLE, pad=3)
    if plot_schema is not None and plot_schema.has_fixed_sif:
        if plot_schema.is_fixed_sif_model and plot_schema.sif_mode == "on":
            annotation = (
                "Fixed: $L_{759}$ = %.5f; $m_{SIF}$ = %.5g (native units)\n"
                "Derived $L_{760}$ = %.5f; no SIF retrieval spread "
                "(paired p01--p10; n=%d)" % (
                    known_anchor[1], plot_schema.fixed_msif,
                    truth_at_reference, count,
                )
            )
        elif plot_schema.sif_mode == "on":
            annotation = (
                "Known exactly: $L_{759}$ = %.5f; $m_{SIF}$ retrieved, "
                "$L_{760}$ derived\n"
                "Derived $L_{760}$: state-space truth %.5f; "
                "unc %.5f $\\pm$ %.5f; "
                "corr %.5f $\\pm$ %.5f  (paired p01--p10; n=%d)" % (
                    known_anchor[1], truth_at_reference,
                    reference_summary["uncorrected"][0],
                    reference_summary["uncorrected"][1],
                    reference_summary["corrected"][0],
                    reference_summary["corrected"][1], count,
                )
            )
        else:
            annotation = (
                "Known exactly: $L_{759}=0$; $m_{SIF}=0$ and "
                "$L_{760}=0$ are fixed (not retrieved)\n"
                "paired p01--p10; n=%d" % count
            )
    else:
        annotation = (
            "Truth: $L_{760}$ = %.5f; $2\\pi L_{760}$ = %.3f\n"
            "unc %.5f $\\pm$ %.5f; corr %.5f $\\pm$ %.5f"
            "  (paired p01--p10; n=%d)" % (
                truth_at_reference,
                sif_reference_angular_integral(card["truth"]),
                reference_summary["uncorrected"][0],
                reference_summary["uncorrected"][1],
                reference_summary["corrected"][0],
                reference_summary["corrected"][1],
                count,
            )
        )
    annotation_axis.text(
        0.5, -0.12, annotation,
        transform=annotation_axis.transAxes, ha="center", va="top",
        fontsize=FONT_METADATA, color="0.22", linespacing=1.16,
        clip_on=False,
    )
    axis.tick_params(labelsize=FONT_COMPACT_TICK)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
    axis.ticklabel_format(axis="y", style="plain", useOffset=False)
    axis.grid(alpha=0.16)

def draw_co2_card(axis, card, atmosphere, altitude, xlimits,
                  show_individual, show_noiseless):
    """Draw model-consistent CO2 molecular concentration against altitude."""
    ensemble = card["ensemble"]
    count = len(ensemble["indices"])
    state_index = card["state_index"]
    surface = card["row"]["surface"]
    if count < 2:
        draw_co2_placeholder(axis, state_index, surface, ensemble)
        return

    if xlimits is None:
        xlimits = co2_profile_limits(
            {"local": card}, atmosphere, altitude,
            include_noiseless=show_noiseless
        )
    truth_record_value = card["truth"]
    truth = co2_molecular_profile(truth_record_value, atmosphere)
    for edge in truth["altitude_interface"][1:-1]:
        if 0.0 < edge < altitude[-1]:
            axis.axhline(
                edge, color="0.58", linewidth=0.48, alpha=0.22, zorder=0
            )
    fixed_boundary = truth["altitude_interface"][4]
    axis.axhline(
        fixed_boundary, color="0.35", linestyle="--",
        linewidth=0.9, alpha=0.65, zorder=1,
    )
    truth_curve = interpolate_co2_profile(truth, altitude)
    # These three profiles can be nearly coincident.  Keep truth thin, then
    # distinguish the retrieval classes by line geometry as well as colour;
    # otherwise the superposed black/vermilion/green strokes appear grey.
    axis.plot(
        truth_curve, altitude, color=TRUTH_COLOR, linewidth=0.90,
        alpha=CO2_TRUTH_LINE_ALPHA, zorder=4,
    )
    visible_truth_points = truth["altitude_center"] <= altitude[-1]
    axis.plot(
        truth["number_concentration"][visible_truth_points],
        truth["altitude_center"][visible_truth_points],
        linestyle="none", marker="o", markersize=2.4,
        markerfacecolor="none", markeredgecolor=TRUTH_COLOR,
        markeredgewidth=0.7, alpha=CO2_TRUTH_MARKER_ALPHA, zorder=4,
    )

    psurf_summary = []
    column_summary = []
    xco2_summary = {}
    mean_marker_heights = truth["altitude_center"][visible_truth_points]
    for retrieval_class, color, mean_width, mean_zorder, linestyle, marker, marker_offset in (
            ("uncorrected", CO2_UNCORRECTED_COLOR, 2.20, 7,
             "-", "s", 0),
            ("corrected", CO2_CORRECTED_COLOR, 2.20, 8,
             (0, (3.0, 1.8)), "^", 1)):
        records = ensemble[retrieval_class]
        profiles = [
            co2_molecular_profile(record, atmosphere) for record in records
        ]
        curves = np.asarray([
            interpolate_co2_profile(profile, altitude) for profile in profiles
        ])
        if show_individual:
            for profile in profiles:
                axis.plot(
                    interpolate_co2_profile(profile, altitude), altitude,
                    color=color, linewidth=0.45,
                    alpha=CO2_MEMBER_ALPHA, zorder=2,
                )
        mean, lower, upper = curve_statistics(curves)
        curve_sigma = np.std(curves, axis=0, ddof=1)
        if count >= 3:
            axis.fill_betweenx(
                altitude, lower, upper, color=color,
                alpha=CO2_ENVELOPE_ALPHA,
                linewidth=0.0, zorder=3,
            )
        axis.plot(
            mean, altitude, color=color, linewidth=mean_width,
            linestyle=linestyle, alpha=CO2_MEAN_ALPHA,
            zorder=mean_zorder,
        )
        # Alternate open markers between the two classes at actual layer-centre
        # heights.  This tags coincident curves without displacing either one.
        marker_heights = mean_marker_heights[marker_offset::2]
        marker_values = np.interp(marker_heights, altitude, mean)
        marker_sigma = np.interp(marker_heights, altitude, curve_sigma)
        # Number concentration is plotted logarithmically.  Guard against an
        # extreme partial ensemble producing a lower error-bar endpoint at or
        # below zero while preserving the requested symmetric sigma whenever
        # it is physically admissible.
        marker_sigma = np.minimum(marker_sigma, 0.95 * marker_values)
        axis.errorbar(
            marker_values, marker_heights, xerr=marker_sigma,
            fmt=marker, linestyle="none", markersize=4.2,
            markerfacecolor="white", markeredgecolor=color,
            markeredgewidth=1.15, ecolor=color, elinewidth=0.85,
            capsize=1.8, capthick=0.85, alpha=CO2_ERRORBAR_ALPHA,
            zorder=10,
        )

        psurfs = np.asarray([record["psurf"] for record in records])
        mean_psurf = float(np.mean(psurfs))
        psurf_lower, psurf_upper = np.percentile(psurfs, [16.0, 84.0])
        short_label = "unc" if retrieval_class == "uncorrected" else "corr"
        psurf_summary.append(
            "%s %.2f [%.2f, %.2f]" % (
                short_label, mean_psurf, psurf_lower, psurf_upper
            )
        )
        total_columns = np.asarray([
            profile["total_co2_column"] for profile in profiles
        ]) / 1.0e21
        column_summary.append(
            "%s %.3f" % (short_label, np.mean(total_columns))
        )
        xco2_values = np.asarray([
            record.get("saved_xco2_ppm", profile["xco2_ppm"])
            for record, profile in zip(records, profiles)
        ])
        xco2_summary[retrieval_class] = (
            float(np.mean(xco2_values)),
            float(np.std(xco2_values, ddof=1)),
        )
        if show_noiseless and ensemble["noiseless_paired"]:
            noiseless = co2_molecular_profile(
                ensemble["noiseless"][retrieval_class], atmosphere
            )
            noiseless_curve = interpolate_co2_profile(noiseless, altitude)
            axis.plot(
                noiseless_curve, altitude, color=color, linestyle=":",
                linewidth=0.85, alpha=CO2_NOISELESS_ALPHA, zorder=5,
            )
            axis.plot(
                noiseless_curve[0], altitude[0], marker="s",
                markerfacecolor="white",
                markeredgecolor=color, markersize=3.2,
                alpha=CO2_NOISELESS_ALPHA, zorder=9,
            )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(xlimits)
    axis.set_ylim(altitude[0], altitude[-1])
    axis.set_xlabel(
        "Layer-mean CO$_2$ number density (molecules cm$^{-3}$)",
        fontsize=FONT_AXIS_LABEL,
    )
    axis.set_ylabel(
        "Altitude (km; log scale)", fontsize=FONT_AXIS_LABEL,
    )
    visible_ticks = [
        value for value in CO2_ALTITUDE_TICKS_KM
        if altitude[0] <= value <= altitude[-1]
    ]
    axis.yaxis.set_major_locator(FixedLocator(visible_ticks))
    axis.yaxis.set_major_formatter(FuncFormatter(
        lambda value, _position: "%g" % value
    ))
    axis.tick_params(labelsize=FONT_AXIS_TICK)
    axis.grid(axis="x", which="both", alpha=0.18)
    axis.grid(axis="y", which="major", alpha=0.18)
    axis.grid(axis="y", which="minor", alpha=0.07, linewidth=0.45)
    unc_mean, unc_std = xco2_summary["uncorrected"]
    corr_mean, corr_std = xco2_summary["corrected"]
    xco2_annotation = (
        "$XCO_2$: truth %.4f | unc %.4f $\\pm$ %.4f | "
        "corr %.4f $\\pm$ %.4f ppm  (n=%d)" % (
            truth_record_value.get("saved_xco2_ppm", truth["xco2_ppm"]),
            unc_mean, unc_std, corr_mean, corr_std, count,
        )
    )
    axis.text(
        0.5, 1.018, xco2_annotation,
        transform=axis.transAxes, ha="center", va="bottom",
        fontsize=FONT_METADATA, color="0.22", clip_on=False,
    )
    axis.text(
        0.018, 0.025,
        "$p_{surf}$ (hPa): truth %.2f; %s; %s\n"
        "$\\Sigma N_{CO_2}$ ($10^{21}$ cm$^{-2}$): truth %.3f; %s; %s" % (
            truth_record_value["psurf"], psurf_summary[0], psurf_summary[1],
            truth["total_co2_column"] / 1.0e21,
            column_summary[0], column_summary[1]
        ),
        transform=axis.transAxes, ha="left", va="bottom",
        fontsize=FONT_METADATA,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white",
                  edgecolor="0.82", alpha=0.88),
        zorder=15,
    )
    annotate_status(axis, ensemble)
    axis.set_title(
        "%s — state %03d" % (surface.title(), state_index),
        fontsize=FONT_CARD_TITLE, pad=24,
    )


def draw_effect_card(profile_axis, surface_axes, card, altitude, boundaries,
                     profile_limits, surface_limits_by_band, show_individual,
                     show_noiseless):
    ensemble = card["ensemble"]
    count = len(ensemble["indices"])
    state_index = card["state_index"]
    if count < 2:
        draw_placeholder(profile_axis, surface_axes, state_index, ensemble)
        profile_axis.set_title(
            "%s — state %03d" % (card["row"]["surface"].title(), state_index),
            fontsize=FONT_CARD_TITLE, pad=8,
        )
        return

    model_top = float(boundaries[0])
    if profile_limits is None:
        profile_limits = effect_profile_limits(
            {"local": card}, altitude, model_top,
            include_noiseless=show_noiseless,
        )
    draw_layer_boundaries(profile_axis, boundaries, altitude[-1])
    profile_axis.axvline(0.0, color=TRUTH_COLOR, linewidth=1.1, alpha=0.7)
    for species, _, sigma, linestyle, linewidth in SPECIES:
        corrected = profile_stack(
            ensemble["corrected"], species, sigma, altitude, model_top
        )
        uncorrected = profile_stack(
            ensemble["uncorrected"], species, sigma, altitude, model_top
        )
        differences = corrected - uncorrected
        if show_individual:
            for difference in differences:
                profile_axis.plot(
                    difference, altitude, color=EFFECT_COLOR,
                    linestyle=linestyle, linewidth=0.65, alpha=0.09,
                    zorder=2,
                )
        mean, lower, upper = curve_statistics(differences)
        if count >= 3:
            profile_axis.fill_betweenx(
                altitude, lower, upper, color=EFFECT_COLOR, alpha=0.14,
                linewidth=0.0, zorder=3,
            )
        profile_axis.plot(
            mean, altitude, color=EFFECT_COLOR, linestyle=linestyle,
            linewidth=max(1.35, 0.76 * linewidth), zorder=6,
        )
        if show_noiseless and ensemble["noiseless_paired"]:
            difference = (
                aerosol_profile(
                    ensemble["noiseless"]["corrected"], species, sigma,
                    altitude, model_top,
                ) -
                aerosol_profile(
                    ensemble["noiseless"]["uncorrected"], species, sigma,
                    altitude, model_top,
                )
            )
            profile_axis.plot(
                difference, altitude, color=EFFECT_COLOR,
                linestyle=linestyle, linewidth=0.85, alpha=0.75,
                marker="s", markevery=75, markersize=2.0,
                markerfacecolor="white", zorder=5,
            )

    profile_axis.set_xlim(profile_limits)
    profile_axis.set_ylim(0.0, altitude[-1])
    profile_axis.set_xlabel(
        "$\\Delta(d\\tau_{760}/dz)$: corrected $-$ uncorrected (km$^{-1}$)",
        fontsize=FONT_AXIS_LABEL,
    )
    profile_axis.set_ylabel("Altitude (km)", fontsize=FONT_AXIS_LABEL)
    profile_axis.grid(axis="x", alpha=0.18)
    profile_axis.tick_params(labelsize=FONT_AXIS_TICK)
    annotate_status(profile_axis, ensemble)
    profile_axis.set_title(
        "%s — state %03d" % (card["row"]["surface"].title(), state_index),
        fontsize=FONT_CARD_TITLE, pad=8,
    )

    for band_axis, (state_band, _, band_label, _, _) in zip(
            surface_axes, BANDS):
        if surface_limits_by_band is None:
            surface_limits = effect_surface_limits(
                {"local": card}, state_band,
                include_noiseless=show_noiseless,
            )
        else:
            surface_limits = surface_limits_by_band[state_band]
        wavelength, _ = BAND_GRIDS[state_band]
        corrected = surface_stack(ensemble["corrected"], state_band)
        uncorrected = surface_stack(ensemble["uncorrected"], state_band)
        differences = corrected - uncorrected
        band_axis.axhline(0.0, color=TRUTH_COLOR, linewidth=1.0, alpha=0.7)
        if show_individual:
            for difference in differences:
                band_axis.plot(
                    wavelength, difference, color=EFFECT_COLOR,
                    linewidth=0.55, alpha=0.09, zorder=2,
                )
        mean, lower, upper = curve_statistics(differences)
        if count >= 3:
            band_axis.fill_between(
                wavelength, lower, upper, color=EFFECT_COLOR, alpha=0.14,
                linewidth=0.0, zorder=3,
            )
        band_axis.plot(
            wavelength, mean, color=EFFECT_COLOR, linewidth=1.55, zorder=6
        )
        if show_noiseless and ensemble["noiseless_paired"]:
            difference = (
                surface_curve(
                    ensemble["noiseless"]["corrected"], state_band
                ) -
                surface_curve(
                    ensemble["noiseless"]["uncorrected"], state_band
                )
            )
            band_axis.plot(
                wavelength, difference, color=EFFECT_COLOR,
                linestyle=":", linewidth=0.95, alpha=0.82, zorder=5,
            )
        band_axis.set_title(
            band_label, fontsize=FONT_COMPACT_TITLE, pad=3,
        )
        band_axis.set_ylim(surface_limits)
        band_axis.set_xlim(wavelength[0], wavelength[-1])
        band_axis.set_xlabel(
            "Wavelength (nm)", fontsize=FONT_COMPACT_LABEL,
        )
        band_axis.tick_params(labelsize=FONT_COMPACT_TICK)
        band_axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        band_axis.ticklabel_format(
            axis="y", style="sci", scilimits=(-2, 2), useOffset=False
        )
        band_axis.grid(alpha=0.16)
    surface_axes[0].set_ylabel(
        "$\\Delta$ reflectance", fontsize=FONT_COMPACT_LABEL
    )


def figure_legend_groups(effect=False, show_noiseless=False):
    """Return signal/statistic and aerosol-species legend rows separately."""
    if effect:
        signal_handles = [
            Line2D([], [], color=EFFECT_COLOR, linewidth=2.1,
                   label="Paired corrected $-$ uncorrected mean"),
            Patch(facecolor=EFFECT_COLOR, alpha=0.14,
                  label="16th--84th percentile"),
        ]
    else:
        signal_handles = [
            Line2D([], [], color=TRUTH_COLOR, linewidth=2.2, label="Truth"),
            Line2D([], [], color=UNCORRECTED_COLOR, linewidth=2.0,
                   label="Uncorrected mean"),
            Line2D([], [], color=CORRECTED_COLOR, linewidth=2.0,
                   label="Corrected mean"),
            Patch(facecolor="0.45", alpha=0.14,
                  label="16th--84th percentile"),
        ]
    if show_noiseless:
        signal_handles.append(Line2D(
            [], [], color="0.35", linestyle=":", marker="s",
            markerfacecolor="white", linewidth=1.0,
            label="Perturbation 11 (noiseless; excluded from statistics)",
        ))
    species_handles = [
        Line2D([], [], color="0.15", linestyle=linestyle,
               linewidth=linewidth, label=label)
        for _, label, _, linestyle, linewidth in SPECIES
    ]
    return signal_handles, species_handles


def make_figure(category, cards, co2_selection, boundaries, output, effect,
                scale_mode, show_individual, show_noiseless, dpi):
    altitude = np.linspace(0.001, 16.5, 700)
    model_top = float(boundaries[0])
    if scale_mode == "shared":
        if effect:
            profile_limits = effect_profile_limits(
                cards, altitude, model_top,
                include_noiseless=show_noiseless,
            )
            surface_limits_by_band = {
                state_band: effect_surface_limits(
                    cards, state_band, include_noiseless=show_noiseless
                )
                for state_band, _, _, _, _ in BANDS
            }
        else:
            profile_limits = main_profile_limits(
                cards, altitude, model_top,
                include_noiseless=show_noiseless,
            )
            surface_limits_by_band = {
                state_band: main_surface_limits(
                    cards, state_band, include_noiseless=show_noiseless
                )
                for state_band, _, _, _, _ in BANDS
            }
    else:
        profile_limits = None
        surface_limits_by_band = None

    has_campaign_label = bool(co2_selection["campaign_label"])
    plot_top = 0.845 if has_campaign_label else 0.875
    fig = plt.figure(figsize=(16.8, 12.4))
    outer = GridSpec(
        2, 2, figure=fig, left=0.07, right=0.985, bottom=0.125,
        top=plot_top, hspace=0.25, wspace=0.18,
    )
    for position, surface in enumerate(SURFACE_ORDER):
        row = position // 2
        column = position % 2
        inner = GridSpecFromSubplotSpec(
            2, 3, subplot_spec=outer[row, column],
            height_ratios=(1.22, 1.0), hspace=0.42, wspace=0.24,
        )
        card = cards[surface]
        profile_axis = fig.add_subplot(inner[0, :])
        surface_axes = [
            fig.add_subplot(inner[1, index]) for index in range(3)
        ]
        if effect:
            draw_effect_card(
                profile_axis, surface_axes, card, altitude, boundaries,
                profile_limits, surface_limits_by_band, show_individual,
                show_noiseless,
            )
        else:
            draw_main_card(
                profile_axis, surface_axes, card, altitude, boundaries,
                profile_limits, surface_limits_by_band, show_individual,
                show_noiseless,
            )

    category_label = CATEGORY_LABELS[category]
    if effect:
        title = (
            "Paired correction effect in physical retrieval products — "
            "%s, %s" % (category_label, co2_selection["sif_short"])
        )
        subtitle = (
            "Each curve is corrected minus uncorrected for the same noise draw; "
            "bands use separate wavelength axes"
        )
    else:
        title = (
            "Corrected and uncorrected physical retrieval products — "
            "%s, %s" % (category_label, co2_selection["sif_short"])
        )
        subtitle = (
            "Terminal-state curves; paired perturbations 01--10; "
            "bands use separate canonical wavelength axes"
        )
        if category == "no_aerosol":
            subtitle += "; the truth aerosol-density profile is zero"
    if scale_mode == "shared":
        subtitle += "; aerosol limits are shared by category and reflectance limits by band"
    else:
        subtitle += "; locally zoomed diagnostic axes"
    if co2_selection.get("state_description"):
        subtitle += "\n" + co2_selection["state_description"]
    if has_campaign_label:
        fig.text(
            0.5, 0.990, co2_selection["campaign_label"],
            ha="center", va="top", fontsize=FONT_CAMPAIGN, color="0.34",
            fontweight="semibold",
        )
    title_y = 0.968 if has_campaign_label else 0.988
    display_y = 0.940 if has_campaign_label else 0.957
    subtitle_y = 0.902 if has_campaign_label else 0.924
    fig.suptitle(title, fontsize=17, y=title_y)
    fig.text(
        0.5, display_y, co2_selection["display"], ha="center", va="center",
        fontsize=FONT_SELECTION,
    )
    fig.text(
        0.5, subtitle_y, subtitle, ha="center", va="center",
        fontsize=FONT_HEADER,
        linespacing=1.22,
    )
    signal_handles, species_handles = figure_legend_groups(
        effect=effect, show_noiseless=show_noiseless,
    )
    fig.legend(
        handles=signal_handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.050),
        ncol=len(signal_handles), frameon=False,
        fontsize=FONT_AEROSOL_LEGEND,
        handlelength=2.8, columnspacing=1.35,
    )
    fig.legend(
        handles=species_handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.019),
        ncol=len(species_handles), frameon=False,
        fontsize=FONT_AEROSOL_LEGEND,
        handlelength=3.0, columnspacing=1.7,
    )
    snapshot = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    fig.text(
        0.985, 0.004, "Availability snapshot: %s" % snapshot,
        ha="right", va="bottom", fontsize=FONT_SNAPSHOT, color="0.42",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output), dpi=dpi, facecolor="white")
    plt.close(fig)


def co2_figure_legend_groups(show_noiseless=False, round4=False, round5=False,
                           round6=False):
    """Return physical-state and diagnostic legend rows separately."""
    truth_label = (
        "Truth (CO$_2$); round-4 state-space truth (SIF)"
        if round4 else "Truth (CO$_2$ points; SIF line)"
    )
    if round5:
        truth_label = "Truth (CO$_2$); round-5 fixed state (SIF)"
    if round6:
        truth_label = "Truth (CO$_2$); round-6 fixed state (SIF)"
    state_handles = [
        Line2D([], [], color=TRUTH_COLOR, linewidth=0.90, marker="o",
               markersize=2.8, markerfacecolor="none", markeredgewidth=0.7,
               alpha=CO2_TRUTH_MARKER_ALPHA,
               label=truth_label),
        Line2D([], [], color=CO2_UNCORRECTED_COLOR, linewidth=2.20,
               linestyle="-", marker="s", markersize=4.2,
               markerfacecolor="white", markeredgewidth=1.15,
               label="Uncorrected CO$_2$ mean"),
        Line2D([], [], color=CO2_CORRECTED_COLOR, linewidth=2.20,
               linestyle=(0, (3.0, 1.8)), marker="^", markersize=4.2,
               markerfacecolor="white", markeredgewidth=1.15,
               label="Corrected CO$_2$ mean"),
        Line2D([], [], color=UNCORRECTED_COLOR, linewidth=2.05,
               label="Uncorrected SIF mean"),
        Line2D([], [], color=CORRECTED_COLOR, linewidth=1.35,
               label="Corrected SIF mean"),
    ]
    diagnostic_handles = [
        Patch(facecolor="0.45", alpha=CO2_ENVELOPE_ALPHA,
              label="Noise-ensemble 16th--84th percentile cloud"),
        Line2D([], [], color="0.35", linewidth=0.9, marker="o",
               markersize=3.8, markerfacecolor="white",
               label="Noise-ensemble layer mean $\\pm$ 1 sample SD"),
        Line2D([], [], color="0.35", linestyle="--", linewidth=0.9,
               label="Layers 1--4 / active-state boundary"),
    ]
    if show_noiseless:
        diagnostic_handles.append(Line2D(
            [], [], color="0.35", linestyle=":", marker="s",
            markerfacecolor="white", linewidth=1.0,
            label="Perturbation 11 (noiseless; excluded from statistics)",
        ))
    return state_handles, diagnostic_handles


def make_co2_figure(category, cards, co2_selection, atmosphere, output,
                    scale_mode, show_individual, show_noiseless, dpi):
    """Create a 2x2 surface comparison of CO2 profiles and SIF spectra."""
    altitude = np.geomspace(
        CO2_ALTITUDE_MIN_KM, CO2_ALTITUDE_MAX_KM, 1800
    )
    shared_co2_limits = None
    shared_sif_limits = None
    if scale_mode == "shared":
        shared_co2_limits = co2_profile_limits(
            cards, atmosphere, altitude, include_noiseless=show_noiseless
        )
        shared_sif_limits = sif_limits(
            cards, include_noiseless=show_noiseless
        )

    has_campaign_label = bool(co2_selection["campaign_label"])
    plot_top = 0.805 if has_campaign_label else 0.835
    fig = plt.figure(figsize=(14.4, 15.1))
    outer = GridSpec(
        2, 2, figure=fig, left=0.078, right=0.98, bottom=0.135,
        top=plot_top, hspace=0.26, wspace=0.20,
    )
    for position, surface in enumerate(SURFACE_ORDER):
        inner = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[position // 2, position % 2],
            height_ratios=(2.35, 1.34), hspace=0.58,
        )
        sif_block = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=inner[1, 0],
            height_ratios=(1.0, 0.34), hspace=0.28,
        )
        co2_axis = fig.add_subplot(inner[0, 0])
        sif_axis = fig.add_subplot(sif_block[0, 0])
        sif_annotation_axis = fig.add_subplot(sif_block[1, 0])
        draw_co2_card(
            co2_axis, cards[surface], atmosphere, altitude,
            shared_co2_limits,
            show_individual, show_noiseless
        )
        draw_sif_card(
            sif_axis, sif_annotation_axis, cards[surface], shared_sif_limits,
            show_individual, show_noiseless,
        )

    category_label = CATEGORY_LABELS[category]
    title = (
        "CO$_2$ molecular concentration and O$_2$ A-band SIF state — "
        "%s, %s" % (category_label, co2_selection["sif_truth"])
    )
    if has_campaign_label:
        fig.text(
            0.5, 0.990, co2_selection["campaign_label"],
            ha="center", va="top", fontsize=FONT_CAMPAIGN, color="0.34",
            fontweight="semibold",
        )
    title_y = 0.968 if has_campaign_label else 0.988
    display_y = 0.940 if has_campaign_label else 0.957
    description_y = 0.875 if has_campaign_label else 0.908
    fig.suptitle(title, fontsize=15.5, y=title_y)
    fig.text(
        0.5, display_y, co2_selection["display"], ha="center", va="center",
        fontsize=FONT_SELECTION,
    )
    scale_label = (
        "CO$_2$/SIF limits shared by category" if scale_mode == "shared" else
        "locally zoomed CO$_2$ and SIF axes"
    )
    state_label = co2_selection.get("state_description")
    scale_and_state = scale_label
    if state_label:
        scale_and_state += "; " + state_label
    fig.text(
        0.5, description_y,
        "Layer CO$_2$ column = VMR $\\times$ dry-air VCD; concentration = "
        "column / layer thickness; log-PCHIP through layer centers; "
        "logarithmic altitude emphasizes the lower atmosphere\n"
        "SIF is the wavenumber-linear state converted to radiance "
        "per nm; card annotations use matched perturbations 01--10 and report "
        "$XCO_2$ mean $\\pm$ sample SD\n"
        "Dashed line marks the layers 1--4 / active-state boundary; the fixed "
        "profile continues above the plotted range\n%s" % scale_and_state,
        ha="center", va="center", fontsize=FONT_HEADER, linespacing=1.22,
    )
    state_handles, diagnostic_handles = co2_figure_legend_groups(
        show_noiseless=show_noiseless,
        round4=bool(co2_selection.get("state_description")),
        round5=(co2_selection.get("state_description") or "").startswith("round 5:"),
        round6=(co2_selection.get("state_description") or "").startswith("round 6:"),
    )
    fig.legend(
        handles=state_handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.052),
        ncol=len(state_handles), frameon=False, fontsize=FONT_CO2_LEGEND,
        handlelength=2.5, columnspacing=1.05,
    )
    fig.legend(
        handles=diagnostic_handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.019),
        ncol=len(diagnostic_handles), frameon=False,
        fontsize=FONT_CO2_LEGEND,
        handlelength=2.5, columnspacing=1.15,
    )
    snapshot = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    fig.text(
        0.98, 0.004, "Availability snapshot: %s" % snapshot,
        ha="right", va="bottom", fontsize=FONT_SNAPSHOT, color="0.42",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output), dpi=dpi, facecolor="white")
    plt.close(fig)


def print_inventory(categories):
    print("state surface category paired_noisy corrected_valid uncorrected_valid p11")
    for category in CATEGORY_AEROSOL_CASES:
        if category not in categories:
            continue
        for surface in SURFACE_ORDER:
            card = categories[category][surface]
            ensemble = card["ensemble"]
            print(
                "%03d %-7s %-12s %2d/10 %2d/10 %2d/10 %s" % (
                    card["state_index"], surface, category,
                    len(ensemble["indices"]),
                    ensemble["valid_counts"]["corrected"],
                    ensemble["valid_counts"]["uncorrected"],
                    "paired" if ensemble["noiseless_paired"] else "missing",
                )
            )


def main():
    args = parse_args()
    if args.co2_only and args.product != "all":
        raise ValueError("--co2-only cannot be combined with --product")
    product = "co2-sif" if args.co2_only else args.product
    categories, boundaries, atmosphere, co2_selection = prepare_categories(args)
    print_inventory(categories)
    show_individual = not args.hide_individual
    outputs = []
    co2_value = co2_selection["value"]
    co2_tag = (
        str(int(round(co2_value)))
        if math.isclose(co2_value, round(co2_value), abs_tol=1.0e-9)
        else ("%g" % co2_value).replace(".", "p")
    )
    for category in co2_selection["category_names"]:
        base = "%s_%s_%s_%s" % (
            co2_selection["slug"], co2_tag,
            co2_selection["sif_slug"], category
        )
        if product in ("all", "aerosol-surface"):
            plot_specs = [("profiles_surface", False)]
            if args.include_paired_correction:
                plot_specs.append(("paired_correction_effect", True))
            for suffix, effect in plot_specs:
                output = args.output_dir / ("%s_%s.png" % (base, suffix))
                make_figure(
                    category, categories[category], co2_selection,
                    boundaries, output, effect, args.scale_mode,
                    show_individual, args.show_noiseless, args.dpi,
                )
                outputs.append(output)
        if product in ("all", "co2-sif"):
            co2_output = args.output_dir / ("%s_co2_profiles.png" % base)
            make_co2_figure(
                category, categories[category], co2_selection, atmosphere,
                co2_output,
                args.scale_mode, show_individual, args.show_noiseless, args.dpi,
            )
            outputs.append(co2_output)
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
