#!/usr/bin/env python3
"""Central physical-ensemble plotter for the RRS/XCO2 retrieval regimes.

The required positional index selects a complete retrieval interpretation and
its associated truth, prior family, and output namespace:

1. full-column CO2 retrievals;
2. tightly correlated bottom-layer CO2 retrievals;
3. loosely correlated bottom-layer CO2 retrievals with SIF parameters
   co-retrieved against either no-SIF or corrected SIF-on truth;
4. loosely correlated bottom-layer CO2 retrievals with the round-4 reduced
   SIF state (known SIF at 759 nm), again for either truth category;
5. round-5 fixed SIF at 759 nm and fixed spectral slope, with tighter UTLS
   aerosol priors and the unchanged loose CO2 prior, for either truth category;
6. round-6 fixed SIF, with the original UTLS and CO2 priors restored/preserved.

``--SIF=0`` selects no-SIF truth; ``--SIF=1`` selects corrected SIF-on truth.

The plotting backend remains shared with the established inversion workflow.
This entry point owns campaign routing so a bottom-layer result cannot be
silently interpreted with the full-column truth table.
"""

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import runpy
import sys

from netCDF4 import Dataset


HERE = Path(__file__).resolve().parent
RRS_ROOT = HERE.parent
INVERSION_ROOT = RRS_ROOT / "inversion"
BOTTOM_ROOT = RRS_ROOT / "bottom_layer_XCO2_retrievals"
BACKEND = INVERSION_ROOT / "plot_physical_retrieval_ensembles.py"
VERTICAL_PROFILE_TABLE = (
    RRS_ROOT / "truth_map_aerosols" / "aerosol_vertical_profiles.dat"
)
ATMOSPHERIC_PROFILE = (
    INVERSION_ROOT / "retrieval_setup" / "retrieval_atmosphere_16layer.nc"
)
PRIVATE_RESULTS_ROOT = Path(os.environ.get(
    "RRS_XCO2_PRIVATE_RESULTS_ROOT",
    str(Path.home() / "RRS_XCO2_private" / "results"),
))
ROUND3_SIF_CAMPAIGN = (
    PRIVATE_RESULTS_ROOT /
    "bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1"
)
ROUND4_SIF_CAMPAIGN = (
    PRIVATE_RESULTS_ROOT /
    "bottom_layer_round4_known_sif759_sif_on_"
    "acos_mapped_tapered_vertical_correlation_v1"
)
ROUND5_SIF_CAMPAIGN = (
    PRIVATE_RESULTS_ROOT /
    "bottom_layer_round5_fixed_sif_on_tight_utls_"
    "acos_mapped_tapered_vertical_correlation_v1"
)
ROUND6_SIF_CAMPAIGN = (
    PRIVATE_RESULTS_ROOT /
    "bottom_layer_round6_fixed_sif_on_standard_utls_"
    "acos_mapped_tapered_vertical_correlation_v1"
)
ROUND3_SIF_PLOT_ROOT = (
    BOTTOM_ROOT /
    "retrievals_acos_mapped_tapered_vertical_correlation_sif" /
    "plots" / "physical_ensembles"
)
ROUND4_SIF_PLOT_ROOT = (
    BOTTOM_ROOT / "round4_known_sif759" / "plots_sif" /
    "physical_ensembles"
)
CORRECTED_SIF_TRUTH_TABLE = (
    ROUND3_SIF_CAMPAIGN / "retrieval_setup" /
    "true_states_corrected_sif_v2.dat"
)
SIF_CASE_ON = "angular_integral760_0p5"


@dataclass(frozen=True)
class RetrievalRegime:
    index: int
    name: str
    label: str
    retrieval_root: Path
    truth_table: Path
    scene_components: Path
    output_root: Path
    co2_coordinate: str
    sif_case: str = "off"
    source_prior_basename: str = ""
    state_model: str = "legacy"


REGIMES = {
    1: RetrievalRegime(
        1,
        "full-column",
        "Type 1: full-column CO2",
        INVERSION_ROOT,
        RRS_ROOT / "truth_map" / "true_states.dat",
        RRS_ROOT / "truth_map" / "scene_components.dat",
        INVERSION_ROOT / "physical_ensemble_visualizations" / "full_column",
        "xco2",
    ),
    2: RetrievalRegime(
        2,
        "bottom-layer-tight",
        "Type 2: tightly constrained bottom-layer CO2",
        BOTTOM_ROOT / "retrievals",
        BOTTOM_ROOT / "truth" / "true_states.dat",
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        BOTTOM_ROOT / "retrievals" / "plots" / "physical_ensembles",
        "bottom_co2",
        source_prior_basename="apriori_states.nc",
    ),
    3: RetrievalRegime(
        3,
        "bottom-layer-loose",
        (
            "Type 3: loosely constrained bottom-layer CO2; "
            "SIF parameters co-retrieved"
        ),
        (
            BOTTOM_ROOT /
            "retrievals_acos_mapped_tapered_vertical_correlation_nosif"
        ),
        BOTTOM_ROOT / "truth" / "true_states.dat",
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        (
            BOTTOM_ROOT /
            "retrievals_acos_mapped_tapered_vertical_correlation_nosif" /
            "plots" / "physical_ensembles"
        ),
        "bottom_co2",
        source_prior_basename=(
            "apriori_states_acos_mapped_tapered_vertical_correlation.nc"
        ),
    ),
    4: RetrievalRegime(
        4,
        "bottom-layer-loose-reduced-sif",
        "Type 4: loose bottom-layer CO2; no SIF intercept co-retrieval",
        BOTTOM_ROOT / "round4_known_sif759" / "retrievals_nosif",
        BOTTOM_ROOT / "truth" / "true_states.dat",
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        (
            BOTTOM_ROOT / "round4_known_sif759" / "plots_nosif" /
            "physical_ensembles"
        ),
        "bottom_co2",
        source_prior_basename=(
            "apriori_states_round4_known_sif759_off_"
            "acos_mapped_tapered_vertical_correlation.nc"
        ),
        state_model="round4_known_sif759",
    ),
    5: RetrievalRegime(
        5,
        "bottom-layer-fixed-sif-tight-utls",
        "Round 5: loose bottom-layer CO2; fixed SIF; tight UTLS aerosol prior",
        BOTTOM_ROOT / "round5_fixed_sif" / "retrievals_nosif",
        BOTTOM_ROOT / "truth" / "true_states.dat",
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        BOTTOM_ROOT / "round5_fixed_sif" / "plots_nosif" / "physical_ensembles",
        "bottom_co2",
        source_prior_basename=(
            "apriori_states_round5_fixed_sif_off_tight_utls_"
            "acos_mapped_tapered_vertical_correlation.nc"
        ),
        state_model="round5_fixed_sif",
    ),
    6: RetrievalRegime(
        6,
        "bottom-layer-fixed-sif-standard-utls",
        "Round 6: loose bottom-layer CO2; fixed SIF; original UTLS aerosol prior",
        BOTTOM_ROOT / "round6_fixed_sif" / "retrievals_nosif",
        BOTTOM_ROOT / "truth" / "true_states.dat",
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        BOTTOM_ROOT / "round6_fixed_sif" / "plots_nosif" / "physical_ensembles",
        "bottom_co2",
        source_prior_basename=(
            "apriori_states_round6_fixed_sif_off_standard_utls_"
            "acos_mapped_tapered_vertical_correlation.nc"
        ),
        state_model="round6_fixed_sif",
    ),
}


SIF_ON_REGIMES = {
    3: RetrievalRegime(
        3,
        "bottom-layer-loose-sif-on",
        "Type 3: loosely constrained bottom-layer CO2; SIF-on truth",
        ROUND3_SIF_CAMPAIGN / "retrievals",
        CORRECTED_SIF_TRUTH_TABLE,
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        ROUND3_SIF_PLOT_ROOT,
        "bottom_co2",
        sif_case=SIF_CASE_ON,
        source_prior_basename=(
            "apriori_states_acos_mapped_tapered_vertical_correlation.nc"
        ),
    ),
    4: RetrievalRegime(
        4,
        "bottom-layer-loose-reduced-sif-on",
        "Type 4: loose bottom-layer CO2; known SIF at 759 nm",
        ROUND4_SIF_CAMPAIGN / "retrievals",
        CORRECTED_SIF_TRUTH_TABLE,
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        ROUND4_SIF_PLOT_ROOT,
        "bottom_co2",
        sif_case=SIF_CASE_ON,
        source_prior_basename=(
            "apriori_states_round4_known_sif759_on_"
            "acos_mapped_tapered_vertical_correlation.nc"
        ),
        state_model="round4_known_sif759",
    ),
    5: RetrievalRegime(
        5,
        "bottom-layer-fixed-sif-tight-utls-on",
        "Round 5: loose bottom-layer CO2; fixed SIF; tight UTLS aerosol prior",
        ROUND5_SIF_CAMPAIGN / "retrievals",
        CORRECTED_SIF_TRUTH_TABLE,
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        BOTTOM_ROOT / "round5_fixed_sif" / "plots_sif" / "physical_ensembles",
        "bottom_co2",
        sif_case=SIF_CASE_ON,
        source_prior_basename=(
            "apriori_states_round5_fixed_sif_on_tight_utls_"
            "acos_mapped_tapered_vertical_correlation.nc"
        ),
        state_model="round5_fixed_sif",
    ),
    6: RetrievalRegime(
        6,
        "bottom-layer-fixed-sif-standard-utls-on",
        "Round 6: loose bottom-layer CO2; fixed SIF; original UTLS aerosol prior",
        ROUND6_SIF_CAMPAIGN / "retrievals",
        CORRECTED_SIF_TRUTH_TABLE,
        BOTTOM_ROOT / "truth" / "scene_components.dat",
        BOTTOM_ROOT / "round6_fixed_sif" / "plots_sif" / "physical_ensembles",
        "bottom_co2",
        sif_case=SIF_CASE_ON,
        source_prior_basename=(
            "apriori_states_round6_fixed_sif_on_standard_utls_"
            "acos_mapped_tapered_vertical_correlation.nc"
        ),
        state_model="round6_fixed_sif",
    ),
}


LOCKED_INPUT_OPTIONS = (
    "--inversion-root",
    "--truth-table",
    "--scene-components",
    "--vertical-profile-table",
    "--atmospheric-profile",
    "--campaign-label",
    "--sif-case",
)


def resolve_regime(index, sif_enabled):
    """Select one protected campaign from the retrieval type and SIF flag."""
    if not sif_enabled:
        return REGIMES[index]
    if index not in SIF_ON_REGIMES:
        raise ValueError(
            "--SIF=1 is available for retrieval types 3 through 6 only; "
            "types 1 and 2 do not have registered corrected-SIF campaigns"
        )
    return SIF_ON_REGIMES[index]


def has_option(arguments, option):
    return any(
        value == option or value.startswith(option + "=")
        for value in arguments
    )


def first_completed_retrieval(root):
    for retrieval_class in ("corrected", "uncorrected"):
        for path in sorted(
                (root / retrieval_class).glob(
                    "retrieval_state*_perturbation*.nc"
                )):
            try:
                with Dataset(path) as dataset:
                    if int(dataset.getncattr("retrieval_complete")) == 1:
                        return path
            except (OSError, AttributeError, ValueError):
                continue
    raise RuntimeError(
        "no completed retrieval was found under %s" % root
    )


def validate_regime(regime):
    """Fail closed if an index points at the wrong campaign family."""
    required = (
        regime.truth_table,
        regime.scene_components,
        VERTICAL_PROFILE_TABLE,
        ATMOSPHERIC_PROFILE,
        BACKEND,
    )
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(
            "retrieval type %d is missing required input(s): %s" % (
                regime.index, ", ".join(str(path) for path in missing)
            )
        )
    sample = first_completed_retrieval(regime.retrieval_root)
    with Dataset(sample) as dataset:
        has_bottom_truth = "truth_bottom_co2_ppm" in dataset.ncattrs()
        if regime.co2_coordinate == "xco2" and has_bottom_truth:
            raise RuntimeError(
                "type %d unexpectedly points at bottom-layer retrievals" %
                regime.index
            )
        if regime.co2_coordinate == "bottom_co2" and not has_bottom_truth:
            raise RuntimeError(
                "type %d does not point at bottom-layer retrievals" %
                regime.index
            )
        actual_model = (
            str(dataset.getncattr("retrieval_state_model"))
            if "retrieval_state_model" in dataset.ncattrs() else "legacy"
        )
        if actual_model != regime.state_model:
            raise RuntimeError(
                "type %d expected state model %s but %s uses %s" % (
                    regime.index, regime.state_model, sample, actual_model
                )
            )
        actual_sif_case = (
            str(dataset.getncattr("sif_case"))
            if "sif_case" in dataset.ncattrs() else "off"
        )
        if actual_sif_case != regime.sif_case:
            raise RuntimeError(
                "type %d expected SIF case %s but %s records %s" % (
                    regime.index, regime.sif_case, sample, actual_sif_case,
                )
            )
        if regime.source_prior_basename:
            if "source_apriori" not in dataset.ncattrs():
                raise RuntimeError(
                    "%s is missing source_apriori metadata" % sample
                )
            actual_prior = Path(
                str(dataset.getncattr("source_apriori"))
            ).name
            if actual_prior != regime.source_prior_basename:
                raise RuntimeError(
                    "type %d expected prior %s but %s records %s" % (
                        regime.index, regime.source_prior_basename,
                        sample, actual_prior,
                    )
                )
    return sample


def parse_frontend(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "All remaining options are passed to the physical-ensemble "
            "backend.\nUse INDEX --backend-help to display those options."
        ),
    )
    parser.add_argument(
        "retrieval_type", type=int, choices=tuple(REGIMES),
        help=(
            "1=full column, 2=tight bottom layer, 3=loose bottom layer, "
            "4=loose bottom layer with reduced/fixed SIF state, "
            "5=fixed SIF and tight UTLS aerosol prior, "
            "6=fixed SIF and original UTLS aerosol prior"
        ),
    )
    parser.add_argument(
        "--backend-help", action="store_true",
        help="Show all physical-ensemble selection and styling options",
    )
    parser.add_argument(
        "--SIF", "--sif", dest="sif", type=int, choices=(0, 1), default=0,
        help=(
            "Select the truth SIF category: 0=no SIF (default), "
            "1=corrected SIF-on truth"
        ),
    )
    frontend, remaining = parser.parse_known_args(arguments)
    return parser, frontend, remaining


def main():
    parser, frontend, remaining = parse_frontend(sys.argv[1:])
    if frontend.backend_help:
        sys.argv = [str(BACKEND), "--help"]
        sys.path.insert(0, str(INVERSION_ROOT))
        runpy.run_path(str(BACKEND), run_name="__main__")
        return

    forbidden = [
        option for option in LOCKED_INPUT_OPTIONS
        if has_option(remaining, option)
    ]
    if forbidden:
        parser.error(
            "retrieval type controls campaign inputs; remove %s" %
            ", ".join(forbidden)
        )

    selection = argparse.ArgumentParser(add_help=False)
    selection.add_argument("--bottom-co2-ppm")
    selection.add_argument("--xco2-ppm")
    selection.add_argument("--output-dir")
    selected, _ = selection.parse_known_args(remaining)
    try:
        regime = resolve_regime(frontend.retrieval_type, frontend.sif)
    except ValueError as error:
        parser.error(str(error))
    if regime.co2_coordinate == "xco2":
        if selected.bottom_co2_ppm is not None:
            parser.error(
                "type 1 is full-column; use --xco2-ppm, not "
                "--bottom-co2-ppm"
            )
        co2_default = ("--xco2-ppm", 400)
    else:
        if selected.xco2_ppm is not None:
            parser.error(
                "types 2--6 are bottom-layer experiments; use "
                "--bottom-co2-ppm, not --xco2-ppm"
            )
        co2_default = ("--bottom-co2-ppm", 400)

    sample = validate_regime(regime)
    defaults = (
        ("--inversion-root", regime.retrieval_root),
        ("--truth-table", regime.truth_table),
        ("--scene-components", regime.scene_components),
        ("--vertical-profile-table", VERTICAL_PROFILE_TABLE),
        ("--atmospheric-profile", ATMOSPHERIC_PROFILE),
        ("--output-dir", regime.output_root),
        ("--campaign-label", regime.label),
        co2_default,
        ("--sif-case", regime.sif_case),
    )
    injected = []
    for option, value in defaults:
        if not has_option(remaining, option):
            injected.extend((option, str(value)))

    effective_output = (
        Path(selected.output_dir) if selected.output_dir is not None
        else regime.output_root
    )
    print(
        "retrieval_type=%d (%s)\nSIF=%d\ninput_root=%s\n"
        "sample=%s\noutput_dir=%s" % (
            regime.index, regime.name, frontend.sif, regime.retrieval_root,
            sample, effective_output,
        )
    )
    sys.argv = [str(BACKEND)] + injected + remaining
    sys.path.insert(0, str(INVERSION_ROOT))
    runpy.run_path(str(BACKEND), run_name="__main__")


if __name__ == "__main__":
    main()
