#!/usr/bin/env python3
"""Central physical-ensemble table writer for the RRS/XCO2 regimes.

This is the tabular twin of ``plot_physical_retrieval_ensembles.py`` in this
directory and takes the same positional retrieval-type index, with the same
protection: the index selects the retrieval files, truth table, prior family,
state-vector interpretation, campaign label, and output namespace as one unit.

Tables are written to a ``tables*`` directory beside the corresponding
``plots*`` figure directory, so a Markdown table always sits next to the figure
drawn from the same retrieval files.

Unlike the plotter, which defaults to one CO2 case, this tool sweeps every
truth CO2 case of the selected campaign unless one case is named explicitly.
Cases with no paired retrievals yet are reported as pending and skipped, so the
command can simply be re-run as a campaign fills in.
"""

import argparse
import importlib.util
from pathlib import Path
import runpy
import sys


HERE = Path(__file__).resolve().parent
RRS_ROOT = HERE.parent
INVERSION_ROOT = RRS_ROOT / "inversion"
BACKEND = INVERSION_ROOT / "tabulate_physical_retrieval_ensembles.py"
PLOT_FRONTEND = HERE / "plot_physical_retrieval_ensembles.py"


def load_plot_frontend():
    """Reuse the plotting launcher's campaign registry and guards verbatim."""
    if not PLOT_FRONTEND.is_file():
        raise RuntimeError(
            "the physical-ensemble plotting launcher is missing: %s" %
            PLOT_FRONTEND
        )
    spec = importlib.util.spec_from_file_location(
        "_physical_ensemble_plot_frontend", PLOT_FRONTEND
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


router = load_plot_frontend()
REGIMES = router.REGIMES
LOCKED_INPUT_OPTIONS = router.LOCKED_INPUT_OPTIONS


def table_root(regime):
    """Return the ``tables*`` directory beside a regime's figure directory.

    Types 2--6 write figures to ``<root>/plots[_suffix]/physical_ensembles``,
    so the sibling becomes ``<root>/tables[_suffix]/physical_ensembles``.  The
    historical type-1 figure directory is not named ``plots``; its tables get
    the matching ``physical_ensemble_tables`` name instead.
    """
    output_root = regime.output_root
    parent = output_root.parent
    if parent.name == "plots" or parent.name.startswith("plots_"):
        return parent.parent / ("tables" + parent.name[len("plots"):]) / \
            output_root.name
    if parent.name == "physical_ensemble_visualizations":
        return parent.parent / "physical_ensemble_tables" / output_root.name
    return parent.parent / (parent.name + "_tables") / output_root.name


def parse_frontend(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "All remaining options are passed to the table backend.\n"
            "Use INDEX --backend-help to display those options."
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
        help="Show all table selection options",
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
        if router.has_option(remaining, option)
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
        regime = router.resolve_regime(frontend.retrieval_type, frontend.sif)
    except ValueError as error:
        parser.error(str(error))
    if regime.co2_coordinate == "xco2" and selected.bottom_co2_ppm is not None:
        parser.error("type 1 is full-column; use --xco2-ppm, not "
                     "--bottom-co2-ppm")
    if regime.co2_coordinate == "bottom_co2" and selected.xco2_ppm is not None:
        parser.error("types 2--6 are bottom-layer experiments; use "
                     "--bottom-co2-ppm, not --xco2-ppm")

    sample = router.validate_regime(regime)
    output_root = table_root(regime)
    defaults = [
        ("--inversion-root", regime.retrieval_root),
        ("--truth-table", regime.truth_table),
        ("--scene-components", regime.scene_components),
        ("--vertical-profile-table", router.VERTICAL_PROFILE_TABLE),
        ("--atmospheric-profile", router.ATMOSPHERIC_PROFILE),
        ("--output-dir", output_root),
        ("--campaign-label", regime.label),
        ("--sif-case", regime.sif_case),
    ]
    injected = []
    for option, value in defaults:
        if not router.has_option(remaining, option):
            injected.extend((option, str(value)))
    # Sweeping every truth CO2 case is the useful default for a table set; an
    # explicit case selector narrows it back to one case.
    sweep = (
        selected.bottom_co2_ppm is None and selected.xco2_ppm is None and
        not router.has_option(remaining, "--all-co2-cases")
    )
    if sweep:
        injected.append("--all-co2-cases")

    effective_output = (
        Path(selected.output_dir) if selected.output_dir is not None
        else output_root
    )
    print(
        "retrieval_type=%d (%s)\nSIF=%d\ninput_root=%s\nsample=%s\noutput_dir=%s\n"
        "co2_cases=%s" % (
            regime.index, regime.name, frontend.sif, regime.retrieval_root,
            sample, effective_output, "all" if sweep else "explicit",
        )
    )
    sys.argv = [str(BACKEND)] + injected + remaining
    sys.path.insert(0, str(INVERSION_ROOT))
    runpy.run_path(str(BACKEND), run_name="__main__")


if __name__ == "__main__":
    main()
