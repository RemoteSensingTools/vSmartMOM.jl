#!/usr/bin/env python3
"""Plot one explicitly selected bottom-layer retrieval campaign/category.

The shared plotting implementation also supports full-column retrievals.  This
entry point removes that ambiguity by supplying bottom-layer truth/profile
inputs and named presets for the local round-3 and round-4 no-SIF campaigns.
Imported SIF-on campaigns use ``--campaign custom`` with explicit retrieval
and corrected-v2 truth paths.
"""

import argparse
import sys

from _adapter import (
    ATMOSPHERIC_PROFILE,
    CAMPAIGN_ROOT,
    SCENE_COMPONENTS,
    TRUTH_TABLE,
    VERTICAL_PROFILE_TABLE,
    run_inversion_script,
)


CAMPAIGN_PRESETS = {
    "round3-nosif": {
        "inversion_root": (
            CAMPAIGN_ROOT /
            "retrievals_acos_mapped_tapered_vertical_correlation_nosif"
        ),
        "output_dir": (
            CAMPAIGN_ROOT /
            "retrievals_acos_mapped_tapered_vertical_correlation_nosif" /
            "plots" / "physical_ensembles"
        ),
        "label": "Round 3 bottom-layer",
    },
    "round4-nosif": {
        "inversion_root": (
            CAMPAIGN_ROOT / "round4_known_sif759" / "retrievals_nosif"
        ),
        "output_dir": (
            CAMPAIGN_ROOT / "round4_known_sif759" / "plots_nosif" /
            "physical_ensembles"
        ),
        "label": "Round 4 bottom-layer",
    },
}


def wrapper_arguments(arguments):
    """Remove the wrapper-only campaign option and validate its preset."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--campaign",
        choices=tuple(CAMPAIGN_PRESETS) + ("custom",),
        default="round4-nosif",
    )
    wrapper, remaining = parser.parse_known_args(arguments)

    guard = argparse.ArgumentParser(add_help=False)
    guard.add_argument("--sif-case", default="off")
    guard.add_argument("--inversion-root")
    guard.add_argument("--truth-table")
    guard.add_argument("--output-dir")
    guard.add_argument("--campaign-label")
    selected, _ = guard.parse_known_args(remaining)

    if wrapper.campaign in CAMPAIGN_PRESETS and selected.sif_case != "off":
        parser.error(
            "%s contains no-SIF retrievals only; use --campaign custom with "
            "explicit --inversion-root, --truth-table, and --output-dir for "
            "an imported SIF-on campaign" % wrapper.campaign
        )
    if wrapper.campaign == "custom":
        missing = [
            option for option, value in (
                ("--inversion-root", selected.inversion_root),
                ("--truth-table", selected.truth_table),
                ("--output-dir", selected.output_dir),
                ("--campaign-label", selected.campaign_label),
            ) if value is None
        ]
        if missing:
            parser.error(
                "--campaign custom requires %s" % ", ".join(missing)
            )
    return wrapper.campaign, remaining


if __name__ == "__main__":
    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print(
            "Bottom-layer wrapper option:\n"
            "  --campaign {round3-nosif,round4-nosif,custom}\n"
            "      Select a local campaign preset (default: round4-nosif), or\n"
            "      use custom with explicit retrieval/truth/output paths.\n"
        )
    campaign, passthrough = wrapper_arguments(sys.argv[1:])
    sys.argv = [sys.argv[0]] + passthrough
    campaign_defaults = ()
    if campaign in CAMPAIGN_PRESETS:
        preset = CAMPAIGN_PRESETS[campaign]
        campaign_defaults = (
            ("--inversion-root", preset["inversion_root"]),
            ("--output-dir", preset["output_dir"]),
            ("--campaign-label", preset["label"]),
        )
    run_inversion_script(
        "plot_physical_retrieval_ensembles.py",
        campaign_defaults + (
            ("--truth-table", TRUTH_TABLE),
            ("--scene-components", SCENE_COMPONENTS),
            ("--vertical-profile-table", VERTICAL_PROFILE_TABLE),
            ("--atmospheric-profile", ATMOSPHERIC_PROFILE),
            ("--bottom-co2-ppm", 400),
        ),
    )
