#!/usr/bin/env python3
"""Shared path adapter for bottom-layer retrieval visualizations."""

from pathlib import Path
import runpy
import sys


HERE = Path(__file__).resolve().parent
CAMPAIGN_ROOT = HERE.parent
RRS_ROOT = CAMPAIGN_ROOT.parent
RETRIEVAL_ROOT = CAMPAIGN_ROOT / "retrievals"
TRUTH_ROOT = CAMPAIGN_ROOT / "truth"
TRUTH_TABLE = TRUTH_ROOT / "true_states.dat"
SCENE_COMPONENTS = TRUTH_ROOT / "scene_components.dat"
INVERSION_ROOT = RRS_ROOT / "inversion"
VERTICAL_PROFILE_TABLE = (
    RRS_ROOT / "truth_map_aerosols" / "aerosol_vertical_profiles.dat"
)
ATMOSPHERIC_PROFILE = (
    INVERSION_ROOT / "retrieval_setup" / "retrieval_atmosphere_16layer.nc"
)
PLOT_ROOT = RETRIEVAL_ROOT / "plots"


def has_option(arguments, option):
    """Return true when an argparse-style option is already present."""
    return any(
        argument == option or argument.startswith(option + "=")
        for argument in arguments
    )


def run_inversion_script(script_name, defaults=(), subdirectory=None):
    """Run a shared inversion script with overridable bottom-layer defaults.

    Defaults are inserted before the user's arguments.  Python's argparse uses
    the last occurrence of a scalar option, so an explicit user value remains
    authoritative.  Avoiding copied plotting implementations keeps the
    full-column and bottom-layer diagnostics on one tested code path.
    """
    original_arguments = list(sys.argv[1:])
    injected = []
    for option, value in defaults:
        if not has_option(original_arguments, option):
            injected.extend([option, str(value)])
    sys.argv = [sys.argv[0]] + injected + original_arguments

    source_root = INVERSION_ROOT
    if subdirectory is not None:
        source_root = source_root / subdirectory
    source = source_root / script_name
    if not source.is_file():
        raise RuntimeError("shared plotting implementation is missing: %s" % source)
    sys.path.insert(0, str(INVERSION_ROOT))
    runpy.run_path(str(source), run_name="__main__")
