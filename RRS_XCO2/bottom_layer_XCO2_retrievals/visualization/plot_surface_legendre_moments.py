#!/usr/bin/env python3
"""Plot bottom-layer truth, prior, and retrieved surface moments."""

from _adapter import TRUTH_TABLE, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "plot_surface_legendre_moments.py",
        (("--truth-table", TRUTH_TABLE),),
    )
