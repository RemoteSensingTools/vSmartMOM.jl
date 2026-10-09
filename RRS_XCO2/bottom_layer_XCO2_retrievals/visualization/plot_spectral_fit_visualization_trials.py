#!/usr/bin/env python3
"""Create the bottom-layer spectral-fit visualization suite."""

from _adapter import PLOT_ROOT, RETRIEVAL_ROOT, TRUTH_TABLE, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "plot_spectral_fit_visualization_trials.py",
        (
            ("--inversion-dir", RETRIEVAL_ROOT),
            ("--truth-table", TRUTH_TABLE),
            ("--state", 13),
            ("--output-dir", PLOT_ROOT / "spectral_fit_visualizations"),
        ),
    )
