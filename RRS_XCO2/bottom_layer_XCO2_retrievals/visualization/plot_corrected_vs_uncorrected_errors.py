#!/usr/bin/env python3
"""Plot corrected versus uncorrected bottom-layer retrieval displacements."""

from _adapter import (
    RETRIEVAL_ROOT,
    SCENE_COMPONENTS,
    TRUTH_TABLE,
    run_inversion_script,
)


if __name__ == "__main__":
    run_inversion_script(
        "plot_corrected_vs_uncorrected_errors.py",
        (
            ("--inversion-root", RETRIEVAL_ROOT),
            ("--truth-table", TRUTH_TABLE),
            ("--scene-components", SCENE_COMPONENTS),
        ),
    )
