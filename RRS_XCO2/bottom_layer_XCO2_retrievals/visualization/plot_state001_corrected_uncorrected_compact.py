#!/usr/bin/env python3
"""Compare terminal corrected/uncorrected bottom-layer states (default 001)."""

from _adapter import (
    RETRIEVAL_ROOT,
    SCENE_COMPONENTS,
    TRUTH_TABLE,
    run_inversion_script,
)


if __name__ == "__main__":
    run_inversion_script(
        "plot_state001_corrected_uncorrected_compact.py",
        (
            ("--inversion-root", RETRIEVAL_ROOT),
            ("--truth-table", TRUTH_TABLE),
            ("--scene-components", SCENE_COMPONENTS),
        ),
    )
