#!/usr/bin/env python3
"""Plot bottom-layer state evolution and OE convergence diagnostics."""

from _adapter import SCENE_COMPONENTS, TRUTH_TABLE, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "plot_retrieval_state_convergence.py",
        (
            ("--truth-table", TRUTH_TABLE),
            ("--scene-components", SCENE_COMPONENTS),
        ),
    )
