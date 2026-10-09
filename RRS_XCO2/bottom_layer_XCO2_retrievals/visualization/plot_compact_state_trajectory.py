#!/usr/bin/env python3
"""Compare bottom-layer compact state trajectories for several retrievals."""

from _adapter import SCENE_COMPONENTS, TRUTH_TABLE, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "plot_compact_state_trajectory.py",
        (
            ("--truth-table", TRUTH_TABLE),
            ("--scene-components", SCENE_COMPONENTS),
        ),
    )
