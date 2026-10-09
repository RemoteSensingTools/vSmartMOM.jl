#!/usr/bin/env python3
"""Export bottom-layer retrieval trajectories in physical units."""

from _adapter import TRUTH_TABLE, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "export_retrieval_state_trajectory.py",
        (("--truth-table", TRUTH_TABLE),),
    )
