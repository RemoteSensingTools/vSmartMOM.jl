#!/usr/bin/env python3
"""Compare all 80 bottom-layer synthetic radiances with OCO-2 observations."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _adapter import PLOT_ROOT, TRUTH_ROOT, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "compare_oco2_observed_radiances.py",
        (
            ("--synthetic-root", TRUTH_ROOT / "OCO_radiances"),
            (
                "--output-root",
                PLOT_ROOT / "validation_against_OCO2",
            ),
            ("--expected-scene-count", 80),
        ),
        subdirectory="instrument",
    )
