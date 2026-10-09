#!/usr/bin/env python3
"""Fit representative OCO-2 EOFs to bottom-layer retrieval residuals."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _adapter import PLOT_ROOT, RETRIEVAL_ROOT, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "fit_oco2_eofs_to_retrieval_residuals.py",
        (
            ("--inversion-root", RETRIEVAL_ROOT),
            ("--output-root", PLOT_ROOT / "eof_residual_fits"),
        ),
        subdirectory="instrument",
    )
