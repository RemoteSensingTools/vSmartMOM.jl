#!/usr/bin/env python3
"""Round-3 versus round-6 four-panel pressure/XCO2 comparison.

Use --SIF=0, --SIF=1, or --SIF=both. Earlier-round outputs remain untouched.
"""

from plot_round3_vs_round4_xco2_psurf import main


if __name__ == "__main__":
    main(comparison_round=6)
