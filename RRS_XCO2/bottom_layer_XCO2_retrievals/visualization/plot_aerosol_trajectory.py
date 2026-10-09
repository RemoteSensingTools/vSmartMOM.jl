#!/usr/bin/env python3
"""Plot aerosol trial trajectories using bottom-layer truth metadata."""

from _adapter import SCENE_COMPONENTS, run_inversion_script


if __name__ == "__main__":
    run_inversion_script(
        "plot_aerosol_trajectory.py",
        (("--scene-components", SCENE_COMPONENTS),),
    )
