#!/usr/bin/env python3
"""Plot the forest-scene RRS discontinuity before and after the surface fix."""

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_rrs_jump_before_after")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OBSOLETE_DIAGNOSTICS = (
    ROOT / "obsolete" / "truth_map" / "chunk_surface_error_diagnostics"
)
DATA = OBSOLETE_DIAGNOSTICS / "rrs_chunk42_43_before_after.dat"
OUTPUT = OBSOLETE_DIAGNOSTICS / "rrs_chunk42_43_before_after.png"


def main():
    OBSOLETE_DIAGNOSTICS.mkdir(parents=True, exist_ok=True)
    source_wavelength, _, before_source, after_source = np.loadtxt(DATA, unpack=True)
    order = np.argsort(source_wavelength)
    wavelength = source_wavelength[order]
    before = before_source[order]
    after = after_source[order]

    # The source table is in increasing wavenumber. The boundary is therefore
    # between source rows 64 and 65, irrespective of wavelength plot ordering.
    boundary = 0.5 * (source_wavelength[63] + source_wavelength[64])

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), sharex=True)
    cases = (
        (before, before_source, "Before: chunk-local surface normalization", "#d62728"),
        (after, after_source, "After: full-band surface normalization", "#1f77b4"),
    )
    for ax, (values, source_values, title, color) in zip(axes, cases):
        ax.plot(wavelength, values, color=color, linewidth=1.5)
        ax.axvline(boundary, color="black", linestyle="--", linewidth=1.0,
                   label="chunk 42/43 boundary")
        ax.scatter(source_wavelength[63:65], source_values[63:65],
                   color="black", s=22, zorder=3)
        ax.set_title(title, fontsize=11, pad=10)
        ax.set_xlabel("Wavelength (nm)")
        ax.grid(alpha=0.22)
        ax.legend(loc="best")

    axes[0].set_ylabel(r"RRS Stokes I [mW m$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$]")
    # Use common limits so the disappearance of the artificial jump is visual,
    # while allowing for the small nstreams=9 versus nstreams=8 baseline shift.
    ymin = min(before.min(), after.min())
    ymax = max(before.max(), after.max())
    pad = 0.06 * (ymax - ymin)
    for ax in axes:
        ax.set_ylim(ymin - pad, ymax + pad)

    old_jump = before_source[64] - before_source[63]
    new_jump = after_source[64] - after_source[63]
    fig.suptitle(
        "Forest aerosol scene 57: O₂ A-band RRS chunk-boundary correction\n"
        f"adjacent-point jump: {old_jump:+.3e} → {new_jump:+.3e} "
        f"({abs(old_jump/new_jump):.0f}× smaller)",
        fontsize=13, y=0.97,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    fig.savefig(OUTPUT, dpi=180, bbox_inches="tight")
    print(OUTPUT)


if __name__ == "__main__":
    main()
