#!/usr/bin/env python3
"""Compare δBGE runs with the Float64 50-stream, m=0:99 reference."""

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_truncation_50stream")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "truth_map_aerosols"
TRUNCATED = ROOT / "truncation_iqU.dat"
REFERENCE = ROOT / "untruncated_stokes_by_m_float64_nstreams50_mmax99.dat"
OUTPUT = ROOT / "truncation_vs_untruncated_50stream_IQ.png"
CAPS = (8, 12, 16, 24, 32)
COLORS = dict(zip(CAPS, plt.cm.viridis(np.linspace(0.05, 0.85, len(CAPS)))))


def load_truncated():
    dtype = [("vza", float), ("vaz", float), ("case", "U12"),
             ("I", float), ("Q", float), ("U", float)]
    return np.loadtxt(TRUNCATED, dtype=dtype, comments="#")


def signed_truncated(rows, case, component):
    """Map relative azimuth 180° to negative VZA and 0° to positive VZA."""
    pos = rows[(rows["case"] == case) & (rows["vaz"] == 0.0)]
    neg = rows[(rows["case"] == case) & (rows["vaz"] == 180.0)]
    neg = neg[neg["vza"] > 0]
    x = np.concatenate((-neg["vza"], pos["vza"]))
    y = np.concatenate((neg[component], pos[component]))
    order = np.argsort(x)
    return x[order], y[order]


def load_reference():
    data = np.loadtxt(REFERENCE, comments="#")
    orders = np.unique(data[:, 0]).astype(int)
    final = data[data[:, 0] == orders[-1]]
    # Columns 5 and 6 are cumulative I and Q, already on the signed-VZA grid.
    order = np.argsort(final[:, 1])
    return final[order, 1], final[order, 5:7], orders[-1]


def main():
    rows = load_truncated()
    xref, reference, m_max = load_reference()
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex="col")

    for ic, component in enumerate(("I", "Q")):
        ax, dax = axes[ic]
        ref = reference[:, ic]
        ax.plot(xref, ref, color="black", linewidth=2.5,
                label=f"untruncated Float64, 50 streams, m≤{m_max}")

        for cap in CAPS:
            xcur, cur = signed_truncated(rows, str(cap), component)
            if not np.array_equal(xcur, xref):
                raise ValueError(f"Signed-VZA grid mismatch for l_trunc={cap}")
            ax.plot(xcur, cur, color=COLORS[cap], linewidth=1.25,
                    label=f"l_trunc={cap}")
            difference = (100.0 * (cur - ref) / ref
                          if component == "I" else cur - ref)
            dax.plot(xref, difference, color=COLORS[cap], linewidth=1.4)

        ax.set_ylabel(component)
        dax.set_ylabel("100·(Iₗ−Iᵣᵉᶠ)/Iᵣᵉᶠ  [%]" if component == "I"
                       else "Qₗ − Qᵣᵉᶠ  [absolute]")
        for panel in (ax, dax):
            panel.grid(alpha=0.25)
            panel.axvline(0, color="0.65", linewidth=0.8)
        dax.axhline(0, color="0.35", linewidth=0.8)

    axes[0, 0].set_title("TOA Stokes radiance at 757.001655 nm")
    axes[0, 1].set_title("Difference from 50-stream untruncated reference")
    for ax in axes[-1]:
        ax.set_xlabel("Signed viewing zenith angle [deg]")
        ax.set_xlim(-70, 70)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9,
               bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(
        "Aerosol angular-truncation convergence; SZA=30°\n"
        "Float64 50-stream untruncated reference, m=0…99; "
        "negative VZA: Δφ=180°, positive VZA: Δφ=0°",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.1, 1, 0.965))
    fig.savefig(OUTPUT, dpi=220)
    print(OUTPUT)


if __name__ == "__main__":
    main()
