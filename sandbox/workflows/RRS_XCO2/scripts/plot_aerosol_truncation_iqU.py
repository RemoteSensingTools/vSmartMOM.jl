#!/usr/bin/env python3
"""Plot Stokes I/Q/U and truncation errors against the full-Mie reference."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "truth_map_aerosols" / "truncation_iqU.dat"
OUTPUT = ROOT / "truth_map_aerosols" / "truncation_vs_untruncated_IQU.png"
CAPS = (8, 12, 16, 24, 32)
COLORS = dict(zip(CAPS, plt.cm.viridis(np.linspace(0.05, 0.85, len(CAPS)))))


def load_rows():
    dtype = [("vza", float), ("vaz", float), ("case", "U12"),
             ("I", float), ("Q", float), ("U", float)]
    return np.loadtxt(DATA, dtype=dtype, comments="#")


def main():
    rows = load_rows()
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex="col")
    components = ("I", "Q")

    def signed_series(case, component):
        """RelAz=180° → negative VZA; RelAz=0° → positive VZA."""
        pos = rows[(rows["case"] == case) & (rows["vaz"] == 0.0)]
        neg = rows[(rows["case"] == case) & (rows["vaz"] == 180.0)]
        # Drop the duplicate negative-side nadir point.
        neg = neg[neg["vza"] > 0]
        x = np.concatenate((-neg["vza"], pos["vza"]))
        y = np.concatenate((neg[component], pos[component]))
        order = np.argsort(x)
        return x[order], y[order]

    for i, component in enumerate(components):
        ax, dax = axes[i]
        xref, ref = signed_series("untruncated", component)
        ax.plot(xref, ref, color="black", linewidth=2.4,
                label="untruncated")
        for cap in CAPS:
            xcur, cur = signed_series(str(cap), component)
            if not np.array_equal(xcur, xref):
                raise ValueError(f"Signed-VZA grid mismatch for l_trunc={cap}")
            ax.plot(xcur, cur, color=COLORS[cap], linewidth=1.25,
                    label=f"l_trunc={cap}")
            if component == "I":
                difference = 100.0 * (cur - ref) / ref
            else:
                difference = cur - ref
            dax.plot(xref, difference, color=COLORS[cap], linewidth=1.4,
                     label=f"l_trunc={cap}")

        ax.set_ylabel(component)
        if component == "I":
            dax.set_ylabel("100·(Iₗ−Iᵣᵉᶠ)/Iᵣᵉᶠ  [%]")
        else:
            dax.set_ylabel(f"{component}ₗ − {component}ᵣᵉᶠ  [absolute]")
        ax.grid(alpha=0.25)
        dax.grid(alpha=0.25)
        ax.axvline(0, color="0.65", linewidth=0.8)
        dax.axvline(0, color="0.65", linewidth=0.8)
        dax.axhline(0, color="0.35", linewidth=0.8)

    axes[0, 0].set_title("TOA Stokes radiance at 757.001655 nm")
    axes[0, 1].set_title("Difference from untruncated")
    axes[-1, 0].set_xlabel("Signed viewing zenith angle [deg]")
    axes[-1, 1].set_xlabel("Signed viewing zenith angle [deg]")
    for ax in axes[-1]:
        ax.set_xlim(-70, 70)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=9,
               bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(
        "Aerosol angular-truncation convergence; SZA=30°\n"
        "negative VZA: Δφ=180°; positive VZA: Δφ=0°", y=0.995)
    fig.tight_layout(rect=(0, 0.09, 1, 0.975))
    fig.savefig(OUTPUT, dpi=220)
    print(OUTPUT)


if __name__ == "__main__":
    main()
