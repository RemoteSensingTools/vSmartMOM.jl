from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "so4_fixed_vza0to80_alt1p67_spectra.csv"
OUTPUT = ROOT / "so4_fixed_vza0to80_aod0_over_aod003_TOA_and_alt1p67km.png"
OUTPUT_LOG = ROOT / "so4_fixed_vza0to80_aod0_over_aod003_TOA_and_alt1p67km_log.png"
OUTPUT_LOGDIFF = ROOT / "so4_fixed_vza0to80_logI_aod0_minus_logI_aod003_TOA_and_alt1p67km.png"
OUTPUT_LOGDIFF_ZOOM = ROOT / "so4_fixed_logI_difference_with_vza0to60_zoom.png"

data = np.genfromtxt(DATA, delimiter=",", names=True)
vzas = np.unique(data["vza_deg"])
cmap = plt.get_cmap("viridis")
norm = Normalize(vmin=float(vzas.min()), vmax=float(vzas.max()))

fig, axes = plt.subplots(
    2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True
)
fields = (
    ("toa_aod0", "toa_aod003", "TOA"),
    ("alt1p67_aod0", "alt1p67_aod003", "1.67 km above BOA"),
)

for ax, (numerator, denominator, title) in zip(axes, fields):
    all_ratios = []
    print(title)
    for vza in vzas:
        selected = data["vza_deg"] == vza
        subset = data[selected]
        ratio = subset[numerator] / subset[denominator]
        all_ratios.append(ratio[np.isfinite(ratio)])
        ax.plot(
            subset["wavenumber_cm1"], ratio,
            color=cmap(norm(vza)), lw=0.95,
        )
        print(
            f"  VZA={vza:4.0f}: min={np.min(ratio):.6g} "
            f"median={np.median(ratio):.6g} max={np.max(ratio):.6g}"
        )
    limits = np.concatenate(all_ratios)
    padding = 0.035 * (np.max(limits) - np.min(limits))
    ax.set_ylim(np.min(limits) - padding, np.max(limits) + padding)
    ax.axhline(1.0, color="0.25", lw=0.8, ls=":")
    ax.grid(alpha=0.20)
    ax.set_ylabel("upwelling ratio\n(AOD 0 / AOD 0.03)")
    ax.set_title(title)

axes[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
fig.suptitle(
    "SO$_4$ fixed optics — 8 streams — wet atmosphere — VZA sweep"
)
scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
scalar_mappable.set_array([])
colorbar = fig.colorbar(scalar_mappable, ax=axes, pad=0.015)
colorbar.set_label("Viewing zenith angle (degrees)")
colorbar.set_ticks(vzas[::2])
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)

for ax in axes:
    positive = np.concatenate([
        line.get_ydata()[line.get_ydata() > 0] for line in ax.lines[:-1]
    ])
    ax.set_yscale("log")
    ax.set_ylim(np.min(positive) / 1.25, np.max(positive) * 1.04)
fig.suptitle(
    "SO$_4$ fixed optics — 8 streams — wet atmosphere — VZA sweep (log ratio)"
)
fig.savefig(OUTPUT_LOG, dpi=180)
print(OUTPUT_LOG)

fig_logdiff, axes_logdiff = plt.subplots(
    2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True
)
for ax, (numerator, denominator, title) in zip(axes_logdiff, fields):
    for vza in vzas:
        selected = data["vza_deg"] == vza
        subset = data[selected]
        log_difference = np.log(subset[numerator]) - np.log(subset[denominator])
        ax.plot(
            subset["wavenumber_cm1"], log_difference,
            color=cmap(norm(vza)), lw=0.95,
        )
    ax.axhline(0.0, color="0.25", lw=0.8, ls=":")
    ax.grid(alpha=0.20)
    ax.set_ylabel(r"$\ln I_{\rm AOD=0}-\ln I_{\rm AOD=0.03}$")
    ax.set_title(title)

axes_logdiff[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
fig_logdiff.suptitle(
    "SO$_4$ fixed optics — 8 streams — wet atmosphere — VZA sweep"
)
scalar_mappable_logdiff = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
scalar_mappable_logdiff.set_array([])
colorbar_logdiff = fig_logdiff.colorbar(
    scalar_mappable_logdiff, ax=axes_logdiff, pad=0.015
)
colorbar_logdiff.set_label("Viewing zenith angle (degrees)")
colorbar_logdiff.set_ticks(vzas[::2])
fig_logdiff.savefig(OUTPUT_LOGDIFF, dpi=180)
print(OUTPUT_LOGDIFF)

fig_zoom, axes_zoom = plt.subplots(
    2, 2, figsize=(14, 8), sharex=True, constrained_layout=True
)
for row, (numerator, denominator, level_title) in enumerate(fields):
    for column, (vza_limit, range_title) in enumerate(
        ((80.0, "full VZA range"), (60.0, "VZA 0–60° zoom"))
    ):
        ax = axes_zoom[row, column]
        for vza in vzas[vzas <= vza_limit]:
            selected = data["vza_deg"] == vza
            subset = data[selected]
            log_difference = np.log(subset[numerator]) - np.log(subset[denominator])
            ax.plot(
                subset["wavenumber_cm1"], log_difference,
                color=cmap(norm(vza)), lw=0.95,
            )
        ax.axhline(0.0, color="0.25", lw=0.8, ls=":")
        ax.grid(alpha=0.20)
        ax.set_title(f"{level_title} — {range_title}")
        if column == 0:
            ax.set_ylabel(r"$\ln I_{\rm AOD=0}-\ln I_{\rm AOD=0.03}$")
        if row == 1:
            ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")

fig_zoom.suptitle(
    "SO$_4$ fixed optics — 8 streams — wet atmosphere — VZA sweep"
)
scalar_mappable_zoom = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
scalar_mappable_zoom.set_array([])
colorbar_zoom = fig_zoom.colorbar(
    scalar_mappable_zoom, ax=axes_zoom, pad=0.012
)
colorbar_zoom.set_label("Viewing zenith angle (degrees)")
colorbar_zoom.set_ticks(vzas[::2])
fig_zoom.savefig(OUTPUT_LOGDIFF_ZOOM, dpi=180)
print(OUTPUT_LOGDIFF_ZOOM)
