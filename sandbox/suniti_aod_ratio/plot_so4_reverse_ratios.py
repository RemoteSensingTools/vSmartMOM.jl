from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def load(kind, altitude=False):
    suffix = "_alt1p67" if altitude else ""
    return np.genfromtxt(
        ROOT / f"so4_{kind}{suffix}_spectra.csv", delimiter=",", names=True
    )


datasets = {
    "TOA": {kind: load(kind) for kind in ("original", "fixed")},
    "1.67 km above BOA": {
        kind: load(kind, altitude=True) for kind in ("original", "fixed")
    },
}
colors = {"original": "#c44e52", "fixed": "#2a6fbb"}


def reverse_ratio(data):
    return data["up_aod0"] / data["up_aod003"]


fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True,
                         constrained_layout=True)
for ax, (location, cases) in zip(axes, datasets.items()):
    reference_wn = None
    for kind, data in cases.items():
        wn = data["wavenumber_cm1"]
        if reference_wn is None:
            reference_wn = wn
        else:
            assert np.array_equal(reference_wn, wn)
        ratio = reverse_ratio(data)
        ax.plot(wn, ratio, lw=1.1, color=colors[kind], label=kind)
        finite = ratio[np.isfinite(ratio)]
        print(location, kind, "ratio quantiles", np.quantile(finite, [0, 0.5, 1]))

    original_ratio = reverse_ratio(cases["original"])
    fixed_ratio = reverse_ratio(cases["fixed"])
    delta = fixed_ratio - original_ratio
    finite = np.isfinite(delta)
    print(
        location,
        "fixed-original max abs", np.max(np.abs(delta[finite])),
        "RMS", np.sqrt(np.mean(delta[finite] ** 2)),
    )

    ax.axhline(1.0, color="0.25", lw=0.8, ls="--")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    ax.set_ylabel("upwelling radiance ratio\n(AOD 0 / AOD 0.03)")
    ax.set_title(f"SO$_4$ — {location}")

axes[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
combined_output = ROOT / "so4_aod0_over_aod003_upwelling_ratios.png"
fig.savefig(combined_output, dpi=180)
print(combined_output)


fig_alt, ax_alt = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
for kind, data in datasets["1.67 km above BOA"].items():
    ax_alt.plot(data["wavenumber_cm1"], reverse_ratio(data), lw=1.15,
                color=colors[kind], label=kind)
ax_alt.axhline(1.0, color="0.25", lw=0.8, ls="--")
ax_alt.grid(alpha=0.22)
ax_alt.legend(frameon=False)
ax_alt.set_ylabel("upwelling radiance ratio (AOD 0 / AOD 0.03)")
ax_alt.set_title("SO$_4$ upwelling at 1.67 km above BOA")
ax_alt.set_xlabel(r"Wavenumber (cm$^{-1}$)")
alt_output = ROOT / "so4_aod0_over_aod003_upwelling_ratios_alt1p67km.png"
fig_alt.savefig(alt_output, dpi=180)
print(alt_output)


# Two panels keep the altitude curves solid and give them an independent,
# tighter radiance-ratio scale.
fig_overlay, overlay_axes = plt.subplots(
    2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True
)
for ax, (location, cases) in zip(overlay_axes, datasets.items()):
    plotted = []
    for kind, data in cases.items():
        values = reverse_ratio(data)
        plotted.append(values[np.isfinite(values)])
        ax.plot(
            data["wavenumber_cm1"], values,
            color=colors[kind], lw=1.15, label=kind,
        )
    limits = np.concatenate(plotted)
    padding = 0.04 * (np.max(limits) - np.min(limits))
    ax.set_ylim(np.min(limits) - padding, np.max(limits) + padding)
    ax.axhline(1.0, color="0.25", lw=0.8, ls=":")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    ax.set_ylabel("upwelling ratio\n(AOD 0 / AOD 0.03)")
    ax.set_title(location)
overlay_axes[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
fig_overlay.suptitle("SO$_4$ upwelling — wet atmosphere")
overlay_output = ROOT / "so4_aod0_over_aod003_upwelling_TOA_and_alt1p67km.png"
fig_overlay.savefig(overlay_output, dpi=180)
print(overlay_output)
