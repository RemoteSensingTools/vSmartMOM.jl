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

fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True,
                         constrained_layout=True)
colors = {"original": "#c44e52", "fixed": "#2a6fbb"}

for ax, (location, cases) in zip(axes, datasets.items()):
    reference_wn = None
    for kind, data in cases.items():
        wn = data["wavenumber_cm1"]
        if reference_wn is None:
            reference_wn = wn
        else:
            assert np.array_equal(reference_wn, wn)
        ratio = data["up_aod003"] / data["up_aod0"]
        ax.plot(wn, ratio, lw=1.1, color=colors[kind], label=kind)
        finite = ratio[np.isfinite(ratio)]
        print(location, kind, "ratio quantiles", np.quantile(finite, [0, 0.5, 1]))

    original_ratio = cases["original"]["up_aod003"] / cases["original"]["up_aod0"]
    fixed_ratio = cases["fixed"]["up_aod003"] / cases["fixed"]["up_aod0"]
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
    ax.set_ylabel("upwelling radiance ratio\n(AOD 0.03 / AOD 0)")
    ax.set_title(f"SO$_4$ — {location}")

axes[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
output = ROOT / "so4_aod003_over_aod0_upwelling_ratios.png"
fig.savefig(output, dpi=180)
print(output)

# Direct counterpart of the earlier BC altitude-only figure.
fig_alt, ax_alt = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
alt_cases = datasets["1.67 km above BOA"]
for kind, data in alt_cases.items():
    ratio = data["up_aod003"] / data["up_aod0"]
    ax_alt.plot(data["wavenumber_cm1"], ratio, lw=1.15,
                color=colors[kind], label=kind)
ax_alt.axhline(1.0, color="0.25", lw=0.8, ls="--")
ax_alt.grid(alpha=0.22)
ax_alt.legend(frameon=False)
ax_alt.set_ylabel("upwelling radiance ratio (AOD 0.03 / AOD 0)")
ax_alt.set_title("SO$_4$ upwelling at 1.67 km above BOA")
ax_alt.set_xlabel(r"Wavenumber (cm$^{-1}$)")
alt_output = ROOT / "so4_aod003_over_aod0_upwelling_ratios_alt1p67km.png"
fig_alt.savefig(alt_output, dpi=180)
print(alt_output)
