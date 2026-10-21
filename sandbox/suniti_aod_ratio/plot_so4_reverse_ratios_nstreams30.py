from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
COLORS = {"original": "#c44e52", "fixed": "#2a6fbb"}


def load(kind, altitude=False, dry=False, nstreams30=True):
    altitude_suffix = "_alt1p67" if altitude else ""
    dry_suffix = "_dry" if dry else ""
    streams_suffix = "_nstreams30" if nstreams30 else ""
    path = ROOT / (
        f"so4_{kind}{altitude_suffix}{dry_suffix}{streams_suffix}_spectra.csv"
    )
    return np.genfromtxt(path, delimiter=",", names=True)


def ratio(data):
    return data["up_aod0"] / data["up_aod003"]


def make_plot(dry):
    datasets = {
        "TOA": {
            kind: load(kind, dry=dry) for kind in ("original", "fixed")
        },
        "1.67 km above BOA": {
            kind: load(kind, altitude=True, dry=dry)
            for kind in ("original", "fixed")
        },
    }
    fig, axes = plt.subplots(
        2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True
    )
    for ax, (location, cases) in zip(axes, datasets.items()):
        plotted = []
        for kind, data in cases.items():
            values = ratio(data)
            finite = values[np.isfinite(values)]
            plotted.append(finite)
            ax.plot(
                data["wavenumber_cm1"], values,
                color=COLORS[kind], lw=1.15, label=kind,
            )
            baseline = ratio(load(
                kind, altitude=(location != "TOA"), dry=dry, nstreams30=False
            ))
            difference = values - baseline
            print(
                "dry" if dry else "wet", location, kind,
                "ratio min median mean max", *map(float, (
                    np.min(finite), np.median(finite),
                    np.mean(finite), np.max(finite),
                )),
            )
            print(
                "  30-stream minus 8-stream: median_abs max_abs rms", *map(float, (
                    np.median(np.abs(difference)),
                    np.max(np.abs(difference)),
                    np.sqrt(np.mean(difference**2)),
                )),
            )
        limits = np.concatenate(plotted)
        padding = 0.04 * (np.max(limits) - np.min(limits))
        ax.set_ylim(np.min(limits) - padding, np.max(limits) + padding)
        ax.axhline(1.0, color="0.25", lw=0.8, ls=":")
        ax.grid(alpha=0.22)
        ax.legend(frameon=False)
        ax.set_ylabel("upwelling ratio\n(AOD 0 / AOD 0.03)")
        ax.set_title(location)

    axes[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
    atmosphere = "dry atmosphere (q = 0)" if dry else "wet atmosphere"
    fig.suptitle(
        f"SO$_4$ upwelling — VZA = 79.1° — 30 streams — {atmosphere}"
    )
    dry_suffix = "_dry" if dry else ""
    output = ROOT / (
        f"so4{dry_suffix}_nstreams30_aod0_over_aod003_"
        "upwelling_TOA_and_alt1p67km.png"
    )
    fig.savefig(output, dpi=180)
    print(output)


for dry_atmosphere in (False, True):
    make_plot(dry_atmosphere)
