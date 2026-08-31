#!/usr/bin/env python3
"""Compare compact terminal states for matched state-001 retrievals.

The compact state replaces the twelve retrieved layer CO2 VMRs with the
derived dry-air-column XCO2 diagnostic.  Every other retrieval quantity is
retained.  Aerosol optical depths and profile-center heights are transformed
from their logarithmic retrieval coordinates back to physical values.
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


HERE = Path(__file__).resolve().parent
WORKFLOW_ROOT = HERE.parent
DATA_ROOT = Path(os.environ.get("RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
INVERSION_ROOT = DATA_ROOT / "inversion"
OUTPUT = INVERSION_ROOT / "state001_corrected_vs_uncorrected_compact_state.png"
TABLE = INVERSION_ROOT / "state001_corrected_vs_uncorrected_compact_state.dat"

PARAMETERS = (
    ("XCO2", "XCO$_2$ (ppm)"),
    ("psurf", "Surface pressure (hPa)"),
    ("sulfate_aod760", "Sulphate AOD$_{760}$"),
    ("organic_carbon_aod760", "Organic-carbon AOD$_{760}$"),
    ("utls_sulfate_aod760", "UTLS sulphate AOD$_{760}$"),
    ("sulfate_z0", "Sulphate height (km)"),
    ("organic_carbon_z0", "Organic-carbon height (km)"),
    ("utls_sulfate_z0", "UTLS sulphate height (km)"),
    ("o2a_surface_P0", "O$_2$ A surface $P_0$"),
    ("o2a_surface_P1", "O$_2$ A surface $P_1$"),
    ("o2a_surface_P2", "O$_2$ A surface $P_2$"),
    ("weak_co2_surface_P0", "Weak CO$_2$ surface $P_0$"),
    ("weak_co2_surface_P1", "Weak CO$_2$ surface $P_1$"),
    ("weak_co2_surface_P2", "Weak CO$_2$ surface $P_2$"),
    ("strong_co2_surface_P0", "Strong CO$_2$ surface $P_0$"),
    ("strong_co2_surface_P1", "Strong CO$_2$ surface $P_1$"),
    ("strong_co2_surface_P2", "Strong CO$_2$ surface $P_2$"),
    ("SIF760", "SIF760 (native units)"),
    ("mSIF", "mSIF (native units)"),
)

TRUTH = {
    "XCO2": 380.0,
    "psurf": 1000.0,
    "sulfate_aod760": 0.0,
    "organic_carbon_aod760": 0.0,
    "utls_sulfate_aod760": 0.0,
    # Nominal truth-map profile centers.  They are spectrally non-identifiable
    # for state 001 because all three corresponding truth AODs are zero.
    "sulfate_z0": 1.525651538,
    "organic_carbon_z0": 2.112319568,
    "utls_sulfate_z0": 12.120602005,
    "o2a_surface_P0": 2.715374708583e-1,
    "o2a_surface_P1": -1.200504636272e-3,
    "o2a_surface_P2": -2.299576457171e-4,
    "weak_co2_surface_P0": 2.522486169619e-1,
    "weak_co2_surface_P1": -2.142459711455e-3,
    "weak_co2_surface_P2": -1.317753254698e-4,
    "strong_co2_surface_P0": 2.176911056292e-1,
    "strong_co2_surface_P1": -1.859347500215e-3,
    "strong_co2_surface_P2": -1.138927333851e-4,
    "SIF760": 0.0,
    "mSIF": 0.0,
}


def physical_state(dataset, variable):
    names = str(dataset.getncattr("parameter_names")).split()
    native = dict(zip(names, np.asarray(dataset[variable][:], dtype=float)))
    state = {
        "psurf": native["psurf"],
        "SIF760": native["SIF760"],
        "mSIF": native["mSIF"],
    }
    for aerosol in ("sulfate", "organic_carbon", "utls_sulfate"):
        state[f"{aerosol}_aod760"] = np.exp(
            native[f"ln_{aerosol}_aod760"]
        )
        state[f"{aerosol}_z0"] = np.exp(native[f"ln_{aerosol}_z0"])
    for band in ("o2a", "weak_co2", "strong_co2"):
        for order in range(3):
            name = f"{band}_surface_P{order}"
            state[name] = native[name]
    return state


def load_class(retrieval_class):
    records = []
    for path in sorted(
        (INVERSION_ROOT / retrieval_class).glob("retrieval_state001_perturbation*.nc")
    ):
        with Dataset(path) as dataset:
            if int(dataset.getncattr("retrieval_complete")) != 1:
                continue
            if int(dataset.getncattr("truth_state_index")) != 1:
                continue
            perturbation = int(dataset.getncattr("perturbation_index"))
            final = physical_state(dataset, "final_state")
            final["XCO2"] = float(dataset["XCO2"][:])
            prior = physical_state(dataset, "a_priori_state")
            prior["XCO2"] = float(dataset["a_priori_XCO2"][:])
        records.append((perturbation, final, prior))
    if not records:
        raise RuntimeError(f"no completed state-001 {retrieval_class} retrievals")
    return records


def write_table(corrected, uncorrected):
    corrected_by_pert = {pert: state for pert, state, _ in corrected}
    uncorrected_by_pert = {pert: state for pert, state, _ in uncorrected}
    perturbations = sorted(set(corrected_by_pert) & set(uncorrected_by_pert))
    with TABLE.open("w", encoding="utf-8") as stream:
        stream.write(
            "# State 001 compact terminal-state comparison. Layer CO2 VMRs are "
            "replaced by XCO2.\n"
        )
        stream.write(
            "# Standard deviations are sample standard deviations across matched "
            "noise perturbations.\n"
        )
        stream.write(
            "# parameter truth corrected_mean corrected_sd uncorrected_mean "
            "uncorrected_sd paired_difference_mean paired_difference_sd\n"
        )
        for key, _ in PARAMETERS:
            corr = np.array([corrected_by_pert[p][key] for p in perturbations])
            uncorr = np.array([uncorrected_by_pert[p][key] for p in perturbations])
            delta = uncorr - corr
            stream.write(
                f"{key} {TRUTH[key]:.12e} "
                f"{corr.mean():.12e} {corr.std(ddof=1):.12e} "
                f"{uncorr.mean():.12e} {uncorr.std(ddof=1):.12e} "
                f"{delta.mean():.12e} {delta.std(ddof=1):.12e}\n"
            )


def main():
    corrected = load_class("corrected")
    uncorrected = load_class("uncorrected")
    corrected_by_pert = {pert: state for pert, state, _ in corrected}
    uncorrected_by_pert = {pert: state for pert, state, _ in uncorrected}
    perturbations = sorted(set(corrected_by_pert) & set(uncorrected_by_pert))
    if len(perturbations) != 10:
        raise RuntimeError(
            f"expected ten matched perturbations, found {len(perturbations)}"
        )
    prior = corrected[0][2]

    fig, axes = plt.subplots(5, 4, figsize=(18, 20))
    axes = axes.ravel()
    colors = {"corrected": "#1f77b4", "uncorrected": "#d62728"}
    markers = {"corrected": "o", "uncorrected": "s"}

    for axis, (key, label) in zip(axes, PARAMETERS):
        for retrieval_class, states in (
            ("corrected", corrected_by_pert),
            ("uncorrected", uncorrected_by_pert),
        ):
            values = np.array([states[p][key] for p in perturbations])
            axis.plot(
                perturbations,
                values,
                color=colors[retrieval_class],
                marker=markers[retrieval_class],
                markersize=4.5,
                linewidth=1.4,
                label=retrieval_class.capitalize(),
            )
            mean = values.mean()
            sigma = values.std(ddof=1)
            axis.axhspan(
                mean - sigma,
                mean + sigma,
                color=colors[retrieval_class],
                alpha=0.07,
                linewidth=0,
            )

        if np.isfinite(TRUTH[key]):
            axis.axhline(TRUTH[key], color="black", linewidth=1.4, label="Truth")
        axis.axhline(
            prior[key], color="0.45", linestyle=":", linewidth=1.5, label="Prior"
        )
        axis.set_title(label, fontsize=11)
        axis.set_xticks((1, 3, 5, 7, 9))
        axis.grid(alpha=0.22)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))

    for axis in axes[:16]:
        axis.set_xlabel("Noise perturbation")
    axes[16].set_xlabel("Noise perturbation")
    axes[17].set_xlabel("Noise perturbation")
    axes[18].set_xlabel("Noise perturbation")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].axis("off")
    axes[-1].legend(handles, labels, frameon=False, loc="center", fontsize=12)
    axes[-1].text(
        0.5,
        0.24,
        "Shading: ensemble mean $\\pm1\\sigma$\n"
        "State 001 aerosol heights are nominal profile centers;\n"
        "they are non-identifiable because truth AOD = 0",
        ha="center",
        va="center",
        fontsize=10,
        transform=axes[-1].transAxes,
    )
    fig.suptitle(
        "State 001: corrected versus uncorrected terminal retrieval states\n"
        "Layer CO$_2$ VMRs replaced by XCO$_2$",
        fontsize=16,
        y=0.985,
    )
    fig.subplots_adjust(
        left=0.06, right=0.985, bottom=0.045, top=0.935, hspace=0.52, wspace=0.30
    )
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)
    write_table(corrected, uncorrected)
    print(OUTPUT)
    print(TABLE)


if __name__ == "__main__":
    main()
