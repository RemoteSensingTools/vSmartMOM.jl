#!/usr/bin/env python3
"""Compare round-3 and round-4/5/6 XCO2 and pressure truth errors.

Each point is one matched, complete corrected/uncorrected retrieval pair.  The
horizontal coordinate is corrected retrieval minus truth and the vertical
coordinate is uncorrected retrieval minus truth.  Perturbations 01--10 are
circles; noiseless perturbation 11 is a square.  Color identifies surface type;
filled markers have no aerosol and open markers have AOD760 = 0.28.  No
complete pair is removed when its uncorrected retrieval fails to converge;
such pairs receive an x overlay.

Use --SIF both to pool no-SIF and SIF-on retrievals in the same four panels.
The circle/square and open/filled conventions remain unchanged; a small plus
inside a marker identifies SIF-on. Combined summaries use all 880 pairs per
round, including the noiseless members, and the DAT retains the SIF category.
"""

import argparse
from datetime import datetime
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


HERE = Path(__file__).resolve().parent
RRS_ROOT = HERE.parent
INVERSION_ROOT = RRS_ROOT / "inversion"
BOTTOM_ROOT = RRS_ROOT / "bottom_layer_XCO2_retrievals"
TRUTH_TABLE = BOTTOM_ROOT / "truth" / "true_states.dat"
SCENE_COMPONENTS = BOTTOM_ROOT / "truth" / "scene_components.dat"
DEFAULT_OUTPUT_ROOT = (
    BOTTOM_ROOT / "comparisons" / "round3_vs_round4_nosif"
)
PRIVATE_RESULTS_ROOT = Path(os.environ.get(
    "RRS_XCO2_PRIVATE_RESULTS_ROOT",
    str(Path.home() / "RRS_XCO2_private" / "results"),
))

sys.path.insert(0, str(INVERSION_ROOT))
from plot_all_corrected_vs_uncorrected_errors import (  # noqa: E402
    candidate_states,
    collect_pairs,
)


NO_SIF_CAMPAIGNS = {
    "Round 3": (
        BOTTOM_ROOT /
        "retrievals_acos_mapped_tapered_vertical_correlation_nosif",
        TRUTH_TABLE,
        ("legacy",),
    ),
    "Round 4": (
        BOTTOM_ROOT / "round4_known_sif759" / "retrievals_nosif",
        TRUTH_TABLE,
        ("off",),
    ),
}
ROUND3_SIF_ROOT = (
    PRIVATE_RESULTS_ROOT /
    "bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1"
)
SIF_ON_CAMPAIGNS = {
    "Round 3": (
        ROUND3_SIF_ROOT / "retrievals",
        ROUND3_SIF_ROOT / "retrieval_setup" /
        "true_states_corrected_sif_v2.dat",
        ("legacy",),
    ),
    "Round 4": (
        PRIVATE_RESULTS_ROOT /
        "bottom_layer_round4_known_sif759_sif_on_"
        "acos_mapped_tapered_vertical_correlation_v1" / "retrievals",
        ROUND3_SIF_ROOT / "retrieval_setup" /
        "true_states_corrected_sif_v2.dat",
        ("on",),
    ),
}
ROUND5_NO_SIF_CAMPAIGNS = {
    "Round 3": NO_SIF_CAMPAIGNS["Round 3"],
    "Round 5": (
        BOTTOM_ROOT / "round5_fixed_sif" / "retrievals_nosif",
        TRUTH_TABLE,
        ("off",),
    ),
}
ROUND5_SIF_ON_CAMPAIGNS = {
    "Round 3": SIF_ON_CAMPAIGNS["Round 3"],
    "Round 5": (
        PRIVATE_RESULTS_ROOT /
        "bottom_layer_round5_fixed_sif_on_tight_utls_"
        "acos_mapped_tapered_vertical_correlation_v1" / "retrievals",
        SIF_ON_CAMPAIGNS["Round 3"][1],
        ("on",),
    ),
}
EXPECTED_PAIRS = 440
ROUND6_NO_SIF_CAMPAIGNS = {
    "Round 3": NO_SIF_CAMPAIGNS["Round 3"],
    "Round 6": (
        BOTTOM_ROOT / "round6_fixed_sif" / "retrievals_nosif",
        TRUTH_TABLE,
        ("off",),
    ),
}
ROUND6_SIF_ON_CAMPAIGNS = {
    "Round 3": SIF_ON_CAMPAIGNS["Round 3"],
    "Round 6": (
        PRIVATE_RESULTS_ROOT /
        "bottom_layer_round6_fixed_sif_on_standard_utls_"
        "acos_mapped_tapered_vertical_correlation_v1" / "retrievals",
        SIF_ON_CAMPAIGNS["Round 3"][1],
        ("on",),
    ),
}
COMPARISON_CAMPAIGNS = {
    4: (NO_SIF_CAMPAIGNS, SIF_ON_CAMPAIGNS),
    5: (ROUND5_NO_SIF_CAMPAIGNS, ROUND5_SIF_ON_CAMPAIGNS),
    6: (ROUND6_NO_SIF_CAMPAIGNS, ROUND6_SIF_ON_CAMPAIGNS),
}
METRICS = (
    ("psurf", "Surface-pressure error", "hPa"),
    ("XCO2", "Column-averaged XCO$_2$ error", "ppm"),
)
SURFACE_COLORS = {
    "urban": "#0072B2",
    "rural": "#CC79A7",
    "desert": "#D55E00",
    "forest": "#009E73",
}


def parse_args(comparison_round=4):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--compare-round", type=int, choices=tuple(COMPARISON_CAMPAIGNS), default=comparison_round,
        help="Compare round 3 with this round (default: %(default)s)",
    )
    parser.add_argument(
        "--SIF", dest="sif_enabled", choices=("0", "1", "both"), default="0",
        help="Select no-SIF (0, default), SIF-on (1), or pool both categories",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output PNG (default depends on --SIF)",
    )
    parser.add_argument(
        "--table-output", type=Path, default=None,
        help="Whitespace-separated table (default depends on --SIF)",
    )
    parser.add_argument(
        "--allow-partial-comparison", "--allow-partial-round4",
        "--allow-partial-round5", "--allow-partial-round6", dest="allow_partial_comparison",
        action="store_true",
        help=(
            "Plot currently available matched comparison-round SIF pairs before the "
            "440-pair campaign is complete; otherwise retain placeholders"
        ),
    )
    parser.add_argument("--dpi", type=int, default=200)
    return parser.parse_args()


def load_campaign(label, root, truth_table, accepted_sif_modes,
                  allow_incomplete=False, expected_state_model=None):
    if not root.is_dir():
        return [], 0
    states = candidate_states(root)
    pairs, unmatched = collect_pairs(
        root, states, truth_table, SCENE_COMPONENTS
    )
    if unmatched and not allow_incomplete:
        details = "; ".join(
            "state %03d corrected-only=%s uncorrected-only=%s" % (
                state, corrected, uncorrected,
            )
            for state, corrected, uncorrected in unmatched
        )
        raise RuntimeError("%s has unmatched complete products: %s" % (
            label, details,
        ))
    unexpected = sorted(set(
        pair["metadata"]["sif_mode"] for pair in pairs
        if pair["metadata"]["sif_mode"] not in accepted_sif_modes
    ))
    if unexpected:
        raise RuntimeError(
            "%s has unexpected SIF metadata: %s" % (label, unexpected)
        )
    if expected_state_model is not None and any(
            pair["metadata"]["state_model"] != expected_state_model
            for pair in pairs):
        raise RuntimeError("%s requires state model %s" % (label, expected_state_model))
    return pairs, len(pairs)


def errors(pairs, key):
    corrected = np.asarray([
        pair["corrected"]["final"][key] - pair["truth"][key]
        for pair in pairs
    ])
    uncorrected = np.asarray([
        pair["uncorrected"]["final"][key] - pair["truth"][key]
        for pair in pairs
    ])
    return corrected, uncorrected


def shared_limit(campaign_pairs, key):
    values = [np.asarray([0.0])]
    for pairs in campaign_pairs.values():
        values.extend(errors(pairs, key))
    maximum = float(np.max(np.abs(np.concatenate(values))))
    return 1.12 * maximum if maximum > 0.0 else 1.0


def write_table(path, campaign_pairs, snapshot, sif_enabled, combined=False,
                comparison_round=4):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# Round-3 versus round-%d %s truth-displacement comparison.\n"
            % (comparison_round, "combined no-SIF and SIF-on" if combined else
               "SIF-on" if sif_enabled else "no-SIF")
        )
        stream.write("# Snapshot %s\n" % snapshot)
        stream.write(
            "# Panel summaries use the arithmetic mean and sample standard "
            "deviation (N-1) over every plotted matched pair, including "
            "noiseless perturbation 11 and failed retrievals.\n"
        )
        stream.write(
            "# round state perturbation noise_case surface aerosol_case "
            "truth_bottom_co2_ppm truth_xco2_ppm parameter units "
            "corrected_minus_truth uncorrected_minus_truth "
            "corrected_converged uncorrected_converged corrected_file "
            "uncorrected_file" + (" sif_category" if combined else "") + "\n"
        )
        for round_label, pairs in campaign_pairs.items():
            for pair in pairs:
                metadata = pair["metadata"]
                for key, _, unit in METRICS:
                    corrected, uncorrected = errors([pair], key)
                    row = (
                        "%s %03d %02d %s %s %s %.12g %.12g %s %s "
                        "%.12g %.12g %d %d %s %s" % (
                            round_label.replace(" ", "_"),
                            pair["state"], pair["perturbation"],
                            pair["noise_case"], metadata["surface"],
                            metadata["aerosol_case"],
                            metadata["truth_bottom_co2"],
                            metadata["truth_xco2"], key, unit,
                            corrected[0], uncorrected[0],
                            int(pair["corrected"]["converged"]),
                            int(pair["uncorrected"]["converged"]),
                            pair["corrected"]["path"],
                            pair["uncorrected"]["path"],
                        )
                    )
                    if combined:
                        row += " " + ("on" if pair["sif_enabled"] else "off")
                    stream.write(row + "\n")


def scatter_group(axis, pairs, key, color, filled, noiseless, combined):
    """Retain circle/square and aerosol fills; mark SIF-on with a central +."""
    if not pairs:
        return
    x, y = errors(pairs, key)
    marker = "s" if noiseless else "o"
    size = (34 if noiseless else 24) if combined else (42 if noiseless else 24)
    alpha = (0.66 if combined else 0.72) if noiseless else (
        (0.36 if filled else 0.62) if combined else
        (0.43 if filled else 0.68)
    )
    linewidth = (0.7 if filled else 1.2) if noiseless else (
        0.4 if filled else 0.9
    )
    zorder = 5 if noiseless else 3
    axis.scatter(
        x, y, marker=marker, s=size,
        facecolors=color if filled else "none", edgecolors=color,
        alpha=alpha, linewidths=linewidth, zorder=zorder,
    )
    if combined:
        sif_on = np.asarray([pair["sif_enabled"] for pair in pairs], dtype=bool)
        axis.scatter(
            x[sif_on], y[sif_on], marker="+", s=11 if noiseless else 8,
            color="0.20", alpha=0.72, linewidths=0.6, zorder=zorder + 0.1,
        )


def draw_panel(axis, pairs, key, title, unit, limit, combined=False):
    axis.plot(
        [-limit, limit], [-limit, limit], color="0.35", linestyle="--",
        linewidth=1.0, zorder=1,
    )
    axis.axhline(0.0, color="0.72", linewidth=0.8, zorder=0)
    axis.axvline(0.0, color="0.72", linewidth=0.8, zorder=0)

    for surface, color in SURFACE_COLORS.items():
        for aerosol_case in ("none", "aod760_0p28"):
            selected = [
                pair for pair in pairs
                if pair["metadata"]["surface"] == surface and
                pair["metadata"]["aerosol_case"] == aerosol_case
            ]
            noisy = [
                pair for pair in selected if pair["perturbation"] != 11
            ]
            noiseless = [
                pair for pair in selected if pair["perturbation"] == 11
            ]
            filled = aerosol_case == "none"
            scatter_group(axis, noisy, key, color, filled, False, combined)
            scatter_group(axis, noiseless, key, color, filled, True, combined)

    uncorrected_failed = [
        pair for pair in pairs
        if not pair["uncorrected"]["converged"]
    ]
    if uncorrected_failed:
        x, y = errors(uncorrected_failed, key)
        axis.scatter(
            x, y, marker="x", s=56, color="#B2182B", linewidths=1.3,
            zorder=7,
        )

    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_aspect("equal", adjustable="box")
    axis.set_title("%s (%s)" % (title, unit), fontsize=12, pad=7)
    axis.tick_params(labelsize=9)
    axis.grid(alpha=0.16)
    corrected_error, uncorrected_error = errors(pairs, key)
    corrected_std = float(np.std(corrected_error, ddof=1)) \
        if len(corrected_error) > 1 else 0.0
    uncorrected_std = float(np.std(uncorrected_error, ddof=1)) \
        if len(uncorrected_error) > 1 else 0.0
    bias_symbol = (
        r"$\overline{\Delta XCO_2}$" if key == "XCO2" else
        r"$\overline{\Delta p_{\mathrm{surf}}}$"
    )
    summary = (
        "n=%d; %s [%s]\n"
        "corr. = %+.3f ± %.3f\n"
        "uncorr. = %+.3f ± %.3f" % (
            len(pairs), bias_symbol, unit,
            float(np.mean(corrected_error)), corrected_std,
            float(np.mean(uncorrected_error)), uncorrected_std,
        )
    )
    axis.text(
        0.03, 0.97, summary,
        transform=axis.transAxes, ha="left", va="top", fontsize=9,
        bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.82,
                  boxstyle="round,pad=0.25"),
    )


def draw_placeholder(axis, title, unit, available_pairs, expected_pairs=EXPECTED_PAIRS,
                     comparison_round=4):
    axis.set_title("%s (%s)" % (title, unit), fontsize=12, pad=7)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_color("0.82")
    axis.text(
        0.5, 0.5,
        "Awaiting complete round-%d dataset\n%d/%d matched pairs available" % (
            comparison_round, available_pairs, expected_pairs,
        ),
        transform=axis.transAxes, ha="center", va="center",
        fontsize=11, color="0.45",
    )


def make_plot(path, campaign_pairs, available_pairs, snapshot, dpi,
              sif_enabled, combined=False, comparison_round=4):
    limits = {
        key: shared_limit(campaign_pairs, key) for key, _, _ in METRICS
    }
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 10.0))
    for column, (round_label, pairs) in enumerate(campaign_pairs.items()):
        for row, (key, title, unit) in enumerate(METRICS):
            panel_title = "%s — %s" % (round_label, title)
            if pairs:
                draw_panel(
                    axes[row, column], pairs, key, panel_title, unit,
                    limits[key], combined=combined,
                )
            else:
                draw_placeholder(
                    axes[row, column], panel_title, unit,
                    available_pairs[round_label],
                    EXPECTED_PAIRS * (2 if combined else 1),
                    comparison_round=comparison_round,
                )

    # Figure.supxlabel/supylabel are unavailable on older Matplotlib releases
    # still used on the campaign hosts.
    figure.text(
        0.535, 0.090, "Corrected retrieval − truth",
        ha="center", va="center", fontsize=13,
    )
    figure.text(
        0.025, 0.52, "Uncorrected retrieval − truth",
        ha="center", va="center", rotation="vertical", fontsize=13,
    )
    legend = [
        *[
            Line2D(
                [0], [0], marker="o", linestyle="none", markersize=7,
                markerfacecolor=color, markeredgecolor=color, alpha=0.8,
                label=surface.capitalize(),
            )
            for surface, color in SURFACE_COLORS.items()
        ],
        Line2D([0], [0], marker="o", linestyle="none", markersize=7,
               markerfacecolor="0.45", markeredgecolor="0.45",
               alpha=0.65, label="No aerosol (filled)"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=7,
               markerfacecolor="none", markeredgecolor="0.35",
               label="AOD$_{760}=0.28$ (open)"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=7,
               markerfacecolor="0.55", markeredgecolor="0.55",
               label="Perturbations 01–10 (circle)"),
        Line2D([0], [0], marker="s", linestyle="none", markersize=6,
               markerfacecolor="0.55", markeredgecolor="black",
               alpha=0.72,
               label="Perturbation 11 (square)"),
        Line2D([0], [0], marker="x", linestyle="none", markersize=7,
               markeredgewidth=1.3, color="#B2182B",
               label="Uncorr. failed"),
        Line2D([0], [0], color="0.35", linestyle="--", linewidth=1.0,
               label="Equal corrected/uncorrected error"),
    ]
    figure.legend(
        handles=legend, loc="lower center", ncol=5, frameon=False,
        fontsize=9.5, bbox_to_anchor=(0.5, 0.006),
    )
    figure.suptitle(
        "%s retrieval truth errors: round 3 versus round %d" % (
            "Combined SIF-on and no-SIF" if combined else
            "SIF-on" if sif_enabled else "No-SIF", comparison_round,
        ),
        fontsize=15, y=0.975,
    )
    if combined:
        figure.text(
            0.535, 0.937,
            "SIF-on and no-SIF pooled; SIF-on markers contain a small +; "
            "no-SIF markers are plain",
            ha="center", va="center", fontsize=10, color="0.25",
        )
    figure.text(
        0.99, 0.006, "Snapshot %s" % snapshot,
        ha="right", va="bottom", fontsize=7.5, color="0.45",
    )
    figure.subplots_adjust(
        left=0.065, right=0.99, bottom=0.135, top=0.90,
        hspace=0.16, wspace=0.08,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main(comparison_round=4):
    args = parse_args(comparison_round)
    comparison_round = args.compare_round
    snapshot = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    combined = args.sif_enabled == "both"
    sif_enabled = args.sif_enabled == "1"
    sif_modes = (False, True) if combined else (sif_enabled,)
    no_sif_campaigns, sif_on_campaigns = COMPARISON_CAMPAIGNS[comparison_round]
    prefix = "round3_vs_round%d" % comparison_round
    output_root = (
        BOTTOM_ROOT / "comparisons" /
        (prefix + ("_combined_sif" if combined else
                   "_sif" if sif_enabled else "_nosif"))
    )
    stem = (
        prefix + ("_combined_sif" if combined else "_sif" if sif_enabled else "") +
        "_xco2_psurf"
    )
    output = args.output or output_root / (stem + ".png")
    table_output = args.table_output or output_root / (stem + ".dat")

    campaign_pairs = {label: [] for label in no_sif_campaigns}
    available_pairs = {label: 0 for label in no_sif_campaigns}
    for mode in sif_modes:
        campaign_specs = sif_on_campaigns if mode else no_sif_campaigns
        for label, (root, truth_table, accepted_modes) in campaign_specs.items():
            allow_incomplete = mode and label == "Round %d" % comparison_round
            pairs, count = load_campaign(
                label, root, truth_table, accepted_modes,
                allow_incomplete=allow_incomplete,
                expected_state_model={
                    "Round 5": "round5_fixed_sif",
                    "Round 6": "round6_fixed_sif",
                }.get(label),
            )
            available_pairs[label] += count
            if count != EXPECTED_PAIRS:
                if not allow_incomplete:
                    raise RuntimeError(
                        "%s SIF=%d has %d matched pairs; expected %d" % (
                            label, int(mode), count, EXPECTED_PAIRS,
                        )
                    )
                if not args.allow_partial_comparison:
                    pairs = []
            for pair in pairs:
                pair["sif_enabled"] = mode
            campaign_pairs[label].extend(pairs)
    if combined and not args.allow_partial_comparison:
        for label, pairs in campaign_pairs.items():
            if len(pairs) != 2 * EXPECTED_PAIRS:
                campaign_pairs[label] = []

    write_table(table_output, campaign_pairs, snapshot, sif_enabled, combined,
                comparison_round=comparison_round)
    make_plot(
        output, campaign_pairs, available_pairs, snapshot, args.dpi,
        sif_enabled, combined,
        comparison_round=comparison_round,
    )

    for label, pairs in campaign_pairs.items():
        uncorrected_failed = sum(
            not pair["uncorrected"]["converged"] for pair in pairs
        )
        print("%s: %d matched pairs; %d uncorrected failures" % (
            label, len(pairs), uncorrected_failed,
        ))
    print(output)
    print(table_output)


if __name__ == "__main__":
    main()
