#!/usr/bin/env python3
"""Write the physical-ensemble retrieval products as readable Markdown tables.

This is the tabular companion to ``plot_physical_retrieval_ensembles.py``.  It
deliberately reuses that module's loading, validation, pairing, and physical
reconstruction functions rather than reimplementing them, so a table can never
disagree with the figure drawn from the same retrieval files.  In particular:

* the same SIF-provenance and campaign guards reject mismatched inputs;
* the same paired ensemble gate (complete, converged, fit-quality accepted for
  both members of a perturbation index) selects the members;
* perturbation 11 is reported separately and never enters the statistics; and
* ``XCO2``, the dry-air column, and the surface-pressure-dependent geometry are
  recomputed per member exactly as the CO2 figure computes them.

One Markdown file is written per truth CO2 case and aerosol category, named
after the corresponding figure base name.  Every table names its scene by
surface and truth state index, so a partially available campaign produces the
same rows with explicit missing-data markers instead of silently shrinking.

The default ``--all-co2-cases`` sweep is the intended way to run this as a
campaign fills in: cases with no paired retrievals yet are skipped and listed,
and re-running simply refreshes the files that have gained data.
"""

import argparse
import importlib.util
import math
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BACKEND_PATH = HERE / "plot_physical_retrieval_ensembles.py"
TABLE_MARKER_VERSION = 1
MISSING = "n/a"
CLASS_LABELS = {"uncorrected": "Uncorrected", "corrected": "Corrected"}
CLASS_ORDER = ("uncorrected", "corrected")


def load_backend():
    """Import the plotting backend as a module without running its ``main``."""
    if not BACKEND_PATH.is_file():
        raise RuntimeError(
            "the physical-ensemble plotting backend is missing: %s" %
            BACKEND_PATH
        )
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    spec = importlib.util.spec_from_file_location(
        "_physical_ensemble_plot_backend", BACKEND_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


backend = load_backend()

DEFAULT_OUTPUT_DIR = HERE / "physical_ensemble_tables"
CATEGORY_TEXT = {
    "no_aerosol": "no aerosol",
    "with_aerosol": "aerosol AOD760 = 0.28",
}


# ----------------------------------------------------------------------------
# formatting helpers
# ----------------------------------------------------------------------------

def fmt(value, digits=4, signed=False):
    """Format one finite float, or the explicit missing marker."""
    if value is None:
        return MISSING
    value = float(value)
    if not np.isfinite(value):
        return MISSING
    return "%+.*f" % (digits, value) if signed else "%.*f" % (digits, value)


def fmt_significant(value, digits=6):
    if value is None or not np.isfinite(float(value)):
        return MISSING
    return "%.*g" % (digits, float(value))


def fmt_stat(values, digits=4, signed=False):
    """Format ``mean ± sample SD`` for an ensemble, degrading gracefully."""
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return MISSING
    mean = fmt(float(np.mean(array)), digits, signed)
    if array.size < 2:
        return mean
    return "%s ± %s" % (mean, fmt(float(np.std(array, ddof=1)), digits))


def render_table(headers, rows, aligns=None):
    """Return a column-aligned GitHub-flavoured Markdown table."""
    if aligns is None:
        aligns = "l" * len(headers)
    if len(aligns) != len(headers):
        raise ValueError("alignment string does not match the header count")
    cells = [[str(value) for value in headers]]
    cells.extend([str(value) for value in row] for row in rows)
    for row in cells:
        if len(row) != len(headers):
            raise ValueError("a table row does not match the header count")
    widths = [
        max(len(row[column]) for row in cells)
        for column in range(len(headers))
    ]
    widths = [max(width, 3) for width in widths]

    def line(row):
        padded = [
            row[column].rjust(widths[column]) if aligns[column] == "r"
            else row[column].ljust(widths[column])
            for column in range(len(headers))
        ]
        return "| " + " | ".join(padded) + " |"

    separator = "| " + " | ".join(
        ("-" * (widths[column] - 1) + ":") if aligns[column] == "r"
        else ("-" * widths[column])
        for column in range(len(headers))
    ) + " |"
    return "\n".join([line(cells[0]), separator] + [line(r) for r in cells[1:]])


# ----------------------------------------------------------------------------
# physical statistics
# ----------------------------------------------------------------------------

def partial_columns(record, profile):
    """Return the exact partial-column contributions to XCO2 in ppm.

    The three groups follow the retrieved state schema: layers 1--4 are held at
    the campaign's fixed upper-layer VMR, layers 5--15 are retrieved and
    vertically correlated, and layer 16 is the bottom layer that carries the
    truth perturbation in the bottom-layer campaigns.  The three contributions
    sum to XCO2 by construction.
    """
    co2_ppm = np.asarray(record["co2_ppm"], dtype=float)
    dry = np.asarray(profile["dry_air_column"], dtype=float)
    total = float(np.sum(dry))
    groups = {
        "fixed": slice(0, 4),
        "middle": slice(4, 15),
        "bottom": slice(15, 16),
    }
    return {
        name: float(np.sum(co2_ppm[index] * dry[index]) / total)
        for name, index in groups.items()
    }


def record_statistics(record, atmosphere):
    """Reduce one physical record to the scalar quantities the tables report."""
    profile = backend.co2_molecular_profile(record, atmosphere)
    co2_ppm = np.asarray(record["co2_ppm"], dtype=float)
    dry = np.asarray(profile["dry_air_column"], dtype=float)
    stats = {
        "xco2_ppm": float(
            record.get("saved_xco2_ppm", profile["xco2_ppm"])
        ),
        "psurf": float(record["psurf"]),
        "co2_column_1e21": profile["total_co2_column"] / 1.0e21,
        "bottom_ppm": float(co2_ppm[15]),
        "middle_ppm": float(
            np.sum(co2_ppm[4:15] * dry[4:15]) / np.sum(dry[4:15])
        ),
        "fixed_ppm": float(co2_ppm[0]),
        "aod_total": float(sum(record["aod"].values())),
        "sif_radiance760": backend.sif_reference_radiance(record),
        "sif_angular760": backend.sif_reference_angular_integral(record),
    }
    stats.update(
        {"partial_%s" % name: value
         for name, value in partial_columns(record, profile).items()}
    )
    for species, _, _, _, _ in backend.SPECIES:
        stats["aod_%s" % species] = float(record["aod"][species])
        stats["z0_%s" % species] = float(record["z0"][species])
    for state_band, _, _, _, _ in backend.BANDS:
        for order in range(3):
            stats["surface_%s_P%d" % (state_band, order)] = float(
                record["surface"][state_band][order]
            )
    return stats


def scene_statistics(card, atmosphere):
    """Return truth, per-class ensembles, and availability for one scene."""
    ensemble = card["ensemble"]
    scene = {
        "state_index": card["state_index"],
        "count": len(ensemble["indices"]),
        "indices": list(ensemble["indices"]),
        "missing_indices": list(ensemble["missing_indices"]),
        "valid_counts": dict(ensemble["valid_counts"]),
        "noiseless_paired": bool(ensemble["noiseless_paired"]),
        "fit_failures": {
            retrieval_class: list(
                ensemble.get("fit_failure_indices", {}).get(retrieval_class, [])
            )
            for retrieval_class in CLASS_ORDER
        },
        "truth": record_statistics(card["truth"], atmosphere),
        "truth_sif_case": card["truth"]["sif"],
        "classes": {},
        "noiseless": {},
    }
    for retrieval_class in CLASS_ORDER:
        scene["classes"][retrieval_class] = [
            record_statistics(record, atmosphere)
            for record in ensemble[retrieval_class]
        ]
        member = ensemble["noiseless"][retrieval_class]
        scene["noiseless"][retrieval_class] = (
            record_statistics(member, atmosphere) if member is not None
            else None
        )
    return scene


def series(scene, retrieval_class, key):
    return [member[key] for member in scene["classes"][retrieval_class]]


def mean_of(scene, retrieval_class, key):
    values = series(scene, retrieval_class, key)
    return float(np.mean(values)) if values else None


# ----------------------------------------------------------------------------
# table sections
# ----------------------------------------------------------------------------

def xco2_table(scenes):
    rows = []
    for surface, scene in scenes:
        truth_xco2 = scene["truth"]["xco2_ppm"]
        unc = series(scene, "uncorrected", "xco2_ppm")
        corr = series(scene, "corrected", "xco2_ppm")
        unc_mean = float(np.mean(unc)) if unc else None
        corr_mean = float(np.mean(corr)) if corr else None
        rows.append([
            surface.title(),
            "%03d" % scene["state_index"],
            "%d/10" % scene["count"],
            fmt(truth_xco2),
            fmt_stat(unc),
            fmt_stat(corr),
            fmt(None if unc_mean is None else unc_mean - truth_xco2,
                signed=True),
            fmt(None if corr_mean is None else corr_mean - truth_xco2,
                signed=True),
            fmt(None if (unc_mean is None or corr_mean is None)
                else corr_mean - unc_mean, signed=True),
        ])
    return render_table(
        ["Scene", "State", "n", "Truth (ppm)", "Uncorrected (ppm)",
         "Corrected (ppm)", "Unc bias", "Corr bias", "Corr - unc"],
        rows, "llrrrrrrr",
    )


DECOMPOSITION_KEYS = (
    "partial_fixed", "partial_middle", "partial_bottom", "xco2_ppm",
    "fixed_ppm", "middle_ppm", "bottom_ppm",
)


def decomposition_table(scenes):
    rows = []
    for surface, scene in scenes:
        entries = [("Truth", scene["truth"])]
        entries.extend(
            (CLASS_LABELS[retrieval_class],
             {key: mean_of(scene, retrieval_class, key)
              for key in DECOMPOSITION_KEYS})
            for retrieval_class in CLASS_ORDER
        )
        for label, values in entries:
            rows.append([
                surface.title(),
                label,
                fmt(values["partial_fixed"]),
                fmt(values["partial_middle"]),
                fmt(values["partial_bottom"]),
                fmt(values["xco2_ppm"]),
                fmt(values["fixed_ppm"], 3),
                fmt(values["middle_ppm"], 3),
                fmt(values["bottom_ppm"], 3),
            ])
    return render_table(
        ["Scene", "Class", "L1-4 (ppm)", "L5-15 (ppm)", "L16 (ppm)",
         "sum = XCO2", "L1-4 VMR", "L5-15 VMR", "L16 VMR"],
        rows, "llrrrrrrr",
    )


def light_path_table(scenes):
    rows = []
    for surface, scene in scenes:
        truth = scene["truth"]
        rows.append([
            surface.title(), "Truth",
            fmt(truth["aod_total"]),
            fmt(truth["psurf"], 3),
            fmt(truth["co2_column_1e21"], 4),
            fmt(0.0, 4, signed=True),
        ])
        for retrieval_class in CLASS_ORDER:
            xco2 = mean_of(scene, retrieval_class, "xco2_ppm")
            rows.append([
                surface.title(), CLASS_LABELS[retrieval_class],
                fmt_stat(series(scene, retrieval_class, "aod_total")),
                fmt_stat(series(scene, retrieval_class, "psurf"), 3),
                fmt_stat(series(scene, retrieval_class, "co2_column_1e21"), 4),
                fmt(None if xco2 is None else xco2 - truth["xco2_ppm"],
                    signed=True),
            ])
    return render_table(
        ["Scene", "Class", "Total AOD760", "p_surf (hPa)",
         "N_CO2 (1e21 cm-2)", "XCO2 bias (ppm)"],
        rows, "llrrrr",
    )


def aerosol_table(scenes):
    rows = []
    for surface, scene in scenes:
        for species, species_label, _, _, _ in backend.SPECIES:
            aod_key = "aod_%s" % species
            z0_key = "z0_%s" % species
            rows.append([
                surface.title(), species_label,
                fmt(scene["truth"][aod_key]),
                fmt_stat(series(scene, "uncorrected", aod_key)),
                fmt_stat(series(scene, "corrected", aod_key)),
                fmt(scene["truth"][z0_key], 3),
                fmt_stat(series(scene, "uncorrected", z0_key), 3),
                fmt_stat(series(scene, "corrected", z0_key), 3),
            ])
    return render_table(
        ["Scene", "Species", "Truth AOD760", "Unc AOD760", "Corr AOD760",
         "Truth z0 (km)", "Unc z0 (km)", "Corr z0 (km)"],
        rows, "llrrrrrr",
    )


def surface_table(scenes):
    rows = []
    for surface, scene in scenes:
        for state_band, _, band_label, _, _ in backend.BANDS:
            plain_label = (
                band_label.replace("$_2$", "2").replace("$_{2}$", "2")
            )
            for order in range(3):
                key = "surface_%s_P%d" % (state_band, order)
                unc = series(scene, "uncorrected", key)
                corr = series(scene, "corrected", key)
                rows.append([
                    surface.title(), plain_label, "P%d" % order,
                    fmt_significant(scene["truth"][key]),
                    fmt_significant(float(np.mean(unc))) if unc else MISSING,
                    fmt_significant(float(np.mean(corr))) if corr else MISSING,
                ])
    return render_table(
        ["Scene", "Band", "Order", "Truth", "Uncorrected mean",
         "Corrected mean"],
        rows, "lllrrr",
    )


def sif_table(scenes):
    rows = []
    for surface, scene in scenes:
        truth_sif = scene["truth_sif_case"]
        rows.append([
            surface.title(),
            fmt(scene["truth"]["sif_radiance760"], 5, signed=True),
            fmt_stat(series(scene, "uncorrected", "sif_radiance760"), 5,
                     signed=True),
            fmt_stat(series(scene, "corrected", "sif_radiance760"), 5,
                     signed=True),
            fmt(truth_sif.get("angular_integral760"), 3, signed=True),
            "%s / %s" % (
                truth_sif.get("sif760_status", "?"),
                truth_sif.get("msif_status", "?"),
            ),
        ])
    return render_table(
        ["Scene", "Truth L760", "Unc L760", "Corr L760", "Truth 2pi*L760",
         "SIF760 / mSIF state"],
        rows, "lrrrrl",
    )


def availability_table(scenes):
    rows = []
    for surface, scene in scenes:
        missing = (
            ", ".join("%02d" % index for index in scene["missing_indices"])
            or MISSING
        )
        failures = []
        for retrieval_class in CLASS_ORDER:
            indices = scene["fit_failures"][retrieval_class]
            if indices:
                failures.append(
                    "%s %s" % (
                        "unc" if retrieval_class == "uncorrected" else "corr",
                        "/".join("%02d" % index for index in indices),
                    )
                )
        rows.append([
            surface.title(),
            "%03d" % scene["state_index"],
            "%d/10" % scene["count"],
            "%d/10" % scene["valid_counts"]["uncorrected"],
            "%d/10" % scene["valid_counts"]["corrected"],
            "paired" if scene["noiseless_paired"] else MISSING,
            missing,
            "; ".join(failures) or MISSING,
        ])
    return render_table(
        ["Scene", "State", "Paired", "Unc valid", "Corr valid", "p11",
         "Missing pairs", "Fit-quality failures"],
        rows, "llrrrlll",
    )


def noiseless_table(scenes):
    rows = []
    for surface, scene in scenes:
        truth_xco2 = scene["truth"]["xco2_ppm"]
        cells = []
        for retrieval_class in CLASS_ORDER:
            member = scene["noiseless"][retrieval_class]
            if member is None:
                cells.extend([MISSING, MISSING])
            else:
                cells.extend([
                    fmt(member["xco2_ppm"]),
                    fmt(member["xco2_ppm"] - truth_xco2, signed=True),
                ])
        rows.append([surface.title(), fmt(truth_xco2)] + cells)
    return render_table(
        ["Scene", "Truth (ppm)", "Unc p11 (ppm)", "Unc p11 bias",
         "Corr p11 (ppm)", "Corr p11 bias"],
        rows, "lrrrrr",
    )


# ----------------------------------------------------------------------------
# document assembly
# ----------------------------------------------------------------------------

NOTES = """\
## Notes

- Statistics use matched perturbations 01--10 for which both the corrected and
  uncorrected retrieval is complete, state-step converged, and accepted by the
  saved spectral fit-quality gate (unless `--include-fit-failures` was passed,
  which is recorded in the metadata above).  `n` is the number of such pairs.
- Values are the ensemble mean followed by one sample standard deviation
  (`ddof=1`); a single member is reported without a spread, and `""" + MISSING + """`
  marks a quantity with no available member.  Table cells are plain ASCII so
  the raw file stays column-aligned in any editor.
- Perturbation 11 is the noiseless retrieval.  It never enters the statistics
  and appears only in its own table.
- `XCO2` is the dry-air-weighted column average.  The dry-air column and the
  geometric layer grid are recomputed at each member's retrieved surface
  pressure before the average is formed, exactly as in the CO2 figure.
- The partial-column contributions `L1-4`, `L5-15`, and `L16` are dry-air
  weighted and sum to `XCO2` by construction, so they show where a column bias
  is actually produced.  Layers 1--4 are held at the campaign's fixed upper CO2
  VMR; layers 5--16 are retrieved.  `VMR` columns are dry-air-weighted mean
  layer volume mixing ratios in ppm over the same groups.
- Aerosol `z0` is the retrieved lognormal median altitude, not the extinction
  mode; the fixed geometric widths come from the campaign's aerosol table.
- Surface coefficients are the Legendre coefficients on the canonical
  coefficient-definition grid, before the strong-band convolution shoulder.
- `L760` is the BOA SIF spectral radiance at 760 nm in
  mW m^-2 sr^-1 nm^-1, and `2*pi*L760` is the campaign's angular-integral
  normalization coordinate.
- "Corrected" fits the Rayleigh truth simulation; "uncorrected" fits the
  Cabannes+RRS truth simulation with a forward model that does not represent
  RRS.  Both use the linearized `RS_type::noRS` forward model.
"""


def document(category, cards, selection, atmosphere, args, source_root):
    scenes = [
        (surface, scene_statistics(cards[surface], atmosphere))
        for surface in backend.SURFACE_ORDER
        if surface in cards
    ]
    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    co2_text = (
        "bottom-layer CO₂ = %g ppm → column XCO₂ = %.4f ppm" % (
            selection["value"], selection["xco2_ppm"]
        )
        if selection["slug"] == "bottom_co2" else
        "column XCO₂ = %g ppm" % selection["value"]
    )
    ascii_co2_text = co2_text.replace("→", "->").replace("₂", "2")
    campaign = selection["campaign_label"] or "unlabelled campaign"
    title = "%s — %s, %s" % (
        campaign, co2_text, CATEGORY_TEXT.get(category, category)
    )
    marker = (
        "<!-- physical-ensemble-table v%d; selection=%s; co2_ppm=%.6f; "
        "xco2_ppm=%.6f; sif=%s; aerosol=%s; campaign=%s; generated=%s -->" % (
            TABLE_MARKER_VERSION, selection["slug"], selection["value"],
            selection["xco2_ppm"], selection["sif_case"], category,
            campaign.replace("-->", "->").replace(";", ","), generated,
        )
    )
    metadata = render_table(
        ["Field", "Value"],
        [
            ["Campaign", campaign],
            ["Truth CO2 case", ascii_co2_text],
            ["Aerosol category", CATEGORY_TEXT.get(category, category)],
            ["Truth SIF case", selection["sif_case"]],
            ["Retrieval root", "`%s`" % source_root],
            ["Truth table", "`%s`" % args.truth_table],
            ["Atmosphere", "`%s`" % args.atmospheric_profile],
            ["Fit failures included",
             "yes" if args.include_fit_failures else "no"],
            ["Availability snapshot", generated],
        ],
        "ll",
    )
    if selection.get("state_description"):
        metadata += "\n\nState model: %s\n" % selection["state_description"]
    sections = [
        marker,
        "",
        "# %s" % title,
        "",
        metadata,
        "",
        "## 1. XCO₂ ensemble summary",
        "",
        xco2_table(scenes),
        "",
        "## 2. Where the column bias is produced",
        "",
        decomposition_table(scenes),
        "",
        "## 3. Light-path diagnostics",
        "",
        light_path_table(scenes),
        "",
        "## 4. Aerosol state by species",
        "",
        aerosol_table(scenes),
        "",
        "## 5. Surface Legendre coefficients",
        "",
        surface_table(scenes),
        "",
        "## 6. O₂ A-band SIF state",
        "",
        sif_table(scenes),
        "",
        "## 7. Noiseless retrieval (perturbation 11)",
        "",
        noiseless_table(scenes),
        "",
        "## 8. Ensemble availability",
        "",
        availability_table(scenes),
        "",
        NOTES,
    ]
    populated = any(scene["count"] > 0 for _, scene in scenes)
    return "\n".join(sections).rstrip() + "\n", populated


MARKER_PATTERN = re.compile(
    r"<!--\s*physical-ensemble-table v(?P<version>\d+);\s*(?P<fields>.*?)-->"
)


def read_marker(path):
    try:
        with path.open("r", encoding="utf-8") as stream:
            head = stream.read(4096)
    except OSError:
        return None
    match = MARKER_PATTERN.search(head)
    if match is None:
        return None
    fields = {}
    for item in match.group("fields").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            fields[key.strip()] = value.strip()
    fields["version"] = match.group("version")
    return fields


def marker_float(text, digits):
    """Format a numeric field recovered from a file marker."""
    try:
        formatted = "%.*f" % (digits, float(text))
    except (TypeError, ValueError):
        return MISSING
    return formatted.rstrip("0").rstrip(".") if "." in formatted else formatted


def write_index(output_dir):
    """Rebuild ``index.md`` from every table file currently in the directory."""
    entries = []
    for path in sorted(output_dir.glob("*.md")):
        if path.name == "index.md":
            continue
        marker = read_marker(path)
        if marker is None:
            continue
        entries.append((path, marker))
    if not entries:
        return None

    def sort_key(entry):
        marker = entry[1]
        try:
            value = float(marker.get("co2_ppm", "nan"))
        except ValueError:
            value = float("nan")
        return (
            marker.get("campaign", ""),
            marker.get("sif", ""),
            value if np.isfinite(value) else 0.0,
            marker.get("aerosol", ""),
        )

    entries.sort(key=sort_key)
    rows = []
    for path, marker in entries:
        rows.append([
            "[`%s`](%s)" % (path.name, path.name),
            marker.get("campaign", MISSING),
            marker_float(marker.get("co2_ppm"), 6),
            marker_float(marker.get("xco2_ppm"), 4),
            CATEGORY_TEXT.get(marker.get("aerosol", ""),
                              marker.get("aerosol", MISSING)),
            marker.get("sif", MISSING),
            marker.get("generated", MISSING),
        ])
    text = "\n".join([
        "# Physical retrieval ensemble tables",
        "",
        "Markdown companions to the figures in the sibling "
        "`plots*/physical_ensembles` directory.  Regenerate with "
        "`RRS_XCO2/visualization/tabulate_physical_retrieval_ensembles.py`.",
        "",
        render_table(
            ["File", "Campaign", "Truth CO2 (ppm)", "Column XCO2 (ppm)",
             "Aerosol", "SIF", "Generated"],
            rows, "llrrlll",
        ),
        "",
    ])
    index_path = output_dir / "index.md"
    index_path.write_text(text, encoding="utf-8")
    return index_path


# ----------------------------------------------------------------------------
# command line
# ----------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--inversion-root", type=Path, default=HERE,
        help="Directory containing corrected/ and uncorrected/ retrieval files",
    )
    parser.add_argument("--truth-table", type=Path, default=backend.TRUTH_TABLE)
    parser.add_argument(
        "--scene-components", type=Path, default=backend.SCENE_COMPONENTS,
    )
    parser.add_argument(
        "--vertical-profile-table", type=Path,
        default=backend.VERTICAL_PROFILE_TABLE,
    )
    parser.add_argument(
        "--atmospheric-profile", type=Path,
        default=backend.DEFAULT_ATMOSPHERIC_PROFILE,
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--campaign-label", default="")
    parser.add_argument("--xco2-ppm", type=float, default=None)
    parser.add_argument("--bottom-co2-ppm", type=float, default=None)
    parser.add_argument(
        "--all-co2-cases", action="store_true",
        help=(
            "Tabulate every truth CO2 case present in the truth table for the "
            "selected SIF case instead of one explicit value"
        ),
    )
    parser.add_argument(
        "--sif-case", choices=tuple(backend.SIF_CASE_METADATA), default="off",
    )
    parser.add_argument(
        "--aerosol-category", choices=tuple(backend.CATEGORY_SELECTIONS),
        default="both",
    )
    parser.add_argument(
        "--include-fit-failures", action="store_true",
        help=(
            "Include complete, state-step-converged retrievals whose saved "
            "spectral fit-quality flag failed; recorded in each file"
        ),
    )
    parser.add_argument(
        "--include-empty", action="store_true",
        help=(
            "Write a table even when no scene has a paired retrieval yet "
            "(default: skip and report the case as pending)"
        ),
    )
    parser.add_argument(
        "--no-index", action="store_true",
        help="Do not rebuild index.md in the output directory",
    )
    return parser.parse_args(argv)


def _campaign_is_bottom(args):
    rows = backend.read_truth_table(args.truth_table)
    return any("bottom_co2_ppm" in row for row in rows.values())


def selection_values(args):
    """Return ``(uses_bottom, [values])`` for the requested run."""
    uses_bottom = _campaign_is_bottom(args)
    if args.all_co2_cases:
        field = "bottom_co2_ppm" if uses_bottom else "xco2_ppm"
        rows = backend.read_truth_table(args.truth_table)
        values = sorted({
            round(float(row[field]), 9)
            for row in rows.values()
            if row.get("sif_case") == args.sif_case and field in row
        })
        if not values:
            raise RuntimeError(
                "no truth CO2 case in %s matches sif_case=%s" %
                (args.truth_table, args.sif_case)
            )
        return uses_bottom, values
    if uses_bottom:
        if args.xco2_ppm is not None:
            raise SystemExit(
                "this is a bottom-layer campaign; use --bottom-co2-ppm"
            )
        return True, [400.0 if args.bottom_co2_ppm is None
                      else args.bottom_co2_ppm]
    if args.bottom_co2_ppm is not None:
        raise SystemExit("this is a full-column campaign; use --xco2-ppm")
    return False, [400.0 if args.xco2_ppm is None else args.xco2_ppm]


def case_namespace(args, uses_bottom, value):
    """Build the backend selection namespace for one CO2 case."""
    return argparse.Namespace(
        inversion_root=args.inversion_root,
        truth_table=args.truth_table,
        scene_components=args.scene_components,
        vertical_profile_table=args.vertical_profile_table,
        atmospheric_profile=args.atmospheric_profile,
        campaign_label=args.campaign_label,
        sif_case=args.sif_case,
        aerosol_category=args.aerosol_category,
        include_fit_failures=args.include_fit_failures,
        xco2_ppm=None if uses_bottom else value,
        bottom_co2_ppm=value if uses_bottom else None,
    )


def co2_tag(value):
    return (
        str(int(round(value)))
        if math.isclose(value, round(value), abs_tol=1.0e-9)
        else ("%g" % value).replace(".", "p")
    )


def main(argv=None):
    args = parse_args(argv)
    uses_bottom, values = selection_values(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    written = []
    skipped = []
    failed = []
    for value in values:
        case_args = case_namespace(args, uses_bottom, value)
        try:
            categories, _boundaries, atmosphere, selection = (
                backend.prepare_categories(case_args)
            )
        except (RuntimeError, ValueError) as error:
            failed.append((value, str(error)))
            continue
        for category in selection["category_names"]:
            text, populated = document(
                category, categories[category], selection, atmosphere,
                args, args.inversion_root,
            )
            base = "%s_%s_%s_%s" % (
                selection["slug"], co2_tag(selection["value"]),
                selection["sif_slug"], category,
            )
            path = args.output_dir / ("%s.md" % base)
            if not populated and not args.include_empty:
                skipped.append(path.name)
                continue
            path.write_text(text, encoding="utf-8")
            written.append(path)
    for value, reason in failed:
        print("skipped CO2 case %g: %s" % (value, reason), file=sys.stderr)
    if skipped:
        print("pending (no paired retrievals yet): %s" % ", ".join(skipped))
    for path in written:
        print(path)
    if written and not args.no_index:
        index_path = write_index(args.output_dir)
        if index_path is not None:
            print(index_path)
    if not written and not skipped and not failed:
        raise SystemExit("no table was produced")
    return 0


if __name__ == "__main__":
    sys.exit(main())
