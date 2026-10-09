"""Shared CO2 labels for full-column and bottom-layer retrieval figures."""

import math


def co2_case_label(bottom_co2_ppm, xco2_ppm, truth=False, digits=4):
    """Return a visible truth-coordinate label for one retrieval scene."""
    xco2_ppm = float(xco2_ppm)
    prefix = "truth " if truth else ""
    if bottom_co2_ppm is None or not math.isfinite(float(bottom_co2_ppm)):
        return "%sXCO$_2$=%g ppm" % (prefix, xco2_ppm)
    bottom_co2_ppm = float(bottom_co2_ppm)
    return (
        "%sbottom-layer CO$_2$=%g ppm $\\rightarrow$ "
        "column XCO$_2$=%.*f ppm" % (
            prefix, bottom_co2_ppm, int(digits), xco2_ppm
        )
    )


def bottom_layer_mapping_note(pairs):
    """Format the unique bottom-layer-to-column CO2 mapping in ``pairs``.

    ``pairs`` is an iterable of ``(bottom_co2_ppm, xco2_ppm)`` values. Missing
    and non-finite bottom-layer values are ignored so the same plot code can
    continue to serve the historical full-column campaign.
    """
    unique = {}
    for bottom_co2_ppm, xco2_ppm in pairs:
        if bottom_co2_ppm is None:
            continue
        bottom_co2_ppm = float(bottom_co2_ppm)
        xco2_ppm = float(xco2_ppm)
        if not (math.isfinite(bottom_co2_ppm) and math.isfinite(xco2_ppm)):
            continue
        key = round(bottom_co2_ppm, 9)
        if key in unique and not math.isclose(
                unique[key], xco2_ppm, rel_tol=0.0, abs_tol=5.0e-7):
            raise ValueError(
                "one bottom-layer CO2 value maps to inconsistent XCO2 values"
            )
        unique[key] = xco2_ppm
    if not unique:
        return None
    entries = [
        "%g $\\rightarrow$ %.4f" % (bottom, unique[bottom])
        for bottom in sorted(unique)
    ]
    return (
        "Bottom-layer CO$_2$ VMR $\\rightarrow$ column XCO$_2$ (ppm): "
        + ", ".join(entries)
    )


def bottom_layer_mapping_from_rows(rows):
    """Return the mapping note for truth-table row dictionaries."""
    return bottom_layer_mapping_note([
        (row.get("bottom_co2_ppm"), row["xco2_ppm"])
        for row in rows
    ])


def bottom_layer_mapping_from_truth_table(path):
    """Read an ASCII truth table and format its complete bottom-CO2 mapping."""
    columns = None
    rows = []
    with open(str(path), "r", encoding="utf-8") as stream:
        for line in stream:
            text = line.strip()
            if text.startswith("# index "):
                columns = text[2:].split()
            elif text and not text.startswith("#"):
                if columns is None:
                    raise RuntimeError(
                        "truth-table column header was not found in %s" % path
                    )
                values = text.split()
                if len(values) != len(columns):
                    raise RuntimeError("malformed truth-table row: %s" % text)
                rows.append(dict(zip(columns, values)))
    return bottom_layer_mapping_from_rows(rows)
