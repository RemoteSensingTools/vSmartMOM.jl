# HITRAN Spectroscopic Data Management

!!! note "HITRAN data management remains supported"
    The HITRAN data management described here (`artifact`, `fetch_hitran`, edition
    selection and the download cache) remains supported: the RT pipeline reads its
    line lists through `artifact()` and parses them with
    [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl).
    Only the standalone line-by-line cross-section model in `vSmartMOM.Absorption`
    is deprecated in favour of AtmosphericAbsorption.jl.

vSmartMOM.jl requires spectroscopic line parameters from the [HITRAN database](https://hitran.org) to compute gas absorption cross-sections. This page explains the two available data pathways, how to switch between them, and how version provenance is tracked.

## Overview

There are two pathways for obtaining HITRAN data:

| | Legacy Artifacts (default) | Direct Download |
|---|---|---|
| **HITRAN edition** | HITRAN 2016 | Current edition on hitran.org (HITRAN 2024 as of early 2025) |
| **Data source** | Pre-packaged tarballs hosted on a Caltech server | Live download from the [HITRAN API](https://hitran.org/lbl/api) |
| **Isotopologues** | All isotopologues per molecule | All isotopologues per molecule |
| **Wavenumber range** | Full range (0--150,000 cm⁻¹) | Configurable per download |
| **Version tracking** | Implicit (URL path contains `hitran_2016`) | Explicit `.meta.toml` file with SHA-256, download date, source URL |
| **Storage** | Julia Artifacts cache (`~/.julia/artifacts/`) | Scratch space (`~/.julia/scratchspaces/`) |
| **Internet required** | Only on first access (lazy download) | Only on first access (cached thereafter) |

**Backward compatibility:** The default behavior is unchanged. Existing code that calls `artifact("CO2")` or uses YAML configs with molecule names will continue to use the legacy HITRAN 2016 artifacts. You must explicitly opt in to a different edition.

## Pathway 1: Legacy Artifacts (HITRAN 2016) -- Default

This is the original data pathway and remains the default. HITRAN 2016 line-parameter files (`.par` format) are distributed as Julia lazy artifacts. They are downloaded automatically on first use from a Caltech server and cached in `~/.julia/artifacts/`.

No setup is required. This pathway is active by default when the package loads.

### How it works

```julia
using vSmartMOM
using vSmartMOM.Absorption

# artifact() returns the path to the .par file for the requested molecule.
# On first call, the artifact is downloaded and cached locally.
co2_path = artifact("CO2")

# Read the line parameters (all isotopologues, filtered by wavenumber range)
hitran_data = read_hitran(co2_path, mol=2, ν_min=6000, ν_max=6400)

# Check which isotopologues are included
println(sort(unique(hitran_data.iso)))  # e.g., [1, 2, 3, 4, 5, 6, 7, 8]
```

### Limitations

- The data is HITRAN 2016. There is no way to update to a newer edition without switching to the direct download pathway.
- Version provenance is not explicitly tracked -- the only indication of the edition is in the artifact download URL.

## Pathway 2: Direct Download from hitran.org

This pathway downloads data directly from the [HITRAN line-by-line API](https://hitran.org/lbl/api), similar to how the Python [HAPI](https://github.com/hitranonline/hapi) library works. Data is cached locally in a Julia scratch space with full provenance metadata.

### Switching to the direct download pathway

```julia
using vSmartMOM

# Switch the active edition. This affects all subsequent artifact() calls.
set_hitran_edition!("HITRAN2024")

# Now artifact() will use the direct download pathway.
# On first call for each molecule, data is downloaded from hitran.org and cached.
co2_path = artifact("CO2")
```

After this call, the CO2 line-parameter file is stored at:
```
~/.julia/scratchspaces/<package-uuid>/hitran_data/HITRAN2024/CO2.par
```

alongside a metadata file:
```
~/.julia/scratchspaces/<package-uuid>/hitran_data/HITRAN2024/CO2.meta.toml
```

### Downloading with custom wavenumber ranges

For large molecules (e.g., H2O), downloading the full 0--150,000 cm⁻¹ range can be slow. You can restrict the download to the spectral region you need:

```julia
# Download only the O2 A-band region, under its own edition label
fetch_hitran("O2"; numin=12900, numax=13200, edition="HITRAN2024_O2A")

# Download CO2 for the 1.6 um band
fetch_hitran("CO2"; numin=6000, numax=6400, edition="HITRAN2024_1p6um")
```

!!! warning "artifact() requires the full dataset"
    A cached file is reused only when its metadata covers the requested band and
    isotopologues and its SHA-256 matches the file. `artifact()` (and therefore
    `model_from_parameters`) needs the full dataset (0--150,000 cm⁻¹, all
    isotopologues), so a narrower file fetched under the same edition label is
    replaced by a full download on the next `artifact()` call. If that download
    fails, `artifact()` raises an error instead of using the narrower file. Keep
    narrow-band files under a separate edition label if you need them.

### Downloading by explicit isotopologue IDs

If you need only specific isotopologues, use the lower-level `fetch_hitran_by_ids`:

```julia
# Download only the two main CO2 isotopologues (global IDs 7 and 8)
fetch_hitran_by_ids("CO2_main", [7, 8]; numin=6000, numax=6400, edition="HITRAN2024")
```

Global isotopologue IDs can be looked up with:
```julia
using vSmartMOM.Absorption
Absorption.mol_globalID(2, 1)  # CO2, isotopologue 1 -> global ID 7
Absorption.mol_globalID(2, 2)  # CO2, isotopologue 2 -> global ID 8
```

### Re-downloading (force refresh)

To re-download data (e.g., after a HITRAN database update), use `force=true`:

```julia
fetch_hitran("CO2"; edition="HITRAN2024", force=true)
```

## Version Tracking and Provenance

Every molecule downloaded via the direct pathway has a companion `.meta.toml` file that records:

```toml
# HITRAN download metadata
# Auto-generated by vSmartMOM.jl -- do not edit

molecule = "CO2"
edition = "HITRAN2024"
download_date = "2026-04-01T13:16:47.011"
numin = 0
numax = 150000
iso_ids = "7,8,9,10,11,12,13,14,121,15,120,122"
source_url = "https://hitran.org/lbl/api?iso_ids_list=7,8,9,10,11,12,13,14,121,15,120,122&numin=0&numax=150000"
sha256 = "ad576a2ac2f32619910ceb6af4e56cdcf78c3656e6883d3ebe705bf96b837874"
file_size_bytes = 755251
```

This metadata is what makes reuse verifiable. A cached file is used only if
the metadata names the same molecule label and edition, covers the requested
wavenumber range, lists every requested isotopologue, and records a size and
SHA-256 that match the file. Missing, unreadable or inconsistent metadata
(for example, after an interrupted update or a hand-edited file) is treated as
"not verified": `fetch_hitran` downloads again, and `artifact()` downloads
again or fails with an explanation. A cache never silently falls back to an
unverified file.

Downloads are written to a temporary file in the same directory and checked
before they replace the cache: the body must be non-empty, every line must be a
well-formed 160-column HITRAN record, and every line must belong to a requested
isotopologue inside the requested range. A transport error is a failure even
if the server had already sent a successful HTTP status, with one exception:
hitran.org declares a length that counts two bytes per line ending but sends
one, so a complete response ends with "N bytes missing", N being the number of
lines. That response is accepted only when the declared length equals the
received size plus the number of lines, the body ends with a newline, and every
record passes the checks above. If any step fails, the previously cached file
and its metadata are left unchanged.

You can query this programmatically:

```julia
set_hitran_edition!("HITRAN2024")
info = hitran_info("CO2")
println(info["sha256"])         # file integrity hash
println(info["download_date"])  # when it was downloaded
println(info["source_url"])     # exact API query used
```

### Edition labels

The "edition" is a user-assigned cache/provenance label. The HITRAN API always serves the current database edition and has no version selector, so the label is simply a name for the local cache directory: a label such as `"HITRAN2020"` does not prove that the cached lines come from HITRAN 2020. The SHA-256, download date and source URL in the metadata record what was actually downloaded. The default `"artifact"` edition is different: it is the HITRAN 2016 Julia artifact, whose content is fixed by the package's `Artifacts.toml`. You can use any label you like:

```julia
# These are all valid edition labels
fetch_hitran("CO2"; edition="HITRAN2024")
fetch_hitran("CO2"; edition="test_narrowband", numin=6000, numax=6100)
fetch_hitran("CO2"; edition="production_v3")
```

Use `available_hitran_editions()` to list all editions that exist locally:

```julia
available_hitran_editions()
# ["artifact", "HITRAN2024", "test_narrowband"]
```

## Managing Editions

### Switching between editions

```julia
# Use HITRAN 2016 (legacy artifacts)
set_hitran_edition!("artifact")
path_2016 = artifact("CO2")

# Use HITRAN 2024 (downloaded from hitran.org)
set_hitran_edition!("HITRAN2024")
path_2024 = artifact("CO2")

# Check which edition is active
get_hitran_edition()  # "HITRAN2024"
```

### Checking cache status

```julia
hitran_is_cached("CO2")          # check current edition
hitran_is_cached("H2O", "HITRAN2024")  # check specific edition
```

`hitran_is_cached` returns `true` only when `artifact()` would use the cached
file without downloading, that is, when a verified full dataset is present.

### Available molecules

```julia
using vSmartMOM.Absorption
Absorption.show_molecules()       # list all HITRAN molecule names
Absorption.search_molecules("H")  # search by substring
```

## Integration with the RT Pipeline

The full RT pipeline (`model_from_parameters`) uses `artifact()` internally to obtain HITRAN data. Switching the edition before building the model is all you need:

```julia
using vSmartMOM

# Option A: Use HITRAN 2016 (default, no setup needed)
params = read_parameters("config/my_config.yaml")
model = model_from_parameters(params)
R, T = rt_run(model)

# Option B: Use HITRAN 2024
set_hitran_edition!("HITRAN2024")
params = read_parameters("config/my_config.yaml")
model = model_from_parameters(params)  # downloads HITRAN 2024 data as needed
R, T = rt_run(model)
```

No changes to YAML configuration files are required.

### Reusing parsed spectroscopy

Forward and linearized constructors and `BatchContext` share a process-local
cache of parsed HITRAN line databases, keyed by resolved file path and floating
point type. Rebuilding a scene or changing its wavelength grid reuses these
read-only line arrays. Pressure, temperature, broadening settings and backend
remain properties of each independently constructed absorption model.

The cache retains each database until explicitly cleared or the process exits;
it does not automatically reload data. Switching editions selects a different
path. When `fetch_hitran`, `fetch_hitran_by_ids` or `artifact()` replaces a
downloaded file, the parsed cache is cleared automatically. If you replace a
file at the same path by other means, call `clear_spectroscopy_cache!()` before
constructing a model with the replacement data. Existing models keep their
original data. LUTs loaded into parameters are caller-owned and are unaffected
by this cache.

```@docs
clear_spectroscopy_cache!
```

## HITRAN API Rate Limits

The HITRAN API at hitran.org imposes a daily query limit. If you exceed it, you will see:

```
ERROR: HITRAN API rate limit exceeded. You have exceeded the daily limit of API queries. Try again tomorrow.
```

To minimize API calls:
- Data is cached locally after the first download. Subsequent calls use the cache.
- Use `hitran_is_cached("CO2")` to check whether `artifact()` (and so
  `model_from_parameters`) will need to download the full dataset.
- For standalone work, download only the wavenumber range you need (use
  `numin`/`numax` in `fetch_hitran`) under a separate edition label;
  `artifact()` always requires the full dataset.

## Cleaning Up Cached Data

Downloaded HITRAN data lives in Julia's scratch space system. To remove all cached data:

```julia
using Scratch
Scratch.clear_scratchspaces!()
```

Legacy artifacts are managed by Julia's package manager:
```julia
using Pkg
Pkg.gc()  # removes unused artifacts
```

## Citing HITRAN

The line parameters behind every absorption calculation come from HITRAN. If you
publish results computed with them, please cite the edition you used:

- the default `"artifact"` edition is **HITRAN2016** (Gordon et al. 2017,
  [doi:10.1016/j.jqsrt.2017.06.038](https://doi.org/10.1016/j.jqsrt.2017.06.038));
- data downloaded from hitran.org is the **current edition** at the time of
  download (HITRAN2024 at the time of writing: Gordon et al. 2026,
  [doi:10.1016/j.jqsrt.2026.109807](https://doi.org/10.1016/j.jqsrt.2026.109807)).
  The edition label you choose does not change what the server delivers; the
  `.meta.toml` sidecar records the download date and checksum;
- for your own `.par` files, cite the edition they came from.

Partition-sum, HAPI, collision-induced-absorption and MT_CKD references are on the
[References](../vSmartMOM/References.md) page; BibTeX entries are in `CITATION.bib`.

## Library Reference

The canonical docstrings for `artifact`, `fetch_hitran`,
`fetch_hitran_by_ids`, HITRAN edition preferences, `hitran_info`, and
`hitran_is_cached` are grouped in the [Absorption API](../api/absorption.md).
