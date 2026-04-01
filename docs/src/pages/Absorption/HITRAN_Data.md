# HITRAN Spectroscopic Data Management

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
# Download only the O2 A-band region
fetch_hitran("O2"; numin=12900, numax=13200, edition="HITRAN2024")

# Download CO2 for the 1.6 um band
fetch_hitran("CO2"; numin=6000, numax=6400, edition="HITRAN2024")
```

!!! warning "Wavenumber range is fixed at download time"
    Once a molecule is cached, subsequent `artifact()` calls return the cached file
    regardless of the wavenumber range used during download. If you need a different
    range, either use `force=true` to re-download or use a different edition label.

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

You can query this programmatically:

```julia
set_hitran_edition!("HITRAN2024")
info = hitran_info("CO2")
println(info["sha256"])         # file integrity hash
println(info["download_date"])  # when it was downloaded
println(info["source_url"])     # exact API query used
```

### Edition labels

The "edition" is a user-assigned label. The HITRAN API always serves the current database edition and has no version selector. The edition label is simply a name for the local cache directory. You can use any label you like:

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
params = parameters_from_yaml("config/my_config.yaml")
model = model_from_parameters(params)
R, T = rt_run(model)

# Option B: Use HITRAN 2024
set_hitran_edition!("HITRAN2024")
params = parameters_from_yaml("config/my_config.yaml")
model = model_from_parameters(params)  # downloads HITRAN 2024 data as needed
R, T = rt_run(model)
```

No changes to YAML configuration files are required.

## HITRAN API Rate Limits

The HITRAN API at hitran.org imposes a daily query limit. If you exceed it, you will see:

```
ERROR: HITRAN API rate limit exceeded. You have exceeded the daily limit of API queries. Try again tomorrow.
```

To minimize API calls:
- Data is cached locally after the first download. Subsequent calls use the cache.
- Use `hitran_is_cached("CO2")` to check before triggering a download.
- Download only the wavenumber range you need (use `numin`/`numax` in `fetch_hitran`).

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

## API Reference

```@docs
vSmartMOM.artifact
vSmartMOM.fetch_hitran
vSmartMOM.fetch_hitran_by_ids
vSmartMOM.set_hitran_edition!
vSmartMOM.get_hitran_edition
vSmartMOM.available_hitran_editions
vSmartMOM.hitran_info
vSmartMOM.hitran_is_cached
```
