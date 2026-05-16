# TOMAS sectional aerosol scheme

> Reference page for the TOMAS (TwO-Moment Aerosol Sectional) microphysics
> scheme as consumed by vSmartMOM. The single source of truth for all
> numerical values on this page is
> [`data/aerosols/tomas_species.yaml`](../../../../data/aerosols/tomas_species.yaml) —
> verified against the GCHP / GEOS-Chem upstream code at
> `src/GCHP_GridComp/GEOSChem_GridComp/geos-chem/GeosCore/tomas_mod.F90`.

## What TOMAS is

TOMAS represents the aerosol size distribution as a set of **mass-space
bins**, with two prognostic moments per bin per species:

1. **Particle number** `Nk` (variable `SpeciesConcVV_NK<bin>` in GCHP
   diagnostic output, units `1e3 · num · mol_dry⁻¹`).
2. **Per-species mass mixing ratio** `Mk[bin, species]` (variable
   `SpeciesConcVV_<SPECIES><bin>`, units `mol/mol dry air`).

Both moments evolve under microphysical processes (nucleation, condensation,
coagulation, activation/deposition) plus emissions and chemistry. This page
documents only the **static** scheme metadata that vSmartMOM needs to
ingest a GCHP scene — bin geometry, species set, molar masses, densities,
naming conventions.

## Supported bin counts

GCHP compile-time `#define`s pick one of four variants:

| Variant    | `IBINS` | Bin rule                              | `M₀` (kg)      | Approx Dp range (assuming ρ=1500 kg/m³) |
|------------|--------:|---------------------------------------|---------------:|-----------------------------------------|
| TOMAS12    |      12 | quadrupling, ×32 jump to upper edge   | 1.0 × 10⁻²¹     | 6 nm – 1.6 μm                            |
| **TOMAS15**|      15 | quadrupling, ×32 jump to upper edge   | 1.5625 × 10⁻²³ | 3.5 nm – 6.4 μm                          |
| TOMAS30    |      30 | doubling                              | 1.0 × 10⁻²¹     | 6 nm – 6 μm                              |
| TOMAS40    |      40 | doubling                              | 9.7656 × 10⁻²⁵ | 0.8 nm – 6 μm                            |

**Quadrupling rule** (TOMAS12 / TOMAS15):
```
Xk[k]       = M₀ · 4^(k-1)        for k < IBINS
Xk[IBINS+1] = Xk[IBINS]   · 32    upper-edge jump
```

**Doubling rule** (TOMAS30 / TOMAS40):
```
Xk[k] = M₀ · 2^(k-1)
```

In both cases `Xk` is the **per-particle dry mass** at bin edge `k`,
in kilograms. The bin **center diameter** in metres is

```
Dp = ( mp/ρ · 6/π )^(1/3)
```

where `mp = sqrt(Xk[k] · Xk[k+1])` is the geometric-mean mass and `ρ` is
the bin-averaged density (composition-weighted).

> ⚠ **Bin geometry depends on IBINS.** Do not assume 15. The vSmartMOM
> reader at [`src/Aerosols/schemes/tomas15.jl`](../../../../src/Aerosols/schemes/tomas15.jl)
> dispatches on the variant configured in `TOMASScheme(:tomas15)`,
> `TOMASScheme(:tomas40)`, etc. `TOMAS15Scheme` remains only as a
> compatibility alias for older configs.

## Species set

Mass tracers — these appear as `SpeciesConcVV_<CODE><bin>` per-bin
mixing ratios in GCHP output:

| SRT index | Code  | Long name                            | MOLWT (g/mol) | ρ (kg/m³)              | Hydrophilic |
|----------:|-------|--------------------------------------|--------------:|------------------------|-------------|
|         1 | SF    | Sulfate (SO₄)                        |  98.0         | 1500 (Tang mix)         | yes         |
|         2 | SS    | Sea salt (NaCl)                      |  58.5         | 1500 (Tang mix)         | yes         |
|         3 | ECOB  | Elemental carbon — hydrophobic       |  12.0         | 2200                    | no          |
|         4 | ECIL  | Elemental carbon — hydrophilic       |  12.0         | 2200                    | yes         |
|         5 | OCOB  | Organic carbon — hydrophobic         |  12.0         | 1400                    | no          |
|         6 | OCIL  | Organic carbon — hydrophilic         |  12.0         | 1400                    | yes         |
|         7 | DUST  | Mineral dust                          | 100.0         | 2650                    | no          |

Diagnostics — derived inside GCHP (NOT independent prognostic tracers,
though they may appear in SpeciesConcVV output depending on
configuration):

| SRT index | Code  | Long name                                          | Notes |
|----------:|-------|----------------------------------------------------|-------|
|         8 | NH4   | Ammonium                                            | Apportioned from bulk NH₃ to bins via sulfate mass fraction (capped at (NH₄)₂SO₄ stoichiometry, ratio 0.375). Not stored per-bin in default output — re-derive from SF if needed. |
|         9 | AW    | Aerosol water                                       | Equilibrium with ambient RH via Tang(1997) growth factors. Per-bin AW01..AW{IBINS} *is* in default SpeciesConcVV output. |

> **Density caveat.** GCHP's `INODENS` (tomas_mod.F90:8958-9199)
> computes a per-bin density for the inorganic mixture (SF + SS + NH₄ + AW)
> using Tang(1997) polynomial fits, defaulting to 1500 kg/m³ when the
> total inorganic mass is zero. EC, OC, and dust carry fixed densities.
> The table above lists the **default** density per species — for full
> fidelity, downstream Mie code should call a composition-aware density
> function rather than picking from this column.

## NetCDF variable naming

For each bin index `nn ∈ {01, 02, …, IBINS}`:

```
SpeciesConcVV_NK<nn>        → particle number  (1e3 · num · mol_dry⁻¹)
SpeciesConcVV_<CODE><nn>    → species mass mixing ratio  (mol/mol dry)
```

The vSmartMOM reader probes `SpeciesConcVV_NK01` to detect TOMAS presence
and discovers per-species availability by checking the `<CODE>` prefixes
listed in the species table above. AW is included when present in the
output; NH4 is skipped (derived).

## Unit conversions (per layer)

vSmartMOM converts GCHP raw output to per-layer (`#/cm³`, `μg/m³`) using
the layer-mean `Met_AD` (dry-air mass, kg) and `Met_AIRVOL` (grid-box
volume, m³). Following GCHP's convention:

```
air          = Met_AD / 28.9644e-3                     # mol dry air
N[k]         = NK[k] · air / (Met_AIRVOL · 1e6) / 1000 # #/cm³
M[k, spc]    = MOLWT[spc] · VV[k, spc] · air · 1e6 / Met_AIRVOL  # μg/m³
dN/dlogDp[k] = N[k] / (log10(Dp[k+1]) - log10(Dp[k]))
```

where `MOLWT[spc]` is in g/mol per the table above; the `× 1e6` factor
converts grams to micrograms. (vSmartMOM's `TOMASScheme.molar_masses`
internally stores kg/mol — the YAML loader converts on read — so the
`read_aerosol_cell` implementation uses `× 1e9` instead, which is the
same total scaling expressed in the kg/mol convention.)

## How vSmartMOM consumes this

```julia
# Open a GCHP file and read one cell.
GCHPFile(path; aerosol_scheme=:auto) do f
    scene = scene_at(f, idx, idy, idf)
    # scene.aerosols isa SectionalAerosolData{Float64, TOMASScheme{Float64}}
    # (when f.has_tomas == true)
end
```

Per-layer composition for Mie:

```julia
composition = bin_composition(scene.aerosols, ibin, ilev; basis=:volume)
# → SpeciesComposition{FT}(species, fractions, :volume)
n_eff = effective_ri(composition, ri_database, λ,
                     VolumeWeightedMixing(), scene.aerosols.scheme)
```

Generic accessors (work on any `<: AbstractAerosolBinData`):

| Accessor                            | Returns                                     |
|-------------------------------------|---------------------------------------------|
| `nbins(data)`                       | Number of size bins                         |
| `nlayers(data)`                     | Number of atmospheric layers                |
| `data.size_grid`                    | `AerosolSizeGrid` (centers, edges, units)   |
| `data.number`                       | `Matrix{FT}` of `#/cm³`, `(nbins, nlayers)` |
| `data.species`                      | `Vector{String}` of species codes           |
| `data.species_mass`                 | `Array{FT,3}` of μg/m³, `(nbins, nlayers, nspecies)` |
| `data.scheme`                       | The originating `<: AbstractAerosolScheme`  |

This keeps TOMAS-specific knowledge inside the TOMAS scheme; downstream
RT / Mie code never hard-codes `15` or species names.

## Bin integration for Mie

TOMAS bins are wide enough that evaluating Mie only at the bin center is
not generally adequate. vSmartMOM therefore represents the bin-to-Mie
choice explicitly:

| Integration mode | Meaning |
|------------------|---------|
| `LogNormalFit()` | Fit one lognormal to the binned column. Kept as an explicit future strategy. |
| `ConstantIntegrationPerBin(nquad=160)` | Integrate inside each bin in log-size space with constant `dN/dlogD`; conserves the bin-integrated number. |
| `LinearIntegrationPerBin(nquad=160)` | Integrate inside each bin in log-size space with a conservative piecewise-linear reconstruction of `dN/dlogD`; the first/last bins taper to zero at the outer edge while conserving bin number. Default for TOMAS scene ingest. |
| `DirectBinSum()` | Compatibility/diagnostic mode using the bin center only; not recommended for production Mie over wide bins. |

For the AOD diagnostic, `nquad` is the number of Gauss-Legendre points per
bin in log-size space. Values around 100-200 are the intended production
range; smaller values are useful for unit tests and smoke checks.

## Validation notes

### The `M_air` factor between NK and species totals

A historical cross-check during the GCHP-TOMAS ingest work: per-bin particle
number derived from **NK** and from **summed species mass** (converted to
number via `M = ρ·V_particle`) agreed in *shape* but differed in
magnitude by `M_air ≈ 28.96 g/mol`. The reason:

- **NK formula** divides by `M_air` once (`n_air = Met_AD / M_air`).
- **Species formula** also includes that factor at the `mol/mol → mol`
  step, but the subsequent `M = ρ·V_particle` step removes the
  `M_air` dependence — leaving NK ≈ M_air × species-derived number.

The ingest code in `read_aerosol_cell` uses the NK form directly (per
GCHP convention) and converts the mass channel separately
into `μg/m³`, so the two channels are independent and the factor does
not appear in the produced [`SectionalAerosolData`](@ref). This
section is kept as developer reference because it surfaces immediately
when comparing the two derivations.

### Loading vSmartMOM output in Python

The benchmark writer produces NetCDF-4 files (HDF5 on disk).
`netCDF4`, `xarray`, or `h5py` all work:

```python
import xarray as xr
ds = xr.open_dataset("scene_02_070_021.nc4")
# axes: (sza, view_pair, brdf, pol, spec)
R = ds["R"].values                  # shape (N_sza, N_view, N_brdf, pol_n, nSpec)
sza_deg = ds["sza"].values          # 1-D
vza_deg, vaz_deg = ds["vza"], ds["vaz"]   # paired (same length)
print("scene lat/lon:", ds.attrs["scene_lat"], ds.attrs["scene_lon"])
```

## See also

- [`data/aerosols/tomas_species.yaml`](../../../../data/aerosols/tomas_species.yaml)
  — machine-readable canonical table
- [`src/Aerosols/abstract_types.jl`](../../../../src/Aerosols/abstract_types.jl)
  — generic sectional-aerosol abstraction
- [`src/Aerosols/schemes/tomas15.jl`](../../../../src/Aerosols/schemes/tomas15.jl)
  — concrete TOMAS reader (loads this YAML at module init)
