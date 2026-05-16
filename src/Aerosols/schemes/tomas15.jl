"""
    TOMAS aerosol scheme reader

# Physical model

TOMAS (TwO-Moment Aerosol Sectional) is a size-resolved microphysics scheme that
discretizes the aerosol size distribution into a configurable number of mass
bins. Each bin tracks particle number and mass of each chemical species
(sulfate, organic carbon, dust, sea salt, etc.).

This enables accurate representation of:
- Nucleation, condensation, coagulation, and deposition
- Size-dependent optical properties via Mie theory
- CCN activation and cloud-aerosol interactions

The scheme is commonly used in GEOS-Chem and similar CTMs.
"""

"""
    read_tomas(scheme::TOMASScheme, ...)
    read_tomas15(config::Dict, netcdf_file::String, FT=Float64)

Read TOMAS aerosol data from NetCDF file including NK (number) and species.
`read_tomas15` is a compatibility shim for older configs.

# Arguments
- `config::Dict`: YAML configuration dictionary
- `netcdf_file::String`: Path to NetCDF file
- `FT`: Floating point type (default: Float64)

# Returns
- `AerosolData{TOMASScheme{FT}}`: Container with:
  - NK: Total particle number per bin (#/cm³)
  - dN_dlogD: Size distribution (cm⁻³)
  - Species mass: Mass per bin per species (μg/m³)
  - Species fractions: Fractional composition per bin

# Data Structure
NK variables: SpeciesConcVV_NK01 through NK{n_bins} (1000 × particles/mol_air)
Species variables: SpeciesConcVV_[SPECIES][BIN] (mol/mol dry air)
Meteorology: Met_AD (kg), Met_AIRVOL (m³)

# NK Unit Convention
NK is stored as 1000 × (particles/mol_air), convert using:
  N (#/cm³) = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)
where M_air = 28.9644e-3 kg/mol

# Configuration
The config dict must specify:
- aerosol_scheme/species: Dict of species with properties
- aerosol_scheme/size_bins: Size bin configuration
- netcdf_mapping/concentration_pattern or per-species aerosol_scheme/species/*/nc_prefix
- processing_options/vertical_flip or aerosol_scheme/options/vertical_flip: Whether to flip BOA→TOA to TOA→BOA
"""
function _tomas_concentration_variable(config::Dict, species_name::AbstractString, ibin::Integer)
    bin = @sprintf("%02d", ibin)
    mapping = get(config, "netcdf_mapping", nothing)

    if mapping !== nothing && haskey(mapping, "concentration_pattern")
        pattern = mapping["concentration_pattern"]
        return replace(
            pattern,
            "{species}" => species_name,
            "{bin:02d}" => bin,
            "{bin}" => string(ibin),
        )
    end

    species_config = config["aerosol_scheme"]["species"][species_name]
    prefix = get(species_config, "nc_prefix", "SpeciesConcVV_$(species_name)")
    return string(prefix, bin)
end

# ============================================================================
# TOMAS scheme presets — single source of truth: data/aerosols/tomas_species.yaml
# Verified against GCHP tomas_mod.F90:107-115 (species), 6497-6503 (M0),
# 6676-6688 (bin rule), 6694-6717 (MOLWT), 8924 (densities).
# ============================================================================

# Resolve the path to the YAML data file shipped with the package.
const _TOMAS_DATA_YAML = joinpath(@__DIR__, "..", "..", "..",
                                   "data", "aerosols", "tomas_species.yaml")

# Lazily loaded, cached. Read once at first use.
const _TOMAS_DATA = Ref{Union{Nothing, Dict}}(nothing)
function _tomas_data()
    if _TOMAS_DATA[] === nothing
        _TOMAS_DATA[] = YAML.load_file(_TOMAS_DATA_YAML)
    end
    return _TOMAS_DATA[]::Dict
end

# Bin variants exposed by the GCHP TOMAS family. The YAML key indexes
# `_tomas_data()["bin_schemes"]`.
const _TOMAS_VARIANTS = (:tomas12, :tomas15, :tomas30, :tomas40)

const _TOMAS_DEFAULT_RI_MAP = Dict{String, String}(
    "SF"   => "sulfate_suso",
    "SS"   => "seasalt_sscm",
    "ECOB" => "black_carbon",
    "ECIL" => "black_carbon",
    "OCOB" => "organic_carbon",
    "OCIL" => "organic_carbon",
    "DUST" => "dust_opac",
    "NH4"  => "sulfate_suso",
    "AW"   => "water",
)

"""
    TOMASScheme{FT} <: AbstractAerosolScheme{FT}

TOMAS sectional aerosol scheme metadata. The `variant` field selects the bin
count and bin-mass rule (`:tomas12`, `:tomas15`, `:tomas30`, or `:tomas40`).

`TOMAS15Scheme` is kept as a compatibility alias for code and configs that
still use the older fixed-bin name.
"""
struct TOMASScheme{FT} <: AbstractAerosolScheme{FT}
    variant::Symbol
    size_grid::AerosolSizeGrid{FT}
    species::Vector{String}
    n_bins::Int
    diam_min::FT
    diam_max::FT
    bin_edges::Vector{FT}
    bin_centers::Vector{FT}
    refractive_indices::Dict{String, String}
    ri_map::Dict{String, String}
    densities::Dict{String, FT}
    molar_masses::Dict{String, FT}
    bin_rule::Symbol
end

"""
    TOMAS15Scheme

Compatibility alias for [`TOMASScheme`](@ref). Prefer `TOMASScheme` in new
code.
"""
const TOMAS15Scheme = TOMASScheme

function _tomas_size_grid_from_edges(edges_nm::AbstractVector{FT},
                                     spacing::Symbol) where {FT}
    centers_nm = sqrt.(edges_nm[1:end - 1] .* edges_nm[2:end])
    edges_μm = collect(edges_nm) ./ FT(1000)
    centers_μm = centers_nm ./ FT(1000)
    widths = log10.(edges_μm[2:end]) .- log10.(edges_μm[1:end - 1])
    return AerosolSizeGrid{FT}(edges_μm, centers_μm, widths,
                               :diameter, spacing, "μm")
end

function _tomas_ri_map(species::AbstractVector{<:AbstractString},
                       overrides::Dict{String, String})
    out = Dict{String, String}()
    for sp in species
        code = String(sp)
        out[code] = get(overrides, code, get(_TOMAS_DEFAULT_RI_MAP, code, code))
    end
    return out
end

"""
    _build_xk(variant::Symbol, FT) -> Vector{FT}

Mass bin edges (kg/particle) following the GCHP `INIT_TOMAS` rule for the
given variant. `length(Xk) = n_bins + 1`.
"""
function _build_xk(variant::Symbol, ::Type{FT}) where {FT}
    data = _tomas_data()["bin_schemes"][string(variant)]
    n    = data["n_bins"]::Int
    Mo   = FT(data["M0_kg"])
    rule = Symbol(data["rule"])

    Xk = zeros(FT, n + 1)
    if rule === :quadruple_then_x32
        # tomas_mod.F90:6677-6683 — quadruple for k < IBINS, then ×32
        for k in 1:(n + 1)
            Xk[k] = k < n ? Mo * FT(4)^(k - 1) : Xk[k - 1] * FT(32)
        end
    elseif rule === :double
        # tomas_mod.F90:6685-6687 — mass-doubling for all bins
        for k in 1:(n + 1)
            Xk[k] = Mo * FT(2)^(k - 1)
        end
    else
        throw(ArgumentError("Unknown TOMAS bin rule: $rule (expected :quadruple_then_x32 or :double)"))
    end
    return Xk
end

"""
    _species_table() -> Vector{NamedTuple}

Parsed species rows (canonical GCHP order). Each row carries `code`,
`gchp_var_prefix`, `molwt_g_per_mol`, `density_kg_per_m3`, `component_type`,
`hydrophilic`.
"""
function _species_table()
    rows = _tomas_data()["species"]::Vector
    return [(; code = String(r["code"]),
              gchp_var_prefix = r["gchp_var_prefix"] === nothing ? nothing : String(r["gchp_var_prefix"]),
              molwt_g_per_mol = Float64(r["molwt_g_per_mol"]),
              density_kg_per_m3 = Float64(r["density_kg_per_m3"]),
              component_type = Symbol(r["component_type"]),
              hydrophilic = Bool(r["hydrophilic"]))
            for r in rows]
end

"""
    TOMASScheme(variant::Symbol=:tomas15; FT=Float64,
                include_species=:default,
                refractive_index_keys=Dict{String,String}())
        -> TOMASScheme{FT}

Construct a TOMAS scheme from a named variant. TOMAS12 / TOMAS15 / TOMAS30 /
TOMAS40 are all supported. The default is `:tomas15`.

Static metadata (molar masses, densities, species order) comes from
[`data/aerosols/tomas_species.yaml`](@ref) — the single source of truth
mirroring GCHP upstream.

# Arguments
- `variant`: one of `:tomas12, :tomas15, :tomas30, :tomas40`.
- `FT`: float type for stored numbers.
- `include_species`: `:default` keeps only mass-tracer species + AW (skips
  NH4 which GCHP derives from bulk NH3). Pass `:all` to include everything,
  or a `Vector{String}` to pick explicitly.
- `refractive_index_keys`: optional override mapping species code → key in
  a `RefractiveIndexDatabase`. Defaults map TOMAS codes to the packaged RI
  database keys used by the AOD diagnostic.

# Sensible defaults
- Diameters reported here are computed from mass bin edges (`Xk`) using the
  scheme-default density of 1500 kg/m³ as a placeholder — the actual
  density used in downstream Mie should come from the species composition
  of each bin (via [`bin_composition`](@ref) + [`densities`](@ref)).
"""
function TOMASScheme(variant::Symbol = :tomas15;
                     FT::DataType = Float64,
                     include_species::Union{Symbol, AbstractVector{<:AbstractString}} = :default,
                     refractive_index_keys::Dict{String, String} = Dict{String, String}())
    variant ∈ _TOMAS_VARIANTS || throw(ArgumentError(
        "Unknown TOMAS variant :$variant (expected one of $(_TOMAS_VARIANTS))"))

    rows     = _species_table()
    scheme_data = _tomas_data()["bin_schemes"][string(variant)]
    n_bins   = scheme_data["n_bins"]::Int
    bin_rule = Symbol(scheme_data["rule"])

    # Filter species per the include_species policy.
    selected = if include_species === :default
        # Mass-tracers + AW. Skip NH4 (derived from bulk NH3 inside GCHP).
        filter(r -> r.component_type === :mass || r.code == "AW", rows)
    elseif include_species === :all
        rows
    else
        codes = Set(String.(include_species))
        filter(r -> r.code in codes, rows)
    end

    species          = String[r.code for r in selected]
    densities_dict   = Dict{String, FT}(r.code => FT(r.density_kg_per_m3) for r in selected)
    # IMPORTANT: `TOMASScheme.molar_masses` is documented (and used by the
    # legacy `read_tomas15` reader) in **kg/mol**. The YAML stores values in
    # GCHP-native **g/mol**; convert at load time so the struct invariant
    # holds for callers that mix legacy + new code paths.
    molar_masses     = Dict{String, FT}(r.code => FT(r.molwt_g_per_mol) / FT(1000) for r in selected)
    ri_map           = _tomas_ri_map(species, refractive_index_keys)

    # Bin geometry — diameters in nm to match the legacy struct's units.
    Xk_kg   = _build_xk(variant, FT)
    ρ_ref   = FT(1500)                                              # placeholder
    Dp_m    = (Xk_kg ./ ρ_ref .* FT(6) ./ FT(π)) .^ (FT(1) / FT(3))
    Dp_nm   = Dp_m .* FT(1e9)
    bin_edges   = Dp_nm
    bin_centers = sqrt.(bin_edges[1:end - 1] .* bin_edges[2:end])
    spacing = bin_rule === :double ? :mass_double : :mass_quadruple
    size_grid = _tomas_size_grid_from_edges(bin_edges, spacing)

    return TOMASScheme{FT}(variant, size_grid, species, n_bins,
                           FT(bin_edges[1]), FT(bin_edges[end]),
                           bin_edges, bin_centers,
                           ri_map, ri_map, densities_dict, molar_masses,
                           bin_rule)
end

"""
    TOMASScheme(config::Dict, FT=Float64)

Construct a TOMAS scheme from the legacy YAML configuration format. This
preserves the old logarithmic diameter-bin behavior for `read_aerosol_data`
fixtures while new GCHP scene ingest uses `TOMASScheme(:tomas15)`.
"""
function TOMASScheme(config::Dict, FT=Float64)
    species_config = config["aerosol_scheme"]["species"]
    size_config = config["aerosol_scheme"]["size_bins"]

    species = collect(String, keys(species_config))
    n_bins = size_config["n_bins"]
    diam_min = FT(size_config["diam_min_nm"])
    diam_max = FT(size_config["diam_max_nm"])

    bin_edges = diam_min .* (diam_max / diam_min) .^
                (FT.(collect(0:n_bins)) ./ FT(n_bins))
    bin_centers = sqrt.(bin_edges[1:end - 1] .* bin_edges[2:end])
    size_grid = _tomas_size_grid_from_edges(bin_edges, :log10)

    refractive_indices = Dict{String, String}()
    densities_dict = Dict{String, FT}()
    molar_masses = Dict{String, FT}()

    for (sp, sp_config) in species_config
        code = String(sp)
        refractive_indices[code] = sp_config["refractive_index"]
        densities_dict[code] = FT(sp_config["density"])
        molar_masses[code] = FT(sp_config["molar_mass"])
    end

    return TOMASScheme{FT}(:legacy_config, size_grid, species, n_bins,
                           diam_min, diam_max, bin_edges, bin_centers,
                           refractive_indices, refractive_indices,
                           densities_dict, molar_masses, :log10)
end

# ----------------------------------------------------------------------------
# Generic-API methods on TOMASScheme
# ----------------------------------------------------------------------------

# Per-species density vector aligned with `species_codes`.
function densities(scheme::TOMASScheme{FT},
                   species_codes::AbstractVector{<:AbstractString}) where {FT}
    return FT[scheme.densities[String(s)] for s in species_codes]
end

"""
    refractive_index_key(scheme, species_code) -> String

Resolve an aerosol scheme's species code to a key in
`RefractiveIndexDatabase`.
"""
refractive_index_key(scheme::TOMASScheme, species_code::AbstractString) =
    get(scheme.ri_map, String(species_code), String(species_code))

# ----------------------------------------------------------------------------
# Per-cell reader returning `SectionalAerosolData{FT, TOMASScheme{FT}}`
# ----------------------------------------------------------------------------

const _M_AIR_KG_PER_MOL = 28.9644e-3

"""
    read_aerosol_cell(scheme, ds, idx, idy, idf; ...) -> SectionalAerosolData

Slice one cell `(idx, idy, idf)` from an open `NCDataset` carrying TOMAS
variables (`SpeciesConcVV_NK01..NK{IBINS}` and
`SpeciesConcVV_<CODE>01..{IBINS}` for each per-bin mass species).
Reads `Met_AD` (kg dry air/layer) and `Met_AIRVOL` (m³/layer) for unit
conversion. Returns column-shaped data oriented TOA→BOA.

Unit conventions (GCHP convention, verified against tomas_mod.F90):
- `Nk[bin, lev]      = NK[bin, lev] · air / (Met_AIRVOL · 1e6) / 1000`     [#/cm³]
- `Mk[bin, lev, spc] = MOLWT[spc] · VV[bin, lev, spc] · air · 1e9 / Met_AIRVOL`  [μg/m³]
- `air = Met_AD / 28.9644e-3` [mol dry air]

MOLWT is stored on the scheme in **kg/mol**; the `× 1e9` factor converts
kg/mol to μg/mol for the mass-concentration calculation.
"""
function read_aerosol_cell(scheme::TOMASScheme{FT}, ds,
                           idx::Integer, idy::Integer, idf::Integer;
                           ri_database::Union{Nothing, RefractiveIndexDatabase{FT}} = nothing,
                           mixing_rule::AbstractMixingRule = VolumeWeightedMixing(),
                           integration::AbstractBinIntegration = LinearIntegrationPerBin()) where {FT}
    nlev    = length(ds["lev"])
    n_bins  = scheme.n_bins
    species = scheme.species
    nspc    = length(species)

    Met_AD     = FT.(ds["Met_AD"].var[idx, idy, idf, :, 1])
    Met_AIRVOL = FT.(ds["Met_AIRVOL"].var[idx, idy, idf, :, 1])
    air        = Met_AD ./ FT(_M_AIR_KG_PER_MOL)

    # Number per bin (#/cm³), still BOA→TOA
    N = zeros(FT, n_bins, nlev)
    for ibin in 1:n_bins
        var_name = @sprintf("SpeciesConcVV_NK%02d", ibin)
        haskey(ds, var_name) || continue
        nk = FT.(ds[var_name].var[idx, idy, idf, :, 1])
        @inbounds for ilev in 1:nlev
            vol_cm3 = Met_AIRVOL[ilev] * FT(1e6)
            N[ibin, ilev] = (nk[ilev] / FT(1000)) * air[ilev] / vol_cm3
        end
    end

    # Per-species mass (μg/m³), still BOA→TOA
    mass = zeros(FT, n_bins, nlev, nspc)
    for (ispc, spc) in enumerate(species)
        # GCHP variable prefix matches the species code in standard output.
        prefix = "SpeciesConcVV_" * spc
        # `scheme.molar_masses[spc]` is in kg/mol (struct invariant — see
        # `TOMASScheme(variant::Symbol)`). Convert kg/mol → μg/mol via ×1e9.
        MW_ug_per_mol = scheme.molar_masses[spc] * FT(1e9)
        for ibin in 1:n_bins
            var_name = string(prefix, @sprintf("%02d", ibin))
            haskey(ds, var_name) || continue
            vv = FT.(ds[var_name].var[idx, idy, idf, :, 1])
            @inbounds for ilev in 1:nlev
                mol_spc = vv[ilev] * air[ilev]
                mass[ibin, ilev, ispc] = (mol_spc * MW_ug_per_mol) / Met_AIRVOL[ilev]
            end
        end
    end

    # Reverse to TOA→BOA along the layer axis
    N    = reverse(N;    dims = 2)
    mass = reverse(mass; dims = 2)

    return SectionalAerosolData{FT, TOMASScheme{FT}}(
        scheme, scheme.size_grid, N, copy(species), mass,
        (; variant = scheme.variant, n_bins = n_bins, bin_rule = scheme.bin_rule),
        ri_database, mixing_rule, integration)
end

"""
    read_tomas(scheme::TOMASScheme, config::Dict, netcdf_file::String)
    read_tomas(scheme::TOMASScheme, netcdf_file::String; vertical_flip=true)

Read TOMAS aerosol data with explicit scheme dispatch. The config form keeps
legacy YAML variable-mapping support; the direct form uses standard GCHP
`SpeciesConcVV_<CODE><BIN>` names.
"""
function read_tomas(scheme::TOMASScheme{FT}, config::Dict,
                    netcdf_file::String) where {FT}
    # Open NetCDF file
    ds = NCDataset(netcdf_file)
    
    try
        # Extract coordinates
        coordinates = extract_coordinates(ds, config)
        n_levels = length(coordinates["lev"])
        
        # Get processing options
        vertical_flip = get(aerosol_processing_options(config), "vertical_flip", false)
        
        # Physical constants
        M_air = FT(28.9644e-3)  # kg/mol (molar mass of dry air)
        
        # ====================================================================
        # Read NK variables (particle number)
        # ====================================================================
        NK_raw = zeros(FT, scheme.n_bins, n_levels)  # Raw NK values
        N_concentration = zeros(FT, scheme.n_bins, n_levels)  # #/cm³
        
        # Read meteorology for conversion
        Met_AD_full = Array(ds["Met_AD"])  # kg
        Met_AIRVOL_full = Array(ds["Met_AIRVOL"])  # m³
        
        # Average over horizontal dimensions
        # Shape: (nf=6, Xdim=24, Ydim=24, lev=72, time=1) → (lev,)
        Met_AD_profile = dropdims(
            mean(Met_AD_full, dims=(1, 2, 3)),
            dims=(1, 2, 3)
        )[:, 1]
        
        Met_AIRVOL_profile = dropdims(
            mean(Met_AIRVOL_full, dims=(1, 2, 3)),
            dims=(1, 2, 3)
        )[:, 1]
        
        # Read NK for all bins
        for ibin in 1:scheme.n_bins
            var_name = @sprintf("SpeciesConcVV_NK%02d", ibin)
            
            if haskey(ds, var_name)
                data_full = Array(ds[var_name])
                
                # Average over horizontal
                data_profile = dropdims(
                    mean(data_full, dims=(1, 2, 3)),
                    dims=(1, 2, 3)
                )[:, 1]
                
                NK_raw[ibin, :] = data_profile
            else
                @warn "Variable $var_name not found, setting to zero"
                NK_raw[ibin, :] .= 0.0
            end
        end
        
        # Convert NK to concentration (#/cm³)
        # Formula: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)
        for ilev in 1:n_levels
            n_air = Met_AD_profile[ilev] / M_air  # mol
            vol_cm3 = Met_AIRVOL_profile[ilev] * FT(1e6)  # cm³
            
            for ibin in 1:scheme.n_bins
                N_concentration[ibin, ilev] = (NK_raw[ibin, ilev] / FT(1000)) * 
                                             (n_air / vol_cm3)
            end
        end
        
        # Compute dN/dlogD
        delta_logD = log10.(scheme.bin_edges[2:end]) .- log10.(scheme.bin_edges[1:end-1])
        dN_dlogD = zeros(FT, scheme.n_bins, n_levels)
        
        for ilev in 1:n_levels
            for ibin in 1:scheme.n_bins
                dN_dlogD[ibin, ilev] = N_concentration[ibin, ilev] / FT(delta_logD[ibin])
            end
        end
        
        # Apply vertical flip if requested
        if vertical_flip
            NK_raw = reverse(NK_raw, dims=2)
            N_concentration = reverse(N_concentration, dims=2)
            dN_dlogD = reverse(dN_dlogD, dims=2)
            Met_AD_profile = reverse(Met_AD_profile)
            Met_AIRVOL_profile = reverse(Met_AIRVOL_profile)
        end
        
        # ====================================================================
        # Read mass species and compute particle number fractions
        # ====================================================================
        # ====================================================================
        # Read mass species and compute particle number fractions
        # ====================================================================
        
        species_data = Dict{String, AerosolSpeciesData}()
        
        # Total mass per bin (for computing fractions)
        total_mass = zeros(FT, scheme.n_bins, n_levels)  # μg/m³
        species_mass = Dict{String, Array{FT, 2}}()  # Species → (bins × levels)
        species_particle_number = Dict{String, Array{FT, 2}}()  # Species → (bins × levels)
        
        for species_name in scheme.species
            # Read all bins for this species
            concentrations = zeros(FT, scheme.n_bins, n_levels)  # mol/mol
            mass_conc = zeros(FT, scheme.n_bins, n_levels)  # μg/m³
            part_num = zeros(FT, scheme.n_bins, n_levels)  # #/cm³
            
            for ibin in 1:scheme.n_bins
                var_name = _tomas_concentration_variable(config, species_name, ibin)
                
                if haskey(ds, var_name)
                    data_full = Array(ds[var_name])
                    
                    # Average over horizontal
                    data_profile = dropdims(
                        mean(data_full, dims=(1, 2, 3)),
                        dims=(1, 2, 3)
                    )[:, 1]
                    
                    concentrations[ibin, :] = data_profile
                else
                    @warn "Variable $var_name not found for species $species_name, bin $ibin"
                    concentrations[ibin, :] .= 0.0
                end
            end
            
            # Convert species from mol/mol to mass concentration (μg/m³)
            # M = MW × (mol_species/mol_air) × (Met_AD/M_air) / Met_AIRVOL
            MW_ug_mol = scheme.molar_masses[species_name] * FT(1e9)  # kg/mol → μg/mol
            
            for ilev in 1:n_levels
                n_air = Met_AD_profile[ilev] / M_air  # mol
                vol_m3 = Met_AIRVOL_profile[ilev]  # m³
                
                for ibin in 1:scheme.n_bins
                    mol_species = concentrations[ibin, ilev] * n_air
                    mass_conc[ibin, ilev] = (mol_species * MW_ug_mol) / vol_m3
                end
            end
            
            # Convert mass to particle number using bin size and density
            # N = M / m_particle where m_particle = ρ × (4/3)π r³
            ρ = scheme.densities[species_name]  # kg/m³
            
            for ibin in 1:scheme.n_bins
                # Use bin center diameter
                D_nm = scheme.bin_centers[ibin]
                r_cm = (D_nm * FT(1e-7)) / FT(2)  # nm → cm
                V_particle = (FT(4)/FT(3)) * FT(π) * r_cm^3  # cm³
                m_particle_ug = ρ * V_particle * FT(1e9)  # μg
                
                for ilev in 1:n_levels
                    if m_particle_ug > 0
                        # N (#/m³) = M (μg/m³) / m_particle (μg)
                        N_m3 = mass_conc[ibin, ilev] / m_particle_ug
                        part_num[ibin, ilev] = N_m3 * FT(1e-6)  # #/cm³
                    else
                        part_num[ibin, ilev] = FT(0)
                    end
                end
            end
            
            # Apply vertical flip if requested
            if vertical_flip
                concentrations = reverse(concentrations, dims=2)
                mass_conc = reverse(mass_conc, dims=2)
                part_num = reverse(part_num, dims=2)
            end
            
            # Accumulate totals
            total_mass .+= mass_conc
            species_mass[species_name] = mass_conc
            species_particle_number[species_name] = part_num
            
            # Store species data
            data_dict = Dict{String, Any}(
                "concentration" => concentrations,  # mol/mol
                "mass" => mass_conc,  # μg/m³
                "particle_number" => part_num  # #/cm³
            )
            
            units_dict = Dict{String, String}(
                "concentration" => "mol mol-1 dry air",
                "mass" => "μg m-3",
                "particle_number" => "cm-3"
            )
            
            species_data[species_name] = AerosolSpeciesData(
                data_dict,
                units_dict,
                "TOMAS size-resolved $(species_name) aerosol"
            )
        end
        
        # ====================================================================
        # Compute species fractions per bin
        # ====================================================================
        
        # Fraction by particle number
        species_fractions = Dict{String, Array{FT, 2}}()
        
        for species_name in scheme.species
            fractions = zeros(FT, scheme.n_bins, n_levels)
            
            for ilev in 1:n_levels
                for ibin in 1:scheme.n_bins
                    if N_concentration[ibin, ilev] > 0
                        fractions[ibin, ilev] = species_particle_number[species_name][ibin, ilev] / 
                                               N_concentration[ibin, ilev]
                    else
                        fractions[ibin, ilev] = FT(0)
                    end
                end
            end
            
            species_fractions[species_name] = fractions
        end
        
        # ====================================================================
        # Store NK data
        # ====================================================================
        
        nk_data_dict = Dict{String, Any}(
            "NK_raw" => NK_raw,  # Raw NK values from file
            "concentration" => N_concentration,  # #/cm³
            "dN_dlogD" => dN_dlogD,  # cm⁻³
            "species_fractions" => species_fractions  # Fraction per species
        )
        
        nk_units_dict = Dict{String, String}(
            "NK_raw" => "1000 × (particles/mol_air)",
            "concentration" => "cm-3",
            "dN_dlogD" => "cm-3",
            "species_fractions" => "dimensionless"
        )
        
        species_data["NK"] = AerosolSpeciesData(
            nk_data_dict,
            nk_units_dict,
            "TOMAS total particle number distribution with species fractions"
        )
        
        # Extract metadata
        metadata = extract_metadata(ds)
        
        return AerosolData(scheme, species_data, coordinates, metadata)
        
    finally
        close(ds)
    end
end

function read_tomas(scheme::TOMASScheme, config::Dict,
                    netcdf_file::String, FT::DataType)
    scheme_ft = scheme isa TOMASScheme{FT} ? scheme : TOMASScheme(config, FT)
    return read_tomas(scheme_ft, config, netcdf_file)
end

function read_tomas(scheme::TOMASScheme{FT}, netcdf_file::String;
                    vertical_flip::Bool = true) where {FT}
    species_config = Dict{String, Any}()
    for sp in scheme.species
        species_config[sp] = Dict{String, Any}(
            "refractive_index" => refractive_index_key(scheme, sp),
            "density" => scheme.densities[sp],
            "molar_mass" => scheme.molar_masses[sp],
            "nc_prefix" => "SpeciesConcVV_$sp",
        )
    end
    config = Dict{String, Any}(
        "aerosol_scheme" => Dict{String, Any}(
            "type" => "TOMAS",
            "species" => species_config,
            "size_bins" => Dict{String, Any}(
                "n_bins" => scheme.n_bins,
                "diam_min_nm" => scheme.diam_min,
                "diam_max_nm" => scheme.diam_max,
            ),
        ),
        "processing_options" => Dict{String, Any}("vertical_flip" => vertical_flip),
    )
    return read_tomas(scheme, config, netcdf_file)
end

read_tomas15(config::Dict, netcdf_file::String, FT=Float64) =
    read_tomas(TOMASScheme(config, FT), config, netcdf_file)

"""
    compute_number_concentration(vmr::AbstractVector, pressure::AbstractVector, 
                                  temperature::AbstractVector)

Convert volume mixing ratio to number concentration.

# Arguments
- `vmr::AbstractVector`: Volume mixing ratio (mol/mol)
- `pressure::AbstractVector`: Pressure (Pa)
- `temperature::AbstractVector`: Temperature (K)

# Returns
- `Vector{Float64}`: Number concentration (#/cm³)

# Formula
n = VMR × (P / (k_B × T)) × 10^-6

where k_B = 1.380649e-23 J/K (Boltzmann constant)
"""
function compute_number_concentration(
    vmr::AbstractVector,
    pressure::AbstractVector,
    temperature::AbstractVector
)
    k_B = 1.380649e-23  # J/K
    
    # Number density of air (molecules/m³)
    n_air = pressure ./ (k_B .* temperature)
    
    # Number concentration (#/m³)
    n_aerosol = vmr .* n_air
    
    # Convert to #/cm³
    return n_aerosol .* 1e-6
end

"""
    compute_mass_concentration(vmr::AbstractVector, molar_mass::Float64,
                               pressure::AbstractVector, temperature::AbstractVector)

Convert volume mixing ratio to mass concentration.

# Arguments
- `vmr::AbstractVector`: Volume mixing ratio (mol/mol)
- `molar_mass::Float64`: Molar mass (kg/mol)
- `pressure::AbstractVector`: Pressure (Pa)
- `temperature::AbstractVector`: Temperature (K)

# Returns
- `Vector{Float64}`: Mass concentration (μg/m³)

# Formula
ρ = VMR × (P × M) / (R × T) × 10^9

where R = 8.314462618 J/(mol·K) (universal gas constant)
"""
function compute_mass_concentration(
    vmr::AbstractVector,
    molar_mass::Float64,
    pressure::AbstractVector,
    temperature::AbstractVector
)
    R = 8.314462618  # J/(mol·K)
    
    # Mass concentration (kg/m³)
    mass_conc = (vmr .* pressure .* molar_mass) ./ (R .* temperature)
    
    # Convert to μg/m³
    return mass_conc .* 1e9
end

"""
    bin_volume(diam_nm::Float64)

Calculate volume of a spherical particle from its diameter.

# Arguments
- `diam_nm::Float64`: Particle diameter (nm)

# Returns
- `Float64`: Particle volume (nm³)
"""
function bin_volume(diam_nm::Float64)
    radius_nm = diam_nm / 2.0
    return (4.0/3.0) * π * radius_nm^3
end
