# Aerosol Integration Design for vSmartMOM.jl

**Date**: October 15, 2025  
**Context**: Reading aerosol data from GEOSChem-TOMAS NetCDF files  
**Goal**: Flexible aerosol optical properties from binned size distributions

---

## 1. GEOSChem-TOMAS File Structure

### Discovered Structure
From `GEOSChem.Custom.20190702_0000z.nc4`:

**Dimensions**:
- `Xdim = 24`, `Ydim = 24` (horizontal grid per face)
- `nf = 6` (cubed-sphere faces)
- `lev = 72` (vertical levels, BOA to TOA)
- `time = 1`

**Aerosol Species** (15 bins each):
- `DUST01-15`: Mineral dust (15 size bins)
- `SS01-15`: Sea salt (15 size bins)
- `SF01-15`: Sulfate (15 size bins)
- `ECIL01-15`: Elemental carbon (hydrophilic)
- `ECOB01-15`: Elemental carbon (hydrophobic)
- `OCIL01-15`: Organic carbon (hydrophilic)
- `OCOB01-15`: Organic carbon (hydrophobic)
- `NK01-15`: Nitrate/potassium (?)
- `AW01-15`: Aerosol water (15 bins)

**Other Species**:
- Gas-phase: `CO2`, `CO`, `CH4`, `C2H6`, `H2O`, `N2O`
- Ionic: `NH4`, `NIT`, `NITs`

**Variable Format**:
```
float SpeciesConcVV_<SPECIES>(time, lev, nf, Ydim, Xdim)
  units: "mol mol-1 dry"  (volume mixing ratio)
```

### TOMAS Bin Definitions

TOMAS uses 15 logarithmically-spaced bins. Standard configuration:
- **Bin 1**: 10 - 20 nm diameter
- **Bin 2**: 20 - 40 nm
- **Bin 3**: 40 - 80 nm
- ...
- **Bin 15**: 5.12 - 10.24 μm

**Bin edges** (dry diameter in nm):
```julia
TOMAS_BIN_EDGES = [
    10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 1280.0,
    2560.0, 5120.0, 10240.0  # Continuing pattern...
]
```

For **radius** (μm), divide by 2 and convert nm→μm:
```julia
TOMAS_RADIUS_EDGES_UM = TOMAS_BIN_EDGES ./ 2000.0  # nm→μm, diameter→radius
```

---

## 2. Proposed Type System

### Core Types

```julia
# src/Scattering/aerosol_types.jl (extend existing types.jl)

"""
TOMAS bin configuration (15 bins, logarithmic spacing)
"""
struct TOMASSizeGrid
    n_bins::Int  # 15
    radius_edges::Vector{Float64}  # [μm] - 16 edges for 15 bins
    radius_centers::Vector{Float64}  # [μm] - 15 geometric means
end

"""
Create standard TOMAS-15 grid
"""
function standard_tomas15_grid()
    # Bin edges in dry diameter (nm)
    diam_edges_nm = [
        10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 1280.0,
        2560.0, 5120.0, 10240.0, 20480.0, 40960.0, 81920.0, 163840.0, 327680.0
    ]
    
    # Convert to radius (μm)
    radius_edges = diam_edges_nm ./ 2000.0  # nm→μm, diam→radius
    
    # Geometric mean radius for each bin
    radius_centers = sqrt.(radius_edges[1:end-1] .* radius_edges[2:end])
    
    return TOMASSizeGrid(15, radius_edges, radius_centers)
end

"""
Binned number density distribution (from model bins)
"""
struct BinnedNumberDistribution
    radius_grid::TOMASSizeGrid
    number_density::Vector{Float64}  # [#/m³] per bin
    
    # Constructor with validation
    function BinnedNumberDistribution(grid::TOMASSizeGrid, nd::Vector{Float64})
        @assert length(nd) == grid.n_bins "Number density must match grid bins"
        @assert all(nd .>= 0) "Number density must be non-negative"
        new(grid, nd)
    end
end

"""
Single aerosol species with optical properties
"""
struct AerosolSpecies
    name::String  # "dust", "seasalt", "sulfate", "BC_hydrophilic", etc.
    size_distribution::BinnedNumberDistribution
    refractive_index::ComplexF64  # n - ik (wavelength-dependent later)
    density::Float64  # [kg/m³] bulk density
end

"""
Mixed aerosol layer (multiple species coexisting)
"""
struct MixedAerosolLayer{FT<:AbstractFloat}
    species::Vector{AerosolSpecies}
    layer_top_pressure::FT  # [hPa]
    layer_bottom_pressure::FT  # [hPa]
    temperature::FT  # [K] layer mean
    metadata::Dict{String, Any}
end

"""
Full vertical aerosol profile
"""
struct AerosolProfile{FT<:AbstractFloat}
    layers::Vector{MixedAerosolLayer{FT}}
    reference_wavelength::FT  # [nm] for diagnostics
    source_file::String
    grid_location::Tuple{Int, Int, Int}  # (idx, idy, idf)
end
```

---

## 3. Reader Implementation

### File Structure

```
src/IO/NetCDF/
├── GeosChem.jl          # Existing atmospheric+gas reader
└── GeosChemAerosols.jl  # NEW: Aerosol-specific reader
```

### Core Reader Function

```julia
# src/IO/NetCDF/GeosChemAerosols.jl

"""
Read aerosol data from GEOSChem-TOMAS file
"""
function read_geoschem_aerosols(
    src::GeosChemSource;
    species_to_read::Vector{String} = ["DUST", "SS", "SF", "ECIL", "OCOB"],
    FT::Type{<:AbstractFloat} = Float64
) :: AerosolProfile{FT}
    
    ds = NCDataset(src.path)
    idx, idy, idf = src.idx, src.idy, src.idf
    
    # Get vertical structure
    n_lev = ds.dim["lev"]
    pressure_grid = compute_pressure_levels(ds, idx, idy, idf)  # From existing code
    temperature_profile = reverse(ds["Met_T"][idf, idy, idx, :, 1])
    
    # Standard TOMAS grid
    tomas_grid = standard_tomas15_grid()
    
    # Read species concentrations and build layers
    layers = MixedAerosolLayer{FT}[]
    
    for iz in 1:n_lev
        layer_species = AerosolSpecies[]
        
        for species_type in species_to_read
            # Read all 15 bins for this species type
            bin_concentrations = zeros(FT, 15)
            
            for bin_idx in 1:15
                var_name = "SpeciesConcVV_$(species_type)$(lpad(bin_idx, 2, '0'))"
                
                if haskey(ds, var_name)
                    # VMR [mol/mol dry] - flip BOA→TOA to TOA→BOA
                    vmr = ds[var_name][idf, idy, idx, n_lev - iz + 1, 1]
                    bin_concentrations[bin_idx] = FT(vmr)
                end
            end
            
            # Convert VMR to number density [#/m³]
            number_density = vmr_to_number_density(
                bin_concentrations,
                pressure_grid[iz:iz+1],
                temperature_profile[iz],
                tomas_grid
            )
            
            # Skip if negligible
            if maximum(number_density) < 1e-10
                continue
            end
            
            # Create distribution
            dist = BinnedNumberDistribution(tomas_grid, number_density)
            
            # Get optical properties for this species
            n_complex, rho = get_species_properties(species_type)
            
            # Create aerosol species
            aerosol_spec = AerosolSpecies(
                species_type,
                dist,
                n_complex,
                rho
            )
            
            push!(layer_species, aerosol_spec)
        end
        
        # Create layer if it has aerosols
        if !isempty(layer_species)
            layer = MixedAerosolLayer{FT}(
                layer_species,
                FT(pressure_grid[iz]),
                FT(pressure_grid[iz+1]),
                FT(temperature_profile[iz]),
                Dict{String, Any}()
            )
            push!(layers, layer)
        end
    end
    
    close(ds)
    
    return AerosolProfile{FT}(
        layers,
        FT(550.0),  # Reference wavelength
        src.path,
        (idx, idy, idf)
    )
end

"""
Convert volume mixing ratio to number density

# Arguments
- `vmr`: Volume mixing ratios per bin [mol/mol]
- `p_levels`: Pressure at layer boundaries [hPa]
- `T`: Temperature [K]
- `grid`: TOMAS size grid

# Returns
- Number density per bin [#/m³]
"""
function vmr_to_number_density(
    vmr::Vector{FT},
    p_levels::Vector{FT},
    T::FT,
    grid::TOMASSizeGrid
) where {FT<:AbstractFloat}
    
    # Layer mean pressure [Pa]
    p_mean = mean(p_levels) * 100.0  # hPa → Pa
    
    # Air number density [molecules/m³] from ideal gas
    # n_air = p / (k_B * T)
    k_B = 1.380649e-23  # J/K (Boltzmann)
    n_air = p_mean / (k_B * T)
    
    # Number density per bin [#/m³]
    # Assumes 1:1 molar ratio (molecules per particle varies by size)
    # This is approximate - proper conversion needs particle composition
    number_density = vmr .* n_air
    
    # TODO: More sophisticated conversion accounting for:
    #   - Particle mass per bin
    #   - Species molecular weight
    #   - Hydration state (for AW bins)
    
    return number_density
end

"""
Get optical properties for aerosol species

# Returns
- `(n_complex, density)`: Refractive index (n-ik) and bulk density [kg/m³]
"""
function get_species_properties(species_type::String)
    # Based on OPAC, AERONET, literature values
    properties = Dict(
        "DUST" => (ComplexF64(1.53, 0.008), 2650.0),    # Mineral dust
        "SS" => (ComplexF64(1.50, 1e-8), 2200.0),       # Sea salt (dry)
        "SF" => (ComplexF64(1.43, 1e-8), 1770.0),       # Sulfate (H2SO4)
        "ECIL" => (ComplexF64(1.95, 0.79), 1800.0),     # BC hydrophilic
        "ECOB" => (ComplexF64(1.95, 0.79), 1800.0),     # BC hydrophobic
        "OCIL" => (ComplexF64(1.53, 0.006), 1400.0),    # OC hydrophilic
        "OCOB" => (ComplexF64(1.53, 0.006), 1400.0),    # OC hydrophobic
        "NK" => (ComplexF64(1.45, 1e-7), 1720.0),       # Nitrate
        "AW" => (ComplexF64(1.33, 0.0), 1000.0)         # Water
    )
    
    return get(properties, species_type, (ComplexF64(1.5, 0.01), 1500.0))
end
```

---

## 4. Optical Properties Computation

### Strategy: Direct Mie for Each Bin

```julia
# src/Scattering/binned_mie.jl (NEW)

"""
Compute Mie scattering for binned size distribution

# Arguments
- `species::AerosolSpecies`: Species with binned distribution
- `wavelengths::Vector{FT}`: Wavelengths [nm]
- `temperature::FT`: Layer temperature [K]

# Returns
- `(τ, ϖ, greek)`: Optical depth, single-scattering albedo, Greek coefficients
"""
function compute_binned_mie(
    species::AerosolSpecies,
    wavelengths::Vector{FT},
    temperature::FT = 273.15
) where {FT<:AbstractFloat}
    
    dist = species.size_distribution
    n_bins = dist.radius_grid.n_bins
    nλ = length(wavelengths)
    
    # Preallocate
    τ_total = zeros(FT, nλ)
    scat_total = zeros(FT, nλ)
    greek_total = nothing  # Accumulate weighted Greek coefficients
    
    # Loop over each size bin
    for i_bin in 1:n_bins
        r_center = dist.radius_grid.radius_centers[i_bin]  # [μm]
        n_dens = dist.number_density[i_bin]  # [#/m³]
        
        if n_dens < 1e-10
            continue  # Skip empty bins
        end
        
        # Compute Mie properties at bin center radius
        # (Use existing vSmartMOM Mie code)
        τ_bin, ϖ_bin, greek_bin = compute_mie_monodisperse(
            r_center,
            species.refractive_index,
            wavelengths,
            n_dens,
            temperature
        )
        
        # Accumulate
        τ_total .+= τ_bin
        scat_total .+= τ_bin .* ϖ_bin
        
        # Weight Greek coefficients by scattering contribution
        if isnothing(greek_total)
            greek_total = greek_bin .* (τ_bin .* ϖ_bin)
        else
            greek_total .+= greek_bin .* (τ_bin .* ϖ_bin)
        end
    end
    
    # Effective single-scattering albedo
    ϖ_eff = scat_total ./ (τ_total .+ 1e-30)
    
    # Normalize Greek coefficients
    greek_total ./= (scat_total .+ 1e-30)
    
    return (τ_total, ϖ_eff, greek_total)
end

"""
Compute layer optical properties from mixed aerosol layer
"""
function compute_layer_optics(
    layer::MixedAerosolLayer{FT},
    wavelengths::Vector{FT}
) where {FT<:AbstractFloat}
    
    nλ = length(wavelengths)
    τ_layer = zeros(FT, nλ)
    scat_layer = zeros(FT, nλ)
    greek_layer = nothing
    
    # Sum contributions from all species
    for species in layer.species
        τ_spec, ϖ_spec, greek_spec = compute_binned_mie(
            species,
            wavelengths,
            layer.temperature
        )
        
        τ_layer .+= τ_spec
        scat_layer .+= τ_spec .* ϖ_spec
        
        if isnothing(greek_layer)
            greek_layer = greek_spec .* (τ_spec .* ϖ_spec)
        else
            greek_layer .+= greek_spec .* (τ_spec .* ϖ_spec)
        end
    end
    
    # Effective properties
    ϖ_layer = scat_layer ./ (τ_layer .+ 1e-30)
    greek_layer ./= (scat_layer .+ 1e-30)
    
    return (τ=τ_layer, ϖ=ϖ_layer, greek=greek_layer)
end
```

---

## 5. Integration with vSmartMOM Workflow

### Extend geoschem_to_dict

```julia
# src/IO/NetCDF/GeosChem.jl (extend existing function)

function geoschem_to_dict(
    src::GeosChemSource;
    include_aerosols::Bool = true,
    aerosol_species::Vector{String} = ["DUST", "SS", "SF", "ECIL", "OCOB"]
)
    # ... existing atmospheric profile and gas absorption code ...
    
    config = Dict{String, Any}()
    config["atmospheric_profile"] = # ... existing ...
    config["absorption"] = # ... existing ...
    
    # NEW: Add aerosol scattering
    if include_aerosols
        aerosol_profile = read_geoschem_aerosols(
            src,
            species_to_read = aerosol_species
        )
        
        config["scattering"] = Dict{String, Any}(
            "aerosol_profile" => aerosol_profile,
            "source" => "geoschem_tomas",
            "λ_ref" => 550.0,
            "decomp_type" => "NAI2()"
        )
    end
    
    return config
end
```

### Modify model_from_parameters

```julia
# src/CoreRT/tools/model_from_parameters.jl (extend)

function model_from_parameters(params::vSmartMOM_Parameters)
    # ... existing code ...
    
    # NEW: Handle aerosol profile if present
    if hasfield(typeof(params), :aerosol_profile)
        aerosol_optics = compute_aerosol_optics_from_profile(
            params.aerosol_profile,
            params.spec_bands
        )
    else
        # Existing aerosol handling
        aerosol_optics = # ... existing code ...
    end
    
    # ... rest of model construction ...
end
```

---

## 6. Implementation Roadmap

### Phase 1: Core Infrastructure (Days 1-3)
- [ ] Create `aerosol_types.jl` with new types
- [ ] Implement `standard_tomas15_grid()`
- [ ] Add `BinnedNumberDistribution` constructor
- [ ] Unit tests for type creation

### Phase 2: Reader (Days 4-7)
- [ ] Create `GeosChemAerosols.jl`
- [ ] Implement `read_geoschem_aerosols()`
- [ ] Implement `vmr_to_number_density()` conversion
- [ ] Add `get_species_properties()` lookup table
- [ ] Test with real NetCDF file

### Phase 3: Optical Properties (Days 8-12)
- [ ] Create `binned_mie.jl`
- [ ] Implement `compute_binned_mie()`
- [ ] Extend existing Mie code for binned inputs
- [ ] Implement `compute_layer_optics()`
- [ ] Validate against known benchmarks

### Phase 4: Integration (Days 13-15)
- [ ] Extend `geoschem_to_dict()` for aerosols
- [ ] Modify `parameters_from_dict()` to handle aerosol profiles
- [ ] Update `model_from_parameters()`
- [ ] Create example: `examples/aerosol_rt_geoschem.jl`
- [ ] Documentation

### Phase 5: Advanced Features (Optional, Days 16+)
- [ ] Parametric distribution fitting (LogNormal to bins)
- [ ] Hygroscopic growth (use AW bins to adjust radii)
- [ ] Wavelength-dependent refractive indices
- [ ] Internal vs external mixing options
- [ ] AOD diagnostics and validation

---

## 7. Example Usage

```julia
using vSmartMOM

# Read atmospheric + aerosol data from GEOSChem
src = GeosChemSource("GEOSChem.Custom.20190702_0000z.nc4", 10, 10, 1)

# Get configuration with aerosols
config = geoschem_to_dict(
    src,
    include_aerosols = true,
    aerosol_species = ["DUST", "SS", "SF", "ECIL"]  # Select species
)

# Inspect aerosol profile
aerosol_prof = config["scattering"]["aerosol_profile"]
println("Number of layers with aerosols: $(length(aerosol_prof.layers))")

# Build model and run RT
params = parameters_from_dict(config)
model = model_from_parameters(params)
R, T = rt_run(model)

# Analyze results
using Plots
plot(model.params.spec_bands[1], R[1,1,:], 
     xlabel="Wavenumber [cm⁻¹]", 
     ylabel="TOA Reflectance",
     title="RT with GEOSChem-TOMAS Aerosols")
```

---

## 8. Key Design Decisions

### ✅ Use Bins Directly (No Fitting Initially)
**Rationale**: 
- TOMAS bins are already at appropriate resolution for Mie
- Fitting parametric distributions loses information
- Simpler implementation, easier validation
- **Future**: Add fitting as optional feature

### ✅ Separate Each Species
**Rationale**:
- Different optical properties (n, k, density)
- Allows selective inclusion/exclusion
- Easier to track contributions for diagnostics
- Follows GEOSChem-TOMAS philosophy

### ✅ External Mixture Assumption (Initially)
**Rationale**:
- Simplest optically: sum individual contributions
- Computationally efficient
- **Future**: Add core-shell or internal mixing models

### ✅ Wavelength-Independent n,k (Initially)
**Rationale**:
- Simplifies first implementation
- Good approximation for visible/NIR
- **Future**: Add wavelength-dependent lookup tables

---

## 9. Open Questions for User

1. **VMR to Number Density**: My conversion is approximate. Do you have GEOSChem documentation on proper conversion accounting for molecular weight and particle mass?

2. **Refractive Indices**: Are the n,k values I listed acceptable, or do you have preferred sources (OPAC, AERONET, etc.)?

3. **Hygroscopic Growth**: The `AW01-15` bins contain aerosol water. Should we:
   - Ignore them (use dry radii)?
   - Use them to compute wet radii (more complex)?
   - Mix water as a separate species?

4. **Species Priority**: Which aerosol types are most important for your RT calculations? (Dust, sea salt, sulfate, BC, OC?)

5. **Validation Data**: Do you have benchmark RT calculations with known aerosol profiles to validate against?

6. **Performance**: 15 bins × multiple species × many layers could be slow. Should we add optimization options (e.g., combine bins, skip low-concentration layers)?

---

## Next Steps

Ready to implement when you provide feedback on the open questions! 🚀
