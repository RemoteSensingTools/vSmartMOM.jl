"""
    Two-moment aerosol scheme reader

# Physical model

The two-moment scheme characterizes each aerosol species by two moments of the
size distribution: total optical depth (AOD) and effective radius. It assumes
a lognormal size distribution with fixed geometric standard deviation σ_g.

This bulk approach is computationally efficient compared to size-resolved schemes:
- AOD at reference wavelength (e.g. 550 nm) scales via Ångström exponent
- Effective radius r_eff and σ_g define the full dN/dr distribution
- Optical properties computed by integrating Mie theory over the lognormal

Common in climate models (e.g. CAM, ECHAM) and satellite retrieval algorithms.
"""

"""
    read_two_moment(config::Dict, netcdf_file::String, FT=Float64)

Read two-moment aerosol data from NetCDF file.

# Arguments
- `config::Dict`: YAML configuration dictionary
- `netcdf_file::String`: Path to NetCDF file
- `FT`: Floating point type (default: Float64)

# Returns
- `AerosolData{TwoMomentScheme{FT}}`: Container with AOD and effective radius data

# Data Structure
Each species has:
- AOD at reference wavelength (typically 550 nm)
- Effective radius (μm)
- Fixed geometric standard deviation σ_g

NetCDF variables:
- AODHyg550nm_[SPECIES] or similar
- Chem_AeroRadi[SPECIES] or similar

# Configuration
The config dict must specify:
- aerosol_scheme/species: Dict of species with σ_g and wavelength
- netcdf_mapping: Variable name patterns for AOD and radius
"""
function read_two_moment(config::Dict, netcdf_file::String, FT=Float64)
    # Construct scheme from config
    scheme = TwoMomentScheme(config, FT)
    
    # Open NetCDF file
    ds = NCDataset(netcdf_file)
    
    try
        # Extract coordinates
        coordinates = extract_coordinates(ds, config)
        n_levels = length(coordinates["lev"])
        
        # Get processing options
        vertical_flip = get(aerosol_processing_options(config), "vertical_flip", false)
        
        # Read each species
        species_data = Dict{String, AerosolSpeciesData}()
        
        for species_name in scheme.species
            species_config = config["aerosol_scheme"]["species"][species_name]
            
            # Get variable names
            aod_var = replace(
                species_config["aod_variable"],
                "{species}" => species_name
            )
            
            radius_var = replace(
                species_config["radius_variable"],
                "{species}" => species_name
            )
            
            # Read AOD
            aod_profile = nothing
            if haskey(ds, aod_var)
                data_full = Array(ds[aod_var])
                # Average over horizontal dimensions
                aod_profile = dropdims(
                    mean(data_full, dims=(1, 2, 3)),
                    dims=(1, 2, 3)
                )[:, 1]  # First time step
            else
                @warn "AOD variable $aod_var not found for species $species_name"
                aod_profile = zeros(FT, n_levels)
            end
            
            # Read effective radius
            radius_profile = nothing
            if haskey(ds, radius_var)
                data_full = Array(ds[radius_var])
                # Average over horizontal dimensions
                radius_profile = dropdims(
                    mean(data_full, dims=(1, 2, 3)),
                    dims=(1, 2, 3)
                )[:, 1]  # First time step
            else
                @warn "Radius variable $radius_var not found for species $species_name"
                radius_profile = zeros(FT, n_levels)
            end
            
            # Apply vertical flip if requested
            if vertical_flip
                aod_profile = reverse(aod_profile)
                radius_profile = reverse(radius_profile)
            end
            
            # Store species data
            data_dict = Dict{String, Any}(
                "aod" => aod_profile,
                "radius" => radius_profile
            )
            
            units_dict = Dict{String, String}(
                "aod" => "dimensionless (at $(species_config["aod_reference_wavelength"]) μm)",
                "radius" => "μm"
            )
            
            species_data[species_name] = AerosolSpeciesData(
                data_dict,
                units_dict,
                "Two-moment $(species_name) aerosol (σ_g = $(species_config["sigma_g"]))"
            )
        end
        
        # Extract metadata
        metadata = extract_metadata(ds)
        
        return AerosolData(scheme, species_data, coordinates, metadata)
        
    finally
        close(ds)
    end
end

"""
    scale_aod_wavelength(aod_ref::Float64, λ_ref::Float64, λ_target::Float64, 
                         angstrom_exponent::Float64 = 1.0)

Scale AOD from reference wavelength to target wavelength using Ångström exponent.

# Arguments
- `aod_ref::Float64`: AOD at reference wavelength
- `λ_ref::Float64`: Reference wavelength (μm)
- `λ_target::Float64`: Target wavelength (μm)
- `angstrom_exponent::Float64`: Ångström exponent (default 1.0)

# Returns
- `Float64`: Scaled AOD at target wavelength

# Formula
AOD(λ) = AOD(λ_ref) × (λ / λ_ref)^(-α)

where α is the Ångström exponent
"""
function scale_aod_wavelength(
    aod_ref::Float64,
    λ_ref::Float64,
    λ_target::Float64,
    angstrom_exponent::Float64 = 1.0
)
    return aod_ref * (λ_target / λ_ref)^(-angstrom_exponent)
end

"""
    lognormal_size_distribution(r::AbstractVector, r_eff::Float64, σ_g::Float64)

Compute lognormal number size distribution dN/dr.

# Arguments
- `r::AbstractVector`: Particle radii (μm)
- `r_eff::Float64`: Effective radius (μm)
- `σ_g::Float64`: Geometric standard deviation

# Returns
- `Vector{Float64}`: Number size distribution dN/dr

# Formula
dN/dr = (N_total / (r × sqrt(2π) × ln(σ_g))) × exp(-ln²(r/r_med) / (2 × ln²(σ_g)))

where r_med is related to r_eff by:
r_eff = r_med × exp(2.5 × ln²(σ_g))
"""
function lognormal_size_distribution(
    r::AbstractVector,
    r_eff::Float64,
    σ_g::Float64
)
    # Convert effective radius to median radius
    ln_sigma_g = log(σ_g)
    r_med = r_eff / exp(2.5 * ln_sigma_g^2)
    
    # Compute distribution (normalized)
    dN_dr = @. (1.0 / (r * sqrt(2π) * ln_sigma_g)) * 
               exp(-log(r / r_med)^2 / (2 * ln_sigma_g^2))
    
    return dN_dr
end

"""
    effective_radius_from_moments(r_med::Float64, σ_g::Float64)

Calculate effective radius from median radius and geometric std dev.

# Arguments
- `r_med::Float64`: Median radius (μm)
- `σ_g::Float64`: Geometric standard deviation

# Returns
- `Float64`: Effective radius (μm)

# Formula
r_eff = r_med × exp(2.5 × ln²(σ_g))
"""
function effective_radius_from_moments(r_med::Float64, σ_g::Float64)
    ln_sigma_g = log(σ_g)
    return r_med * exp(2.5 * ln_sigma_g^2)
end

"""
    median_radius_from_effective(r_eff::Float64, σ_g::Float64)

Calculate median radius from effective radius and geometric std dev.

# Arguments
- `r_eff::Float64`: Effective radius (μm)
- `σ_g::Float64`: Geometric standard deviation

# Returns
- `Float64`: Median radius (μm)

# Formula
r_med = r_eff / exp(2.5 × ln²(σ_g))
"""
function median_radius_from_effective(r_eff::Float64, σ_g::Float64)
    ln_sigma_g = log(σ_g)
    return r_eff / exp(2.5 * ln_sigma_g^2)
end
