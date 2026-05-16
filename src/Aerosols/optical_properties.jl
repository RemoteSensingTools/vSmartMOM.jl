"""
    Aerosol optical property calculations
    
    Compute extinction, scattering, absorption, and phase functions
    for different aerosol schemes using Mie theory.
"""

"""
    compute_optical_properties(data::AerosolData{TOMASScheme{FT}}, 
                               wavelengths::AbstractVector,
                               ri_database::RefractiveIndexDatabase{FT}) where FT

Compute optical properties for TOMAS size-resolved aerosols.

Uses Mie theory for each size bin, then integrates over all bins.

# Arguments
- `data::AerosolData{TOMASScheme{FT}}`: TOMAS aerosol data
- `wavelengths::AbstractVector`: Wavelengths (μm) for calculations
- `ri_database::RefractiveIndexDatabase{FT}`: Refractive index database

# Returns
- `Dict`: Optical properties with keys:
  - "extinction": Extinction coefficient (km⁻¹) [n_levels, n_wavelengths]
  - "scattering": Scattering coefficient (km⁻¹) [n_levels, n_wavelengths]
  - "absorption": Absorption coefficient (km⁻¹) [n_levels, n_wavelengths]
  - "ssa": Single scattering albedo [n_levels, n_wavelengths]
  - "asymmetry_parameter": Asymmetry parameter g [n_levels, n_wavelengths]

# Method
For each species and size bin:
1. Get refractive index at wavelength
2. Compute Mie Q_ext, Q_sca, Q_abs, g using particle size and n
3. Convert concentrations to cross-sections
4. Sum over all bins and species
"""
function compute_optical_properties(
    data::AerosolData{TOMASScheme{FT}},
    wavelengths::AbstractVector,
    ri_database::RefractiveIndexDatabase{FT}
) where FT
    scheme = data.scheme
    n_bins = scheme.n_bins
    n_levels = size(first(values(data.species_data)).data["concentration"], 2)
    n_wavelengths = length(wavelengths)
    
    # Initialize arrays
    extinction = zeros(FT, n_levels, n_wavelengths)
    scattering = zeros(FT, n_levels, n_wavelengths)
    absorption = zeros(FT, n_levels, n_wavelengths)
    
    # Loop over wavelengths
    for (iλ, λ) in enumerate(wavelengths)
        # Loop over species
        for species_name in scheme.species
            species_data_obj = data.species_data[species_name]
            concentrations = species_data_obj.data["concentration"]  # [n_bins, n_levels]
            
            # Get refractive index for this species
            ri_key = scheme.refractive_indices[species_name]
            n_complex = get_refractive_index(ri_database, ri_key, λ)
            
            # Get physical properties
            density = scheme.densities[species_name]  # kg/m³
            
            # Loop over size bins
            for ibin in 1:n_bins
                # Bin properties
                diam_nm = scheme.bin_centers[ibin]
                radius_um = diam_nm / 2000.0  # Convert nm to μm
                
                # Size parameter
                x = 2π * radius_um / λ
                
                # Compute Mie efficiencies (placeholder - needs actual Mie implementation)
                # For now, use simple approximations
                Q_ext, Q_sca, Q_abs, g = compute_mie_efficiencies(x, n_complex)
                
                # Geometric cross-section (μm²)
                σ_geom = π * radius_um^2
                
                # Optical cross-sections (μm²)
                σ_ext = Q_ext * σ_geom
                σ_sca = Q_sca * σ_geom
                σ_abs = Q_abs * σ_geom
                
                # Convert concentration to number density
                # concentration is mol/mol, need #/m³
                # This requires meteorological data (P, T) - placeholder for now
                # Assume ~10^12 molecules/cm³ air → concentrations[ibin, :] × 10^18 #/m³
                number_density = concentrations[ibin, :] .* 1e18  # rough estimate
                
                # Contribution to optical depth per unit height (m⁻¹)
                extinction[:, iλ] .+= number_density .* σ_ext .* 1e-12  # Convert μm² to m²
                scattering[:, iλ] .+= number_density .* σ_sca .* 1e-12
                absorption[:, iλ] .+= number_density .* σ_abs .* 1e-12
            end
        end
        
        # Convert from m⁻¹ to km⁻¹
        extinction[:, iλ] .*= 1000.0
        scattering[:, iλ] .*= 1000.0
        absorption[:, iλ] .*= 1000.0
    end
    
    # Compute derived quantities
    ssa = scattering ./ (extinction .+ 1e-30)  # Avoid division by zero
    ssa = clamp.(ssa, 0.0, 1.0)
    
    # Asymmetry parameter (placeholder - needs phase function integration)
    g = fill(0.7, n_levels, n_wavelengths)  # Typical value for aerosols
    
    return Dict(
        "extinction" => extinction,
        "scattering" => scattering,
        "absorption" => absorption,
        "ssa" => ssa,
        "asymmetry_parameter" => g
    )
end

"""
    compute_optical_properties(data::AerosolData{TwoMomentScheme{FT}},
                               wavelengths::AbstractVector,
                               ri_database::RefractiveIndexDatabase{FT}) where FT

Compute optical properties for two-moment bulk aerosols.

Scales AOD from reference wavelength and computes other properties
assuming lognormal distributions.

# Arguments
- `data::AerosolData{TwoMomentScheme{FT}}`: Two-moment aerosol data
- `wavelengths::AbstractVector`: Wavelengths (μm) for calculations
- `ri_database::RefractiveIndexDatabase{FT}`: Refractive index database

# Returns
- `Dict`: Same structure as the TOMAS version
"""
function compute_optical_properties(
    data::AerosolData{TwoMomentScheme{FT}},
    wavelengths::AbstractVector,
    ri_database::RefractiveIndexDatabase{FT}
) where FT
    scheme = data.scheme
    n_levels = length(first(values(data.species_data)).data["aod"])
    n_wavelengths = length(wavelengths)
    
    # Initialize arrays
    extinction = zeros(FT, n_levels, n_wavelengths)
    scattering = zeros(FT, n_levels, n_wavelengths)
    absorption = zeros(FT, n_levels, n_wavelengths)
    
    # Loop over wavelengths
    for (iλ, λ) in enumerate(wavelengths)
        # Loop over species
        for species_name in scheme.species
            species_data_obj = data.species_data[species_name]
            aod_ref = species_data_obj.data["aod"]
            r_eff = species_data_obj.data["radius"]
            
            # Reference wavelength for this species
            λ_ref = scheme.aod_wavelength[species_name]
            
            # Scale AOD to target wavelength (simple Ångström)
            aod_scaled = [scale_aod_wavelength(aod, λ_ref, λ) for aod in aod_ref]
            
            # Get refractive index
            ri_key = scheme.refractive_indices[species_name]
            n_complex = get_refractive_index(ri_database, ri_key, λ)
            
            # Compute single scattering albedo from refractive index
            # Simple approximation: ω = 1 - imag(n) / real(n)
            n_real = real(n_complex)
            n_imag = imag(n_complex)
            ssa_approx = max(0.0, 1.0 - n_imag / (n_real + 1e-10))
            
            # Partition extinction into scattering and absorption
            extinction[:, iλ] .+= aod_scaled
            scattering[:, iλ] .+= aod_scaled .* ssa_approx
            absorption[:, iλ] .+= aod_scaled .* (1.0 - ssa_approx)
        end
    end
    
    # Compute SSA
    ssa = scattering ./ (extinction .+ 1e-30)
    ssa = clamp.(ssa, 0.0, 1.0)
    
    # Asymmetry parameter (placeholder)
    g = fill(0.7, n_levels, n_wavelengths)
    
    return Dict(
        "extinction" => extinction,
        "scattering" => scattering,
        "absorption" => absorption,
        "ssa" => ssa,
        "asymmetry_parameter" => g
    )
end

"""
    compute_mie_efficiencies(x::Float64, n::ComplexF64)

Compute Mie scattering efficiencies for a homogeneous sphere.

# Arguments
- `x::Float64`: Size parameter (2πr/λ)
- `n::ComplexF64`: Complex refractive index

# Returns
- `Tuple{Float64, Float64, Float64, Float64}`: (Q_ext, Q_sca, Q_abs, g)
  - Q_ext: Extinction efficiency
  - Q_sca: Scattering efficiency
  - Q_abs: Absorption efficiency
  - g: Asymmetry parameter

# Note
This is a PLACEHOLDER. Full Mie implementation should use vSmartMOM's
existing Mie code or a dedicated Mie library.

For now, uses simple approximations:
- Rayleigh regime (x << 1): Q ~ x^4
- Geometric optics (x >> 1): Q → 2
"""
function compute_mie_efficiencies(x::Float64, n::ComplexF64)
    # PLACEHOLDER - Replace with actual Mie calculation
    
    n_real = real(n)
    n_imag = imag(n)
    
    if x < 0.1
        # Rayleigh regime
        m_sq = abs2(n)
        Q_sca = (8.0/3.0) * x^4 * abs2((m_sq - 1) / (m_sq + 2))
        Q_abs = 4.0 * x * imag((m_sq - 1) / (m_sq + 2))
        Q_ext = Q_sca + Q_abs
        g = 0.0  # Rayleigh scattering is symmetric
    elseif x > 20.0
        # Geometric optics approximation
        Q_ext = 2.0
        Q_abs = min(2.0, 4.0 * n_imag)
        Q_sca = Q_ext - Q_abs
        g = 0.85  # Forward scattering dominates
    else
        # Intermediate regime - rough interpolation
        f = (x - 0.1) / (20.0 - 0.1)
        Q_ext = 2.0 * f + (1.0 - f) * 0.1
        Q_abs = min(Q_ext, 2.0 * n_imag)
        Q_sca = Q_ext - Q_abs
        g = 0.85 * f  # Asymmetry increases with size
    end
    
    return (Q_ext, Q_sca, Q_abs, g)
end

"""
    integrate_phase_function(data::AerosolData, wavelength::Float64, 
                            ri_database::RefractiveIndexDatabase, 
                            scattering_angles::AbstractVector)

Compute phase function at specified scattering angles.

# Arguments
- `data::AerosolData`: Aerosol data (any scheme type)
- `wavelength::Float64`: Wavelength (μm)
- `ri_database::RefractiveIndexDatabase`: Refractive index database
- `scattering_angles::AbstractVector`: Scattering angles (degrees)

# Returns
- `Array{Float64, 2}`: Phase function [n_levels, n_angles]

# Note
This is a placeholder. Should integrate with vSmartMOM's existing
phase function and scattering matrix calculations.
"""
function integrate_phase_function(
    data::AerosolData,
    wavelength::Float64,
    ri_database::RefractiveIndexDatabase,
    scattering_angles::AbstractVector
)
    # Placeholder - return Henyey-Greenstein approximation
    g = 0.7  # asymmetry parameter
    
    n_levels = length(first(values(data.species_data)).data[first(keys(first(values(data.species_data)).data))])
    n_angles = length(scattering_angles)
    
    phase = zeros(Float64, n_levels, n_angles)
    
    for (iθ, θ) in enumerate(scattering_angles)
        cosθ = cosd(θ)
        # Henyey-Greenstein phase function
        P = (1 - g^2) / (1 + g^2 - 2g * cosθ)^1.5
        phase[:, iθ] .= P
    end
    
    return phase
end
