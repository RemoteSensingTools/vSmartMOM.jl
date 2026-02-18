"""
    TOMAS-15 aerosol scheme reader

# Physical model

TOMAS (TwO-Moment Aerosol Sectional) is a size-resolved microphysics scheme that
discretizes the aerosol size distribution into 15 logarithmically-spaced bins
(typically 10 nm to 10 μm dry diameter). Each bin tracks particle number and
mass of each chemical species (sulfate, organic carbon, dust, sea salt, etc.).

This enables accurate representation of:
- Nucleation, condensation, coagulation, and deposition
- Size-dependent optical properties via Mie theory
- CCN activation and cloud-aerosol interactions

The scheme is commonly used in GEOS-Chem and similar CTMs.
"""

"""
    read_tomas15(config::Dict, netcdf_file::String, FT=Float64)

Read TOMAS-15 aerosol data from NetCDF file including NK (number) and species.

# Arguments
- `config::Dict`: YAML configuration dictionary
- `netcdf_file::String`: Path to NetCDF file
- `FT`: Floating point type (default: Float64)

# Returns
- `AerosolData{TOMAS15Scheme{FT}}`: Container with:
  - NK: Total particle number per bin (#/cm³)
  - dN_dlogD: Size distribution (cm⁻³)
  - Species mass: Mass per bin per species (μg/m³)
  - Species fractions: Fractional composition per bin

# Data Structure
NK variables: SpeciesConcVV_NK01 through NK15 (1000 × particles/mol_air)
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
- netcdf_mapping: Variable name patterns
- processing_options/vertical_flip: Whether to flip BOA→TOA to TOA→BOA
"""
function read_tomas15(config::Dict, netcdf_file::String, FT=Float64)
    # Construct scheme from config
    scheme = TOMAS15Scheme(config, FT)
    
    # Open NetCDF file
    ds = NCDataset(netcdf_file)
    
    try
        # Extract coordinates
        coordinates = extract_coordinates(ds, config)
        n_levels = length(coordinates["lev"])
        
        # Get processing options
        vertical_flip = get(config["processing_options"], "vertical_flip", false)
        
        # Get meteorology variables for conversion
        met_vars = config["netcdf_mapping"]["meteorology"]
        
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
            # NetCDF variable pattern
            var_pattern = replace(
                config["netcdf_mapping"]["concentration_pattern"],
                "{species}" => species_name,
                "{bin:02d}" => ""
            )
            
            # Read all bins for this species
            concentrations = zeros(FT, scheme.n_bins, n_levels)  # mol/mol
            mass_conc = zeros(FT, scheme.n_bins, n_levels)  # μg/m³
            part_num = zeros(FT, scheme.n_bins, n_levels)  # #/cm³
            
            for ibin in 1:scheme.n_bins
                var_name = replace(var_pattern, "{bin:02d}" => @sprintf("%02d", ibin))
                
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
            data_dict = Dict{String, Array}(
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
                "TOMAS-15 size-resolved $(species_name) aerosol"
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
        
        nk_data_dict = Dict{String, Array}(
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
            "TOMAS-15 total particle number distribution with species fractions"
        )
        
        # Extract metadata
        metadata = extract_metadata(ds)
        
        return AerosolData(scheme, species_data, coordinates, metadata)
        
    finally
        close(ds)
    end
end

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
