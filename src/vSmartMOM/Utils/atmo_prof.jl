"Compute pressure levels, vmr, vcd for atmospheric profile, given p_half, T, q"
function compute_atmos_profile_fields(p_half::AbstractArray, T, q; g₀=9.8196)

    # Floating type to use
    FT = eltype(T)
    
    # Calculate full pressure levels
    p_full = (p_half[2:end] + p_half[1:end - 1]) / 2

    # Dry and wet mass
    dry_mass = 28.9647e-3  / Nₐ  # in kg/molec, weighted average for N2 and O2
    wet_mass = 18.01528e-3 / Nₐ  # just H2O
    ratio = dry_mass / wet_mass 
    n_layers = length(T)

    # Also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end

    return p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o

end

"Validate an input atmospheric profile (YAML only)"
function validate_atmos_profile(params::Dict{Any,Any})

    @assert Set(keys(params)) == Set(["p_surf", "T", "q", "ak", "bk", "vmr"]) || 
            Set(keys(params)) == Set(["p_half", "T", "q", "vmr"]) || 
            Set(keys(params)) == Set(["p_half", "T", "vmr"]) "Set of atmospheric profile fields must be one of (p_surf, T, q, ak, bk), (p_half, T, q), or (p_half, T)"

end

"From a yaml file, get the stored fields (psurf, T, q, ak, bk), calculate derived fields, 
and return an AtmosphericProfile object" 
function read_atmos_profile(file_path::String)

    # Make sure file is yaml type
    @assert endswith(file_path, ".yaml") "File must be yaml"

    # Read in the data and pass to compute fields
    params_dict = YAML.load_file(file_path)

    # Validate the parameters before doing anything else
    validate_atmos_profile(params_dict)

    T = convert.(Float64, params_dict["T"])
    
    # Calculate derived fields
    if ("ak" in keys(params_dict))
        psurf = convert(Float64, params_dict["p_surf"])
        q     = convert.(Float64, params_dict["q"])
        ak    = convert.(Float64, params_dict["ak"])
        bk    = convert.(Float64, params_dict["bk"])
        p_half = (ak + bk * psurf)
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o = compute_atmos_profile_fields(p_half, T, q)
    elseif ("q" in keys(params_dict))
        p_half = convert(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q      = convert.(Float64, params_dict["q"])
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o = compute_atmos_profile_fields(p_half, T, q)
    else
        p_half = convert.(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q = zeros(length(T))
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o = compute_atmos_profile_fields(p_half, T, q)
    end

    # Convert vmr to appropriate type
    vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["vmr"])

    # Return the atmospheric profile struct
    return AtmosphericProfile(nothing, nothing, psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr)

end

# "Read atmospheric profile from netCDF file (ONLY works for FraLab test file)"
# function read_atmos_profile(file::String, lat::Real, lon::Real, time_idx::Int)

#     # Time index must be ∈ [1, 2, 3, 4]
#     @assert 1 <= time_idx <= 4 "Time index must be ∈ [1, 2, 3, 4]" 

#     # Make sure file is nc type
#     @assert endswith(file, ".nc4") "File must be nc4"

#     # Load in the atmospheric profile
#     @timeit "loading file" ds = Dataset(file)

#     # See how easy it is to actually extract data? 
#     # Note the [:] in the end reads in ALL the data in one step
#     file_lats, file_lons = ds["YDim"][:], ds["XDim"][:]
    
#     # Convert the input lat/lon to right type
#     FT = eltype(file_lats)
#     lat, lon = FT(lat), FT(lon)
    
#     # Find index (nearest neighbor, one could envision interpolation in space and time!):
#     lat_idx, lon_idx = argmin(abs.(file_lats .- lat)), argmin(abs.(file_lons .- lon))

#     # Temperature profile
#     @timeit "getting T" T = convert(Array{FT,1}, ds["T"][lon_idx, lat_idx,  :, time_idx])

#     # Specific humidity profile
#     q = convert(Array{FT,1}, ds["QV"][lon_idx, lat_idx, :, time_idx])
    
#     # Surface pressure
#     psurf = convert(FT, ds["PS"][lon_idx, lat_idx, time_idx])
    
#     # AK and BK global attributes (important to calculate pressure half-levels)
#     ak, bk = ds.attrib["HDF_GLOBAL.ak"][:], ds.attrib["HDF_GLOBAL.bk"][:]

#     # Close the file
#     close(ds)

#     # Calculate derived fields
#     p_half = (ak + bk * psurf)
#     p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o = compute_atmos_profile_fields(p_half, T, q)

#     # Return the atmospheric profile struct
#     return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o)
# end

"Reduce profile dimensions by re-averaging to near-equidistant pressure grid"
function reduce_profile(n::Int, profile::AtmosphericProfile{FT}) where {FT}

    # Can only reduce the profile, not expand it
    @assert n < length(profile.T)

    # Unpack the profile
    @unpack lat, lon, psurf, vmr = profile

    # New rough half levels (boundary points)
    a = range(0, maximum(profile.p), length=n + 1)

    # Matrices to hold new values
    T = zeros(FT, n);
    q = zeros(FT, n);
    p_full = zeros(FT, n);
    p_levels = zeros(FT, n + 1);
    vmr_h2o  = zeros(FT, n);
    vcd_dry  = zeros(FT, n);
    vcd_h2o  = zeros(FT, n);

    # Loop over target number of layers
    for i = 1:n

        # Get the section of the atmosphere with the i'th section pressure values
        ind = findall(a[i] .< profile.p .<= a[i + 1]);

        # Set the pressure levels accordingly
        p_levels[i] = profile.p_levels[ind[1]]
        p_levels[i + 1] = profile.p_levels[ind[end]]

        # Re-average the other parameters to produce new layers
        p_full[i] = mean(profile.p_levels[ind])
        T[i] = mean(profile.T[ind])
        q[i] = mean(profile.q[ind])
        vmr_h2o[i] = mean(profile.vmr_h2o[ind])
        vcd_dry[i] = sum(profile.vcd_dry[ind])
        vcd_h2o[i] = sum(profile.vcd_h2o[ind])
    end

    new_vmr = Dict{String, Union{Real, Vector}}()

    for molec_i in keys(vmr)
        if profile.vmr[molec_i] isa AbstractArray
            
            pressure_grid = collect(range(minimum(p_full), maximum(p_full), length=length(profile.vmr[molec_i])))
            interp_linear = LinearInterpolation(pressure_grid, vmr[molec_i])
            new_vmr[molec_i] = [interp_linear(x) for x in p_full]
        else
            new_vmr[molec_i] = profile.vmr[molec_i]
        end
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_levels, vmr_h2o, vcd_dry, vcd_h2o, new_vmr)
end;

"""
$(FUNCTIONNAME)(psurf, λ, depol_fct, vcd_dry)
Returns the Rayleigh optical thickness per layer at reference wavelength `λ` (N₂,O₂ atmosphere, i.e. terrestrial)
Input: 
- `psurf` surface pressure in `[hPa]`
- `λ` wavelength in `[μm]`
- `depol_fct` depolarization factor
- `vcd_dry` dry vertical column (no water) per layer
"""
function getRayleighLayerOptProp(psurf, λ, depol_fct, vcd_dry) 
    FT = eltype(λ)
    # Total vertical Rayleigh scattering optical thickness 
    tau_scat = FT(0.00864) * (psurf / FT(1013.25)) * λ^(-FT(3.916) - FT(0.074) * λ - FT(0.05) / λ) 
    tau_scat = tau_scat * (FT(6.0) + FT(3.0) * depol_fct) / (FT(6.0)- FT(7.0) * depol_fct)
    Nz = length(vcd_dry)
    τRayl = zeros(FT,Nz)
    k = tau_scat / sum(vcd_dry)
    for i = 1:Nz
        τRayl[i] = k * vcd_dry[i]
    end

    return convert.(FT, τRayl)
end

"""
$(FUNCTIONNAME)(total_τ, p₀, σp, p_half)
Returns the aerosol optical depths per layer using a Gaussian distribution function with p₀ and σp on a pressure grid
"""
function getAerosolLayerOptProp(total_τ, p₀, σp, p_half)

    # Need to make sure we can also differentiate wrt σp (FT can be Dual!)
    FT = eltype(p₀)
    Nz = length(p_half)
    ρ = zeros(FT,Nz)

    for i = 2:Nz
        dp = p_half[i] - p_half[i - 1]
        ρ[i] = (1 / (σp * sqrt(2π))) * exp(-(p_half[i] - p₀)^2 / (2σp^2)) * dp
    end
    Norm = sum(ρ)
    τAer  =  (total_τ / Norm) * ρ
    return convert.(FT, τAer)
end

"""
$(FUNCTIONNAME)(τRayl, τAer,  aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs, arr_type)
Computes the composite layer single scattering parameters (τ, ϖ, Z⁺⁺, Z⁻⁺)

Returns:
- `τ`, `ϖ`   : only Rayleigh scattering and aerosol extinction, no gaseous absorption (no wavelength dependence)
- `τ_λ`,`ϖ_λ`: Rayleigh scattering + aerosol extinction + gaseous absorption (wavelength dependent)
- `Z⁺⁺`,`Z⁻⁺`: Composite Phase matrix (weighted average of Rayleigh and aerosols)

# Arguments
- `τRay` layer optical depth for Rayleigh
- `τAer` layer optical depth for Aerosol(s) (vector)
- `aerosol_optics` array of aerosol optics struct
- `Rayl𝐙⁺⁺` Rayleigh 𝐙⁺⁺ phase matrix (2D)
- `Rayl𝐙⁻⁺` Rayleigh 𝐙⁻⁺ phase matrix (2D)
- `Aer𝐙⁺⁺` Aerosol 𝐙⁺⁺ phase matrix (3D)
- `Aer𝐙⁻⁺` Aerosol 𝐙⁻⁺ phase matrix (3D)
- `τ_abs` layer absorption optical depth array (per wavelength) by gaseous absorption
"""
function construct_atm_layer(τRayl, τAer,  aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs, arr_type)
    FT = eltype(τRayl)
    nAer = length(aerosol_optics)

    # Fix Rayleigh SSA to 1
    ϖRayl = FT(1)

    @assert length(τAer) == nAer "Sizes don't match"

    τ = FT(0)
    ϖ = FT(0)
    A = FT(0)
    Z⁺⁺ = similar(Rayl𝐙⁺⁺); 
    Z⁻⁺ = similar(Rayl𝐙⁺⁺);

    if (τRayl + sum(τAer)) < eps(FT)
        fill!(Z⁺⁺, 0); fill!(Z⁻⁺, 0);
        return FT(0), FT(1), Z⁺⁺, Z⁻⁺
    end
 
    τ += τRayl
    ϖ += τRayl * ϖRayl
    A += τRayl * ϖRayl

    Z⁺⁺ = τRayl * ϖRayl * Rayl𝐙⁺⁺
    Z⁻⁺ = τRayl * ϖRayl * Rayl𝐙⁻⁺

    for i = 1:nAer
        τ   += τAer[i]
        ϖ   += τAer[i] * aerosol_optics[i].ω̃
        A   += τAer[i] * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ)
        Z⁺⁺ += τAer[i] * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ) * Aer𝐙⁺⁺[:,:,i]
        Z⁻⁺ += τAer[i] * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ) * Aer𝐙⁻⁺[:,:,i]
    end
    
    Z⁺⁺ /= A
    Z⁻⁺ /= A
    A /= ϖ
    ϖ /= τ
    
    # Rescaling composite SSPs according to Eqs. A.3 of Sanghavi et al. (2013) or Eqs.(8) of Sanghavi & Stephens (2015)
    τ *= (FT(1) - (FT(1) - A) * ϖ)
    ϖ *= A / (FT(1) - (FT(1) - A) * ϖ)#Suniti

    # Adding absorption optical depth / albedo:
    τ_λ = τ_abs .+ τ    
    ϖ_λ = (τ .* ϖ) ./ τ_λ
    
    return Array(τ_λ), Array(ϖ_λ), τ, ϖ, Array(Z⁺⁺), Array(Z⁻⁺)
end

"When performing RT_run, this function pre-calculates properties for all layers, before any Core RT is performed"
function construct_all_atm_layers(FT, nSpec, Nz, NquadN, τRayl, τAer, aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs, arr_type, qp_μ, μ₀, m)

    FT_ext   = eltype(τAer)
    FT_phase = eltype(Aer𝐙⁺⁺)

    @show FT_ext
    @show FT_phase

    # Empty matrices to hold all values
    τ_λ_all   = zeros(FT_ext, nSpec, Nz)
    ϖ_λ_all   = zeros(FT_ext, nSpec, Nz)
    τ_all     = zeros(FT_ext, Nz)
    ϖ_all     = zeros(FT_ext, Nz)
    Z⁺⁺_all   = zeros(FT_phase, NquadN, NquadN, Nz)
    Z⁻⁺_all   = zeros(FT_phase, NquadN, NquadN, Nz)
    
    dτ_max_all  = zeros(FT_ext, Nz)
    dτ_all      = zeros(FT_ext, Nz)
    ndoubl_all  = zeros(Int64, Nz)
    dτ_λ_all    = zeros(FT_ext, nSpec, Nz)
    expk_all    = zeros(FT_ext, nSpec, Nz)
    scatter_all = zeros(Bool, Nz)

    for iz=1:Nz
        
        # Construct atmospheric properties
        τ_λ_all[:, iz], ϖ_λ_all[:, iz], τ_all[iz], ϖ_all[iz], Z⁺⁺_all[:,:,iz], Z⁻⁺_all[:,:,iz] = construct_atm_layer(τRayl[iz], τAer[:,iz], aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs[:,iz], arr_type)

        # Compute doubling number
        dτ_max_all[iz] = minimum([τ_all[iz] * ϖ_all[iz], FT(0.001) * minimum(qp_μ)])
        dτ_all[iz], ndoubl_all[iz] = doubling_number(dτ_max_all[iz], τ_all[iz] * ϖ_all[iz]) #Suniti

        # Compute dτ vector
        dτ_λ_all[:, iz] = (τ_λ_all[:, iz] ./ (FT(2)^ndoubl_all[iz]))
        #@show maximum(dτ_λ_all[:,iz])
        expk_all[:, iz] = exp.(-dτ_λ_all[:, iz] /μ₀) #Suniti
        
        # Determine whether there is scattering
        scatter_all[iz] = (  sum(τAer[:,iz]) > 1.e-8 || 
                          (( τRayl[iz] > 1.e-8 ) && (m < 3))) ? 
                            true : false
    end

    # Compute sum of optical thicknesses of all layers above the current layer
    τ_sum_all = accumulate(+, τ_λ_all, dims=2)

    # First start with all zeros
    # At the bottom of the atmosphere, we have to compute total τ_sum (bottom of lowest layer), for the surface interaction
    τ_sum_all = hcat(zeros(FT, size(τ_sum_all[:,1])), τ_sum_all)

    # Starting scattering interface (None for both added and composite)
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = []

    for iz = 1:Nz
        # Whether there is scattering in the added layer, composite layer, neither or both
        scattering_interface = get_scattering_interface(scattering_interface, scatter_all[iz], iz)
        push!(scattering_interfaces_all, scattering_interface)
    end

    return ComputedAtmosphereProperties(τ_λ_all, ϖ_λ_all, τ_all, ϖ_all, Z⁺⁺_all, Z⁻⁺_all, dτ_max_all, dτ_all, ndoubl_all, dτ_λ_all, expk_all, scatter_all, τ_sum_all, scattering_interfaces_all)
end

"Given the CrossSectionModel, the grid, and the AtmosphericProfile, fill up the τ_abs array with the cross section at each layer
(using pressures/temperatures) from the profile" 
function compute_absorption_profile!(τ_abs::Array{FT,2}, 
                                     hitran_data::HitranTable, 
                                     broadening_function::AbstractBroadeningFunction, 
                                     wing_cutoff, 
                                     CEF::AbstractComplexErrorFunction, 
                                     architecture,
                                     vmr,
                                     grid,
                                     profile::AtmosphericProfile,
                                     ) where FT <: AbstractFloat

    # The array to store the cross-sections must be same length as number of layers
    @assert size(τ_abs,2) == length(profile.p)

    @showprogress 1 for iz in 1:length(profile.p)

        # Pa -> hPa
        p = profile.p[iz] / 100
        T = profile.T[iz]

        # Either use the current layer's vmr, or use the uniform vmr
        vmr_curr = vmr isa AbstractArray ? vmr[iz] : vmr

        # Create absorption model with parameters
        absorption_model = make_hitran_model(hitran_data, 
                                             broadening_function, 
                                             wing_cutoff = wing_cutoff, 
                                             CEF = CEF, 
                                             architecture = architecture, 
                                             vmr = vmr_curr)

        # Changed index order
        τ_abs[:,iz] += Array(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
    end

    return nothing
    
end
