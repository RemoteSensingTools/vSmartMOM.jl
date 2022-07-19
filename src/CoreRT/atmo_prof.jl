#=

This file contains functions that are related to atmospheric profile calculations

=#

"Compute pressure levels, vmr, vcd for atmospheric profile, given p_half, T, q"
function compute_atmos_profile_fields(T, p_half::AbstractArray, q, vmr; g₀=9.8032465)
    
    # Floating type to use
    FT = eltype(T)
    Nₐ = FT(6.02214179e+23)
    # Calculate full pressure levels
    p_full = (p_half[2:end] + p_half[1:end-1]) / 2

    # Dry and wet mass
    dry_mass = FT(28.9644e-3)    # in kg/molec, weighted average for N2 and O2
    wet_mass = FT(18.01534e-3)   # just H2O
    n_layers = length(T)

    # Also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = dry_mass/(dry_mass-wet_mass*(1-1/q[i]))
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd = Nₐ * Δp / (M  * g₀ * 100^2) * 100
        vcd_dry[i] = vmr_dry    * vcd   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * vcd
    end

    new_vmr = Dict{String, Union{Real, Vector}}()

    for molec_i in keys(vmr)
        if vmr[molec_i] isa AbstractArray
            
            pressure_grid = collect(range(minimum(p_full), maximum(p_full), length=length(vmr[molec_i])))
            interp_linear = LinearInterpolation(pressure_grid, vmr[molec_i])
            new_vmr[molec_i] = [interp_linear(x) for x in p_full]
        else
            new_vmr[molec_i] = vmr[molec_i]
        end
    end

    return p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr

end

"From a yaml file, get the stored fields (psurf, T, q, ak, bk), calculate derived fields, 
and return an AtmosphericProfile object" 
function read_atmos_profile(file_path::String)

    # Make sure file is yaml type
    @assert endswith(file_path, ".yaml") "File must be yaml"

    # Read in the data and pass to compute fields
    params_dict = YAML.load_file(file_path)

    # Validate the parameters before doing anything else
    # validate_atmos_profile(params_dict)

    T = convert.(Float64, params_dict["T"])
    
    # Calculate derived fields
    if ("ak" in keys(params_dict))
        psurf = convert(Float64, params_dict["p_surf"])
        q     = convert.(Float64, params_dict["q"])
        ak    = convert.(Float64, params_dict["ak"])
        bk    = convert.(Float64, params_dict["bk"])
        p_half = (ak + bk * psurf)
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o = compute_atmos_profile_fields(T, p_half, q, Dict())
    elseif ("q" in keys(params_dict))
        p_half = convert(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q      = convert.(Float64, params_dict["q"])
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o = compute_atmos_profile_fields(T, p_half, q, Dict())
    else
        p_half = convert.(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q = zeros(length(T))
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o = compute_atmos_profile_fields(T, p_half, q, Dict())
    end

    # Convert vmr to appropriate type
    vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["vmr"])

    # Return the atmospheric profile struct
    return AtmosphericProfile(T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr)

end

"Reduce profile dimensions by re-averaging to near-equidistant pressure grid"
function reduce_profile(n::Int, profile::AtmosphericProfile{FT}) where {FT}

    # Can only reduce the profile, not expand it
    @assert n < length(profile.T)

    # Unpack the profile vmr
    @unpack vmr = profile

    # New rough half levels (boundary points)
    a = range(0, maximum(profile.p_half), length=n+1)

    # Matrices to hold new values
    T = zeros(FT, n);
    q = zeros(FT, n);
    p_full = zeros(FT, n);
    p_half = zeros(FT, n+1);
    vmr_h2o  = zeros(FT, n);
    vcd_dry  = zeros(FT, n);
    vcd_h2o  = zeros(FT, n);

    # Loop over target number of layers
    indices = []
    for i = 1:n

        # Get the section of the atmosphere with the i'th section pressure values
        ind = findall(a[i] .< profile.p_full .<= a[i+1]);
        push!(indices, ind)
        @assert length(ind) > 0 "Profile reduction has an empty layer"
        #@show i, ind, a[i], a[i+1]
        # Set the pressure levels accordingly
        p_half[i]   = a[i]   # profile.p_half[ind[1]]
        p_half[i+1] = a[i+1] # profile.p_half[ind[end]]

        # Re-average the other parameters to produce new layers
        p_full[i] = mean(profile.p_full[ind])
        T[i] = mean(profile.T[ind])
        q[i] = mean(profile.q[ind])
        vmr_h2o[i] = mean(profile.vmr_h2o[ind])
        vcd_dry[i] = sum(profile.vcd_dry[ind])
        vcd_h2o[i] = sum(profile.vcd_h2o[ind])
    end
    #@show indices

    new_vmr = Dict{String, Union{Real, Vector}}()

    # need to double check this logic, maybe better to add VCDs?!
    for molec_i in keys(vmr)
        if profile.vmr[molec_i] isa AbstractArray
            
            #pressure_grid = collect(range(minimum(p_full), maximum(p_full), length=length(profile.vmr[molec_i])))
            #interp_linear = LinearInterpolation(pressure_grid, vmr[molec_i])
            new_vmr[molec_i] = [mean(profile.vmr[molec_i][ind]) for ind in indices]
        else
            new_vmr[molec_i] = profile.vmr[molec_i]
        end
    end

    return AtmosphericProfile(T, p_full, q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr)
end

"""
    $(FUNCTIONNAME)(psurf, λ, depol_fct, vcd_dry)

Returns the Rayleigh optical thickness per layer at reference wavelength `λ` (N₂,O₂ atmosphere, i.e. terrestrial)

Input: 
    - `psurf` surface pressure in `[hPa]`
    - `λ` wavelength in `[μm]`
    - `depol_fct` depolarization factor
    - `vcd_dry` dry vertical column (no water) per layer
"""
function getRayleighLayerOptProp(psurf::FT, λ::Union{Array{FT}, FT}, depol_fct::FT, vcd_dry::Array{FT}) where FT
    # TODO: Use noRS/noRS_plus to use n2/o2 molecular constants
    # to compute tau_scat and depol_fct
    Nz = length(vcd_dry)
    τRayl = zeros(FT,size(λ,1),Nz)
    # Total vertical Rayleigh scattering optical thickness, TODO: enable sub-layers and use VCD based taus
    tau_scat = FT(0.00864) * (psurf / FT(1013.25)) *  λ.^(-FT(3.916) .- FT(0.074) * λ .- FT(0.05) ./ λ)  
    tau_scat = tau_scat * (FT(6.0) + FT(3.0) * depol_fct) / (FT(6.0)- FT(7.0) * depol_fct) 
    # @show tau_scat, λ
    k = tau_scat / sum(vcd_dry)
    for i = 1:Nz
        τRayl[:,i] .= k * vcd_dry[i]
    end 
    return τRayl
end

"""
    $(FUNCTIONNAME)(total_τ, p₀, σp, p_half)
    
Returns the aerosol optical depths per layer using a Gaussian distribution function with p₀ and σp on a pressure grid
"""
function getAerosolLayerOptProp(total_τ, p₀, σp, p_half)

    # Need to make sure we can also differentiate wrt σp (FT can be Dual!)
    FT = eltype(p₀)
    Nz = length(p_half)-1
    ρ = zeros(FT,Nz)

    # @show p_half, p₀, σp
    for i = 1:Nz
        dp = p_half[i+1] - p_half[i]
        p  = (p_half[i+1] + p_half[i])/2
        # Use Distributions here later:
        ρ[i] = (1 / (σp * sqrt(2π))) * exp(-(p - p₀)^2 / (2σp^2)) * dp
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
    - `fscattRayl`: Rayleigh scattering fraction (needed for Raman computations) 
Arguments:
    - `τRay` layer optical depth for Rayleigh
    - `τAer` layer optical depth for Aerosol(s) (vector)
    - `aerosol_optics` array of aerosol optics struct
    - `Rayl𝐙⁺⁺` Rayleigh 𝐙⁺⁺ phase matrix (2D)
    - `Rayl𝐙⁻⁺` Rayleigh 𝐙⁻⁺ phase matrix (2D)
    - `Aer𝐙⁺⁺` Aerosol 𝐙⁺⁺ phase matrix (3D)
    - `Aer𝐙⁻⁺` Aerosol 𝐙⁻⁺ phase matrix (3D)
    - `τ_abs` layer absorption optical depth array (per wavelength) by gaseous absorption
"""
function construct_atm_layer(τRayl, τAer,  
    ϖ_Cabannes, #elastic fraction of Rayleigh scattering
    aerosol_optics, 
    Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, 
    Aer𝐙⁺⁺, Aer𝐙⁻⁺, 
    τ_abs, arr_type)
    
    FT = eltype(τRayl)
    nAer = length(aerosol_optics)

    # Fixes Rayleigh SSA to 1 for purely elastic (RS_type = noRS) scattering,
    # and assumes values less than 1 for Raman scattering
    ϖRayl = ϖ_Cabannes #FT(1)
    #@show ϖRayl
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
    #@show τRayl, ϖRayl[1], ϖ
    ϖ += τRayl * ϖRayl[1]
    A += τRayl * ϖRayl[1]

    Z⁺⁺ = τRayl * ϖRayl[1] * Rayl𝐙⁺⁺
    Z⁻⁺ = τRayl * ϖRayl[1] * Rayl𝐙⁻⁺

    for i = 1:nAer
        #@show τ, ϖ , A, τAer[i]
        τ   += τAer[i]
        ϖ   += τAer[i] * aerosol_optics[i].ω̃
        A   += τAer[i] * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ)
        Z⁺⁺ += τAer[i] * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ) * Aer𝐙⁺⁺[:,:,i]
        Z⁻⁺ += τAer[i] * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ) * Aer𝐙⁻⁺[:,:,i]
        #@show τ, ϖ , A
    end
    
    Z⁺⁺ /= A
    Z⁻⁺ /= A
    A /= ϖ
    ϖ /= τ
    
    # Rescaling composite SSPs according to Eqs. A.3 of Sanghavi et al. (2013) or Eqs.(8) of Sanghavi & Stephens (2015)
    #@show τRayl, τ,A,  ϖ
    τ *= (FT(1) - (FT(1) - A) * ϖ)
    ϖ *= A / (FT(1) - (FT(1) - A) * ϖ)#Suniti
    #@show τRayl, τ
    fscattRayl = τRayl/τ
    # Adding absorption optical depth / albedo:
    τ_λ = τ_abs .+ τ    
    ϖ_λ = (τ * ϖ) ./ τ_λ
    
    return Array(τ_λ), Array(ϖ_λ), τ, ϖ, Array(Z⁺⁺), Array(Z⁻⁺), fscattRayl
end

"When performing RT_run, this function pre-calculates properties for all layers, before any Core RT is performed"
function construct_all_atm_layers(
        FT, nSpec, Nz, NquadN, 
        τRayl, τAer, aerosol_optics, 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, 
        τ_abs, 
        ϖ_Cabannes,
        arr_type, qp_μ, μ₀, m)

    FT_ext   = eltype(τAer)
    FT_phase = eltype(τAer)

    # Empty matrices to hold all values
    τ_λ_all   = zeros(FT_ext, nSpec, Nz)
    ϖ_λ_all   = zeros(FT_ext, nSpec, Nz)
    τ_all     = zeros(FT_ext, Nz)
    ϖ_all     = zeros(FT_ext, Nz)
    Z⁺⁺_all   = zeros(FT_phase, NquadN, NquadN, Nz)
    Z⁻⁺_all   = zeros(FT_phase, NquadN, NquadN, Nz)
    
    dτ_max_all  = zeros(FT_ext, Nz)
    dτ_all      = zeros(FT_ext, Nz)
    fscattRayl_all  =  zeros(FT_ext, Nz)
    ndoubl_all  = zeros(Int64, Nz)
    dτ_λ_all    = zeros(FT_ext, nSpec, Nz)
    expk_all    = zeros(FT_ext, nSpec, Nz)
    scatter_all = zeros(Bool, Nz)

    for iz=1:Nz
        
        # Construct atmospheric properties
        τ_λ_all[:, iz], 
        ϖ_λ_all[:, iz], 
        τ_all[iz], 
        ϖ_all[iz], 
        Z⁺⁺_all[:,:,iz], 
        Z⁻⁺_all[:,:,iz], 
        fscattRayl_all[iz] = construct_atm_layer(τRayl[iz], τAer[:,iz], 
            ϖ_Cabannes,
            aerosol_optics, 
            Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, 
            τ_abs[:,iz], arr_type)
        #@show fscattRayl_all[iz]
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

    return ComputedAtmosphereProperties(τ_λ_all, ϖ_λ_all, τ_all, ϖ_all, Z⁺⁺_all, Z⁻⁺_all, dτ_max_all, dτ_all, ndoubl_all, dτ_λ_all, expk_all, scatter_all, τ_sum_all, fscattRayl_all, scattering_interfaces_all)
end



"Given the CrossSectionModel, the grid, and the AtmosphericProfile, fill up the τ_abs array with the cross section at each layer
(using pressures/temperatures) from the profile" 
function compute_absorption_profile!(τ_abs::Array{FT,2}, 
                                     absorption_model, 
                                     grid,
                                     vmr,
                                     profile::AtmosphericProfile,
                                     ) where FT 

    # The array to store the cross-sections must be same length as number of layers
    @assert size(τ_abs,2) == length(profile.p_full)

    @showprogress 1 for iz in 1:length(profile.p_full)

        # Pa -> hPa
        p = profile.p_full[iz]
        T = profile.T[iz]

        # Either use the current layer's vmr, or use the uniform vmr
        vmr_curr = vmr isa AbstractArray ? vmr[iz] : vmr

        # Changed index order
        # @show iz,p,T,profile.vcd_dry[iz], vmr_curr
        #@show typeof(τ_abs), typeof(vmr_curr), typeof(profile.vcd_dry[iz]), typeof(p), typeof(T)
        #@show typeof(absorption_cross_section(absorption_model, grid, p, T))
        #temp = Array(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
        #@show minimum(temp), p, T, profile.vcd_dry[iz] * vmr_curr
        #@show iz, profile.vcd_dry[iz], vmr_curr, p, T
        τ_abs[:,iz] += Array(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
    end
    
end