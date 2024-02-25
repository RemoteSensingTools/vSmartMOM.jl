#=

This file contains functions that are related to atmospheric profile calculations

=#

"Compute pressure levels, vmr, vcd for atmospheric profile, given p_half, T, q"
function lin_compute_atmos_profile_fields(
            T::AbstractArray{FT,1}, 
            p_half::AbstractArray{FT,1}, 
            q, vmr,#,
            x;#,
            #dVMR_CO2,
            #dVMR_H2O; 
            g₀=9.807) where FT
    #@show "Atmos",  FT 
    # Floating type to use
    #FT = eltype(T)
    Nₐ = FT(6.02214179e+23)
    R  = FT(8.3144598)

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
    Δz      = zeros(FT, n_layers)
    z       = zeros(FT, n_layers)

    psurf = x[1] 
    @assert x[1]==p_half[end]
    # Now actually compute the layer VCDs
    M = FT(0.0)
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        a = (i<=65) ? x[8] : x[9]
        vmr_h2o[i] = a*dry_mass/(dry_mass-wet_mass*(1-1/q[i]))
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd = Nₐ * Δp / (M  * g₀ * 100^2) * 100
        vcd_dry[i] = vmr_dry    * vcd   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * vcd
        Δz[i] =  (log(p_half[i + 1]/p_half[i])) / (g₀ * M  / (R * T[i]) )
        z[1:i] = z[1:i] .+ Δz[i]#@show Δz, T[i], M, Δp
    end
    #Δp_surf = p_half[end] - p_half[end-1]
    dΔz0dpsurf = (1 ./ (p_half[end])) / (g₀ * M  / (R * T[end]) )
    dzdpsurf = zeros(length(z)) .+ dΔz0dpsurf
    
    prof  = LogNormal(x[6], x[7])
    vmr["CO2"] = (x[2].+zeros(length(z))) + 
                 (x[3] * exp.(-z./x[5])) +
                 (x[4] * pdf.(prof, z))
    vmr_co2 = vmr["CO2"] 
    #dVMR_H2O[1,:] = 0.0
    #dVMR_H2O[1,end] = dVMR_H2O[end]./Δp_surf
    dVMR_H2O = zeros(2, length(z))
    dVMR_CO2 = zeros(7, length(z))
    dVMR_H2O[1,:] = [vcd_h2o[1:65]/x[8]; vcd_h2o[66:end] * 0.0;] # wrt x[7]
    dVMR_H2O[2,:] = [vcd_h2o[1:65] * 0.0; vcd_h2o[66:end]/x[9];] # wrt x[8]

    dVMR_CO2[1,:] = (x[3] * exp.(-z./x[5]) * (-1/x[5]) .-
                    (pdf.(prof,z)./z) .* (1 .+ log.(z)/x[7]^2)) .* dzdpsurf; # wrt x[1] 
    dVMR_CO2[2,:] = 1.0 .+ zeros(length(z)) # wrt x[2]
    dVMR_CO2[3,:] = exp.(-z./x[5]) # wrt x[3]
    dVMR_CO2[4,:] = pdf.(prof, z) # wrt x[4]
    dVMR_CO2[5,:] = x[3] * exp.(-z./x[5]) .* z./(x[5])^2 # wrt x[5]
    dVMR_CO2[6,:] = x[4] * pdf.(prof, z) .* (log.(z) .- x[6]) / x[7]^2
    dVMR_CO2[7,:] = (x[4] * pdf.(prof, z) / x[7]) .* 
                        (((log.(z) .- x[6]) / x[7]).^2 .- 1)


    #=
    # TODO: This is still a bit clumsy:
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
    =#
    #return p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz, z
    return p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr_co2, Δz, z, dzdpsurf, dVMR_H2O, dVMR_CO2
end
#=
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
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, Δz = compute_atmos_profile_fields(T, p_half, q, Dict())
    elseif ("q" in keys(params_dict))
        p_half = convert(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q      = convert.(Float64, params_dict["q"])
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, Δz = compute_atmos_profile_fields(T, p_half, q, Dict())
    else
        p_half = convert.(Float64, params_dict["p_half"])
        psurf = p_half[end]
        q = zeros(length(T))
        p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, Δz = compute_atmos_profile_fields(T, p_half, q, Dict())
    end

    # Convert vmr to appropriate type
    vmr = convert(Dict{String, Union{Real, Vector}}, params_dict["vmr"])

    # Return the atmospheric profile struct
    return AtmosphericProfile(T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr)

end
=#
"Reduce profile dimensions by re-averaging to near-equidistant pressure grid"
function lin_reduce_profile(n::Int, linprofile::linAtmosphericProfile{FT}) where {FT}

    # Can only reduce the profile, not expand it
    @assert n < length(linprofile.T)

    # Unpack the profile vmr
    #@unpack vmr, Δz = linprofile
    @unpack Δz = linprofile

    # New rough half levels (boundary points)
    a = range(0, maximum(linprofile.p_half), length=n+1)

    # Matrices to hold new values
    T = zeros(FT, n);
    q = zeros(FT, n);
    Δz_ = zeros(FT, n);
    p_full = zeros(FT, n);
    p_half = zeros(FT, n+1);
    z = zeros(FT, n);
    vmr_h2o  = zeros(FT, n);
    vmr_co2  = zeros(FT, n);
    vcd_dry  = zeros(FT, n);
    vcd_h2o  = zeros(FT, n);
    dzdpsurf = zeros(FT, n);
    dVMR_H2O = zeros(FT, 2, n); 
    dVMR_CO2 = zeros(FT, 7, n);
    
    # Loop over target number of layers
    indices = []
    for i = 1:n

        # Get the section of the atmosphere with the i'th section pressure values
        ind = findall(a[i] .< linprofile.p_full .<= a[i+1]);
        push!(indices, ind)
        @assert length(ind) > 0 "Profile reduction has an empty layer"
        #@show i, ind, a[i], a[i+1]
        # Set the pressure levels accordingly
        p_half[i]   = a[i]   # profile.p_half[ind[1]]
        p_half[i+1] = a[i+1] # profile.p_half[ind[end]]

        # Re-average the other parameters to produce new layers

        
        p_full[i] = mean(linprofile.p_full[ind])
        T[i] = mean(linprofile.T[ind])
        q[i] = mean(linprofile.q[ind])
        Δz_[i] = sum(Δz[ind])
        z[i] = maximum(linprofile.z[ind])
        vcd_dry[i] = sum(linprofile.vcd_dry[ind])
        vcd_h2o[i] = sum(linprofile.vcd_h2o[ind])
        vmr_h2o[i] = sum(linprofile.vmr_h2o[ind].*linprofile.p_half[ind]./linprofile.T[ind])/
                sum(linprofile.p_half[ind]./linprofile.T[ind])#vcd_h2o[i]/vcd_dry[i]
        vmr_co2[i] = sum(linprofile.vmr_co2[ind].*linprofile.p_half[ind]./linprofile.T[ind])/
                sum(linprofile.p_half[ind]./linprofile.T[ind])
        dzdpsurf[i] = mean(linprofile.dzdpsurf[ind])
        for j=1:2
            dVMR_H2O[j,i] = sum(linprofile.dVMR_H2O[j,ind].*linprofile.p_half[ind]./linprofile.T[ind])/
            sum(linprofile.p_half[ind]./linprofile.T[ind])
            dVMR_CO2[j,i] = sum(linprofile.dVMR_CO2[j,ind].*linprofile.p_half[ind]./linprofile.T[ind])/
                    sum(linprofile.p_half[ind]./linprofile.T[ind]) 
        end
        for j=3:7
            dVMR_CO2[j,i] = sum(linprofile.dVMR_CO2[j,ind].*linprofile.p_half[ind]./linprofile.T[ind])/
                    sum(linprofile.p_half[ind]./linprofile.T[ind]) 
        end
    end
    #@show indices
#=
    new_vmr = Dict{String, Union{Real, Vector}}()

    # need to double check this logic, maybe better to add VCDs?!
    for molec_i in keys(vmr)
        if profile.vmr[molec_i] isa AbstractArray
            # TODO: This needs a VCD_dry weighted average!
            new_vmr[molec_i] = [mean(profile.vmr[molec_i][ind]) for ind in indices]
        else
            new_vmr[molec_i] = profile.vmr[molec_i]
        end
    end
=#
    return linAtmosphericProfile(T, p_full, q, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr_co2, Δz_, z, dzdpsurf, dVMR_H2O, dVMR_CO2)
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
function getRayleighLayerOptProp_lin(psurf::FT, λ::Union{Array{FT}, FT}, depol_fct::FT, vcd_dry::Array{FT}) where FT
    # TODO: Use noRS/noRS_plus to use n2/o2 molecular constants
    # to compute tau_scat and depol_fct
    Nz = length(vcd_dry)
    τRayl = zeros(FT,size(λ,1),Nz)
    lin_τRayl = zeros(FT,size(λ,1),Nz) # derivative of τRayl wrt psurf
    # Total vertical Rayleigh scattering optical thickness, TODO: enable sub-layers and use VCD based taus
    tau_scat = FT(0.00864) * (psurf / FT(1013.25)) *  λ.^(-FT(3.916) .- FT(0.074) * λ .- FT(0.05) ./ λ)  
    tau_scat = tau_scat * (FT(6.0) + FT(3.0) * depol_fct) / (FT(6.0)- FT(7.0) * depol_fct) 
    # @show tau_scat, λ
    k = tau_scat / sum(vcd_dry)
    for i = 1:Nz
        τRayl[:,i] .= k * vcd_dry[i]
        lin_τRayl[:,i] .= τRayl[:,i]/psurf 
    end 
    return τRayl, lin_τRayl
end


"""
    $(FUNCTIONNAME)(total_τ, p₀, σp, p_half)
    
Returns the aerosol optical depths per layer using a Gaussian distribution function with p₀ and σp on a pressure grid
"""

function getAerosolLayerOptProp_lin(total_τ, z₀, σz, z, dzdpsurf)#, p_half)

    # Need to make sure we can also differentiate wrt σp (FT can be Dual!)
    FT = eltype(z₀)
    Nz = length(z)
    #ρ = zeros(FT,Nz)
    #dρdz₀ = zeros(FT,Nz)
    #dρdσz = zeros(FT,Nz)
    # @show p_half, p₀, σp

    prof = LogNormal(log(z₀), σz)
    τ_aer = total_τ * pdf.(prof, z)
    lin_τ_aer_psurf =  - τ_aer./z .* 
            (1 .+ log.(z)/σz^2) .* dzdpsurf
    lin_τ_aer_z₀ = τ_aer .* (log.(z) .- log(z₀)) / σz^2
    lin_τ_aer_σz = (τ_aer / σz) .* 
                        (((log.(z) .- log(z₀)) / σz).^2 .- 1)

    # return convert(FT, τ_aer, lin_τ_aer_psurf, lin_τ_aer_z₀, lin_τ_aer_σz)
    return τ_aer, lin_τ_aer_psurf, lin_τ_aer_z₀, lin_τ_aer_σz;

end



#=
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
    @show ϖRayl
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
=#

# TODO:
"Given the CrossSectionModel, the grid, and the AtmosphericProfile, fill up the τ_abs array with the cross section at each layer
(using pressures/temperatures) from the profile" 
function compute_absorption_profile_lin!(τ_abs::Array{FT,2},
                                     lin_τ_abs::Array{FT,3},
                                     Δp_surf,
                                     dVMR,
                                     #dVMR_CO2,
                                     absorption_model, 
                                     grid,
                                     vmr,
                                     profile::linAtmosphericProfile,
                                     ) where FT 

    # The array to store the cross-sections must be same length as number of layers
    @assert size(τ_abs,2) == length(profile.p_full)
    @assert length(vmr) ==1 || length(vmr) == length(profile.p_full)  "Length of VMR array has to match profile size or be uniform"
    #@show grid
    @showprogress 1 for iz in 1:length(profile.p_full)

        # Pa -> hPa
        p = profile.p_full[iz]
        T = profile.T[iz]
        # Either use the current layer's vmr, or use the uniform vmr
        vmr_curr = vmr isa AbstractArray ? vmr[iz] : vmr
        Δτ = Array(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
        τ_abs[:,iz] += Δτ   # Array(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
        
        for ipar in 1:9
            lin_τ_abs[ipar,:,iz] += Δτ * (dVMR[ipar,iz]./vmr_curr)            
        end
        if iz==length(profile.p_full)
            lin_τ_abs[1,:,iz] += Δτ/Δp_surf  
        end
    end
    
end