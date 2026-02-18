#=

This file contains functions that are related to atmospheric profile calculations

=#

"Reduce profile dimensions by re-averaging to near-equidistant pressure grid"
function reduce_profile(n::Int, profile::AtmosphericProfile{FT}) where {FT}

    # Can only reduce the profile, not expand it
    @assert n < length(profile.T)

    # Unpack the profile vmr
    @unpack vmr, Δz = profile

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
    Δz_ = zeros(FT, n);

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
        Δz_[i] = sum(Δz[ind])
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

    return AtmosphericProfile(T, p_full, q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz_)
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
function getAerosolLayerOptProp(total_τ, z₀, σ₀, p_half, T)
    FT = eltype(T[1])
    R  = FT(8.3144598) # J/mol.K
    g₀ = 9.807 # m/s^2
    M₀ = FT(28.9644e-3) #kg/mol
    H = R*T/(M₀*g₀)
    Nz = length(p_half)-1
    dz = zeros(Nz)
    z = zeros(Nz)
    dz .= H.*log.(p_half[2:end]./p_half[1:end-1])
    dz .*= 1.e-3 #m->km
    z[end] = 0.0#dz[end]./2
    for i=Nz-1:-1:1
        z[i] = z[i+1]+dz[i+1]#(dz[i+1]+dz[i])./2 #this has been done to prevent dz=Inf resulting from p_half[1]=0
    end
    prof = LogNormal(log(z₀), σ₀)
    τAer = total_τ * pdf.(prof, z)
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

"""
Returns the aerosol optical depths per layer using a Gaussian distribution function with p₀ and σp on a pressure grid
"""
function getAerosolLayerOptProp(lin::LinMode, total_τ, z₀, σ₀, p_half, T)
    FT = eltype(T[1])
    #prof = LogNormal(log(z₀), σ₀)
    R  = FT(8.3144598) # J/mol.K
    g₀ = 9.807 # m/s^2
    M₀ = FT(28.9644e-3) #kg/mol
    H = R*T/(M₀*g₀)
    Nz = length(p_half)-1
    dz = zeros(Nz)
    z = zeros(Nz)
    dz .= H.*log.(p_half[2:end]./p_half[1:end-1])
    dz .*= 1.e-3 #m->km
    z[end] = 1.e-6 #0.0#dz[end]./2
    for i=Nz-1:-1:1
        z[i] = z[i+1]+dz[i+1]#(dz[i+1]+dz[i])./2 #this has been done to prevent dz=Inf resulting from p_half[1]=0
        #@show i, z[i]
        #    τAer[i+1] = total_τ * (cdf(prof, z[i]) - cdf(prof, z[i+1])) #pdf.(prof, z)
    end
    #τAer[1] = total_τ * (1.0 - cdf(prof, z[1]))
    @assert all(z .>= 0) "z must be strictly positive"

    # prepare u and phi(u)
    u = (log.(z) .- log(z₀)) ./ σ₀
    φ = pdf.(Normal(), u)             # standard normal pdf at u
    F = cdf.(Normal(), u)             # F(z) = Φ(u)

    # derivatives of F w.r.t parameters
    dF_dz₀ = - φ ./ (σ₀ .* z₀)                # ∂F/∂z0
    dF_dσ₀ = - φ .* (log.(z) .- log(z₀)) ./ (σ₀.^2)  # ∂F/∂σ0

    Nz = length(z)
    τAer   = zeros(eltype(z), Nz)
    dτdz₀  = similar(τAer)
    dτdσ₀  = similar(τAer)

    # vectorized construction (no explicit loop required)
    if Nz >= 2
        τAer[2:Nz]  .= total_τ .* (F[1:end-1] .- F[2:end])
        dτdz₀[2:Nz] .= total_τ .* (dF_dz₀[1:end-1] .- dF_dz₀[2:end])
        dτdσ₀[2:Nz] .= total_τ .* (dF_dσ₀[1:end-1] .- dF_dσ₀[2:end])
    end

    # top layer
    τAer[1]   = total_τ * (1.0 - F[1])
    dτdz₀[1]  = - total_τ * dF_dz₀[1]
    dτdσ₀[1]  = - total_τ * dF_dσ₀[1]

    return convert.(FT, τAer), convert.(FT, dτdz₀), convert.(FT, dτdσ₀)
end

"""
    getAerosolLayerOptProp(lin::LinMode, total_τ, p₀, σp, p_half)

Pressure-based aerosol vertical profile with analytic Jacobians w.r.t. p₀ (layer center 
pressure) and σp (layer width in pressure).  Returns `(τAer, dτ_dp₀, dτ_dσp)`.
Uses a Gaussian in pressure space, normalized so that `sum(τAer) == total_τ`.
"""
function getAerosolLayerOptProp(lin::LinMode, total_τ, p₀, σp, p_half)
    FT = eltype(p₀)
    Nz = length(p_half) - 1
    ρ      = zeros(FT, Nz)
    dρ_dp₀ = zeros(FT, Nz)
    dρ_dσp = zeros(FT, Nz)

    for i = 1:Nz
        dp = p_half[i+1] - p_half[i]
        p  = (p_half[i+1] + p_half[i]) / 2
        gauss = (1 / (σp * sqrt(2π))) * exp(-(p - p₀)^2 / (2σp^2))
        ρ[i] = gauss * dp
        # ∂ρ/∂p₀  = ρ * (p-p₀)/σp²
        dρ_dp₀[i] = ρ[i] * (p - p₀) / σp^2
        # ∂ρ/∂σp  = ρ * ((p-p₀)²/σp³ - 1/σp)
        dρ_dσp[i] = ρ[i] * ((p - p₀)^2 / σp^3 - 1 / σp)
    end

    Norm   = sum(ρ)
    S_dp₀  = sum(dρ_dp₀)
    S_dσp  = sum(dρ_dσp)

    τAer   = (total_τ / Norm) .* ρ
    dτ_dp₀ = (total_τ / Norm) .* (dρ_dp₀ .- ρ .* (S_dp₀ / Norm))
    dτ_dσp = (total_τ / Norm) .* (dρ_dσp .- ρ .* (S_dσp / Norm))

    return convert.(FT, τAer), convert.(FT, dτ_dp₀), convert.(FT, dτ_dσp)
end

"Given the CrossSectionModel, the grid, and the AtmosphericProfile, fill up the τ_abs array with the cross section at each layer
(using pressures/temperatures) from the profile" 
function compute_absorption_profile!(τ_abs::Array{FT,2}, 
                                    τ̇_abs::Array{FT,3}, 
                                    jac_idx::Integer,
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
        τ̇_abs[jac_idx,:,iz] = Array(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz]
    end
    
end
