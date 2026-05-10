#=
Linearized atmospheric profile functions.
Shared functions (reduce_profile, getRayleighLayerOptProp, construct_atm_layer)
are defined in atmo_prof.jl; this file contains only the lin-mode extensions.
=#

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
        #temp = collect(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
        #@show minimum(temp), p, T, profile.vcd_dry[iz] * vmr_curr
        #@show iz, profile.vcd_dry[iz], vmr_curr, p, T
        τ_abs[:,iz] += collect(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
        τ̇_abs[jac_idx,:,iz] = collect(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz]
    end
    
end
