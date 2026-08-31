#=
Linearized atmospheric profile functions.
Shared functions (reduce_profile, getRayleighLayerOptProp, construct_atm_layer)
are defined in atmo_prof.jl; this file contains only the lin-mode extensions.
=#

"""Profile tangents for a fixed grid whose bottom pressure boundary is p_surf."""
function psurf_profile_tangents(profile::AtmosphericProfile{FT}; g₀=FT(9.8032465)) where FT
    n = length(profile.T)
    vcd_dry_dot = zeros(FT, n)
    vcd_h2o_dot = zeros(FT, n)
    Δz_dot = zeros(FT, n)
    dry_mass = FT(28.9644e-3)
    wet_mass = FT(18.01534e-3)
    Nₐ = FT(6.02214179e23)
    R = FT(8.3144598)
    vmr_h2o = profile.vmr_h2o[end]
    x_dry = inv(one(FT) + vmr_h2o)
    x_h2o = vmr_h2o * x_dry
    M = x_dry * dry_mass + x_h2o * wet_mass
    vcd_dot = Nₐ / (M * g₀ * FT(100)^2) * FT(100)
    vcd_dry_dot[end] = x_dry * vcd_dot
    vcd_h2o_dot[end] = x_h2o * vcd_dot
    Δz_dot[end] = R * profile.T[end] / (g₀ * M * profile.p_half[end])
    return (; vcd_dry_dot, vcd_h2o_dot, Δz_dot)
end

"""
Return a normalized column fraction and its tangent without squaring the
column total. The factored quotient rule is algebraically equivalent to
`(x_dot*sum(x) - x*sum(x_dot))/sum(x)^2`, but remains finite in `Float32` for
atmospheric molecular columns of order `1e25`.
"""
function _normalized_column_fraction_tangent(column, column_dot)
    total = sum(column)
    isfinite(total) && total > zero(total) || throw(ArgumentError(
        "column total must be finite and positive"))
    fraction = column ./ total
    fraction_dot = (column_dot .- fraction .* sum(column_dot)) ./ total
    return fraction, fraction_dot
end

"Derivative of a normalized aerosol allocation with respect to p_surf."
function aerosol_profile_psurf_tangent(dist::Normal, profile::AtmosphericProfile,
                                        Δz_dot)
    p = profile.p_full
    ρ = pdf.(dist, p) .* profile.Δz
    ρdot = zero(ρ)
    pdot = one(eltype(p)) / 2
    ρdot[end] = pdf(dist, p[end]) * Δz_dot[end] +
                ρ[end] * (-(p[end] - dist.μ) / dist.σ^2) * pdot
    return (ρdot .- ρ .* (sum(ρdot) / sum(ρ))) ./ sum(ρ)
end

function aerosol_profile_psurf_tangent(dist::LogNormal, profile::AtmosphericProfile,
                                        Δz_dot)
    z_half = half_level_altitudes(profile)
    FT = eltype(z_half)
    z_half_dot = zeros(FT, length(z_half))
    for i in length(profile.Δz):-1:1
        z_half_dot[i] = z_half_dot[i + 1] + Δz_dot[i] / FT(1000)
    end
    cdf_slope(z) = z <= zero(FT) ? zero(FT) : pdf(dist, z)
    ρ = [_altitude_lognormal_cdf(z_half[i], dist) -
         _altitude_lognormal_cdf(z_half[i + 1], dist)
         for i in eachindex(profile.Δz)]
    ρdot = [cdf_slope(z_half[i]) * z_half_dot[i] -
            cdf_slope(z_half[i + 1]) * z_half_dot[i + 1]
            for i in eachindex(profile.Δz)]
    return (ρdot .- ρ .* (sum(ρdot) / sum(ρ))) ./ sum(ρ)
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

"Apply the quotient rule to two tangents of an unnormalized aerosol profile."
function _normalized_aerosol_profile_tangents(total_τ, dist::Distribution,
                                               profile::AtmosphericProfile,
                                               ρ, dρ₁, dρ₂)
    # Obtain the base column from the production forward helper itself. This is
    # deliberately not reconstructed from distribution moments: LinMode must
    # retain bit-exact base-state parity for every supported Distribution.
    τAer = getAerosolLayerOptProp(total_τ, dist, profile)
    norm_ρ = sum(ρ)
    scale = total_τ / norm_ρ
    dτ₁ = scale .* (dρ₁ .- ρ .* (sum(dρ₁) / norm_ρ))
    dτ₂ = scale .* (dρ₂ .- ρ .* (sum(dρ₂) / norm_ρ))
    return τAer, dτ₁, dτ₂
end

"""
    getAerosolLayerOptProp(lin::LinMode, total_τ, dist, profile)

Return the exact production forward aerosol column and two profile-shape
tangents. The tangent convention depends on the distribution:

- `Normal(μ, σ)`: derivatives with respect to `μ` (`p₀`) and `σ` (`σp`).
- `LogNormal(log(z₀), σ₀)`: derivatives with respect to the user-facing
  altitude parameters `z₀ = exp(μ)` and `σ₀ = σ`.

Other `Distribution` types throw an `ArgumentError` because the fixed
two-column aerosol layout does not define a safe parameter mapping for them.
"""
function getAerosolLayerOptProp(::LinMode, total_τ, dist::Distribution,
                                profile::AtmosphericProfile)
    throw(ArgumentError(
        "LinMode aerosol-profile tangents support Normal and LogNormal; " *
        "got $(typeof(dist))"))
end

function getAerosolLayerOptProp(::LinMode, total_τ, dist::Normal,
                                profile::AtmosphericProfile)
    (; p_full, Δz) = profile
    ρ = pdf.(dist, p_full) .* Δz
    offset = p_full .- dist.μ
    dρ_dμ = ρ .* offset ./ dist.σ^2
    dρ_dσ = ρ .* (offset.^2 ./ dist.σ^3 .- inv(dist.σ))
    return _normalized_aerosol_profile_tangents(
        total_τ, dist, profile, ρ, dρ_dμ, dρ_dσ)
end

function getAerosolLayerOptProp(::LinMode, total_τ, dist::LogNormal,
                                profile::AtmosphericProfile)
    z_half = half_level_altitudes(profile)
    z₀ = exp(dist.μ)
    FT = eltype(z_half)
    function cdf_and_tangents(z)
        z <= zero(FT) && return (zero(FT), zero(FT), zero(FT))
        a = (log(z) - dist.μ) / dist.σ
        φ = exp(-a^2 / FT(2)) / sqrt(FT(2π))
        F = _altitude_lognormal_cdf(z, dist)
        return F, -φ / (dist.σ * z₀), -φ * a / dist.σ
    end
    vals = cdf_and_tangents.(z_half)
    F = first.(vals)
    dF_dz₀ = getindex.(vals, 2)
    dF_dσ₀ = getindex.(vals, 3)
    ρ = F[1:end-1] .- F[2:end]
    dρ_dz₀ = dF_dz₀[1:end-1] .- dF_dz₀[2:end]
    dρ_dσ₀ = dF_dσ₀[1:end-1] .- dF_dσ₀[2:end]
    return _normalized_aerosol_profile_tangents(
        total_τ, dist, profile, ρ, dρ_dz₀, dρ_dσ₀)
end

"""
    getAerosolLayerOptProp(lin::LinMode, total_τ, p₀, σp, profile)

Backward-compatible pressure-profile entry point. This is equivalent to the
`Normal(p₀, σp)` distribution dispatch above.
"""
function getAerosolLayerOptProp(lin::LinMode, total_τ, p₀, σp,
                                profile::AtmosphericProfile)
    return getAerosolLayerOptProp(lin, total_τ, Normal(p₀, σp), profile)
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
                                    ;
                                    self_broadener_vmr=nothing,
                                    ) where FT 

    # The array to store the cross-sections must be same length as number of layers
    @assert size(τ_abs,2) == length(profile.p_full)

    @showprogress 1 for iz in 1:length(profile.p_full)

        # Pa -> hPa
        p = profile.p_full[iz]
        T = profile.T[iz]

        # Either use the current layer's vmr, or use the uniform vmr
        vmr_curr = vmr isa AbstractArray ? vmr[iz] : vmr
        broadener_curr = self_broadener_vmr isa AbstractArray ?
            self_broadener_vmr[iz] : self_broadener_vmr

        # Changed index order
        # @show iz,p,T,profile.vcd_dry[iz], vmr_curr
        #@show typeof(τ_abs), typeof(vmr_curr), typeof(profile.vcd_dry[iz]), typeof(p), typeof(T)
        #@show typeof(absorption_cross_section(absorption_model, grid, p, T))
        #temp = collect(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
        #@show minimum(temp), p, T, profile.vcd_dry[iz] * vmr_curr
        #@show iz, profile.vcd_dry[iz], vmr_curr, p, T
        # An explicit kwarg wins (LBL self-broadening); otherwise the model
        # type decides whether the profile supplies a broadener (ABSCO only)
        # — same `_default_broadener` rule as the forward path in atmo_prof.jl.
        broadener_curr = broadener_curr === nothing ?
            _default_broadener(absorption_model, profile, iz) : broadener_curr
        σ = collect(_layer_absorption_cross_section(
            absorption_model, grid, p, T, broadener_curr))
        τ_abs[:,iz] += σ * profile.vcd_dry[iz] * vmr_curr
        # Species-major flattened ordering: (igas - 1) * Nz + iz.
        # Keeping all other layers zero yields dτ(z)/dVMR(igas, iz).
        gas_layer_idx = (jac_idx - 1) * length(profile.p_full) + iz
        τ̇_abs[gas_layer_idx,:,iz] = σ * profile.vcd_dry[iz]
    end
    
end

"""
    compute_h2o_absorption_profile!(τ_abs, τ̇_abs, jac_idx,
                                    absorption_model, grid, profile)

Linearized H₂O line-absorption accumulation using the same layerwise moist
H₂O mole fraction for self broadening as the forward path.
"""
function compute_h2o_absorption_profile!(τ_abs::Array{FT,2},
                                         τ̇_abs::Array{FT,3},
                                         jac_idx::Integer,
                                         absorption_model,
                                         grid,
                                         profile::AtmosphericProfile) where FT
    x_h2o = _h2o_moist_mole_fraction.(profile.vmr_h2o)
    return compute_absorption_profile!(
        τ_abs, τ̇_abs, jac_idx, absorption_model, grid,
        profile.vmr_h2o, profile; self_broadener_vmr=x_h2o)
end
