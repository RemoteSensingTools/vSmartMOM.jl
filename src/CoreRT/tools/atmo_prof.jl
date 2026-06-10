#=

This file contains functions that are related to atmospheric profile calculations

=#

"""
    compute_atmos_profile_fields(T::AbstractArray{FT,1}, p_half::AbstractArray{FT,1}, q, vmr; g₀=9.807) -> Tuple

Computes atmospheric profile fields, including volume mixing ratios (VMR) of H2O, dry and wet volume column densities (VCDs), and layer thicknesses (Δz).

# Arguments
- `T::AbstractArray{FT,1}`: Temperature profile in Kelvin (K).
- `p_half::AbstractArray{FT,1}`: Pressure at half-levels in hectopascals (hPa).
- `q`: Specific humidity in grams of water vapor per kilogram of moist air (g/kg).
- `vmr`: Dictionary containing volume mixing ratios of various trace gases.
- `g₀=9.807`: Gravitational acceleration (m/s²), default is 9.807 m/s².

# Returns
- `p_full`: Pressure at full levels (hPa).
- `p_half`: Pressure at half levels (hPa), same as input.
- `vmr_h2o`: Volume mixing ratio of H2O (unitless).
- `vcd_dry`: Dry volume column density (molec/cm²).
- `vcd_h2o`: Wet volume column density (molec/cm²).
- `new_vmr`: Interpolated volume mixing ratios of trace gases (Dictionary).
- `Δz`: Layer thicknesses (m).

# Description
This function calculates various atmospheric profile fields given temperature, pressure, specific humidity, and initial volume mixing ratios of trace gases. It computes:
1. Pressure at full levels.
2. Volume mixing ratio of H2O from specific humidity.
3. Dry and wet volume column densities (VCDs).
4. Layer thicknesses (Δz).
5. Interpolated volume mixing ratios for other trace gases.
"""
function compute_atmos_profile_fields(T::AbstractArray{FT,1}, p_half::AbstractArray{FT,1}, q, vmr; g₀=9.8032465) where FT
    # g₀ default aligned with sanghavi (Bodhaine 1999 Eq.30 implementation
    # uses gravitational acceleration ≈ 9.8032 m/s² for column-air mass — see
    # sanghavi atmo_prof.jl). q is treated as kg/kg directly (volume mixing
    # ratio of water vapor by mass), matching sanghavi convention.
    #FT = eltype(T)
    Nₐ = FT(6.02214179e+23)
    R  = FT(8.3144598)
    # Calculate full pressure levels
    p_full = (p_half[2:end] + p_half[1:end-1]) / 2

    # Dry and wet mass
    dry_mass = FT(28.9644e-3)    # in kg/molec, weighted average for N2 and O2
    wet_mass = FT(18.01534e-3)   # just H2O
    ratio = dry_mass / wet_mass
    n_layers = length(T)

    # Also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )
    Δz      = zeros(FT, n_layers)
    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i]/(1-q[i]) * ratio# dry_mass/(dry_mass-wet_mass*(1-1/q[i]))
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd = Nₐ * Δp / (M  * g₀ * 100^2) * 100
        vcd_dry[i] = vmr_dry    * vcd   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * vcd
        Δz[i] =  (log(p_half[i + 1]) - log(p_half[i])) / (g₀ * M  / (R * T[i]) )
        #@show Δz, T[i], M, Δp
    end

    # TODO: This is still a bit clumsy:
    new_vmr = Dict{String, Union{Real, Vector}}()

    for molec_i in keys(vmr)
        if vmr[molec_i] isa AbstractArray
            if length(vmr[molec_i]) == length(p_full)
                new_vmr[molec_i] = vmr[molec_i]
            else
                @info "Warning, make sure that the VMR is interpolated correctly! Right now, it might be tricky"
                pressure_grid = collect(range(minimum(p_full), maximum(p_full), length=length(vmr[molec_i])))
                interp_linear = LinearInterpolation(pressure_grid, vmr[molec_i])
                new_vmr[molec_i] = [interp_linear(x) for x in p_full]
            end
        else
            new_vmr[molec_i] = vmr[molec_i]
        end
    end

    return p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz

end

## IO reader methods moved to src/IO/AtmosProfile.jl to decouple CoreRT from IO.

"""
    reduce_profile(n, profile; binavg=false)

Reduce an `AtmosphericProfile` to `n` layers.

Default (binavg=false): linear interpolation onto uniform pressure half-levels
spanning [profile.p_half[1], profile.p_half[end]]. VCDs are recomputed from the
new pressure grid (consistent with `compute_atmos_profile_fields`), Δz is
re-derived via the hydrostatic relation for the new T and vmr_h2o. This is the
physics-forward default: it preserves profile shape across the full column
rather than averaging within coarse bins.

Pass `binavg=true` (or call `reduce_profile_binavg`) for the legacy bin-averaging
method.
"""
function reduce_profile(n::Int, profile::AtmosphericProfile{FT}; binavg::Bool=false) where {FT}

    if binavg
        return reduce_profile_binavg(n, profile)
    end

    @assert n < length(profile.T)

    (; vmr) = profile

    Nₐ       = FT(6.02214179e+23)
    R        = FT(8.3144598)
    dry_mass = FT(28.9644e-3)
    wet_mass = FT(18.01534e-3)
    g₀       = FT(9.8032465)  # aligned with sanghavi compute_atmos_profile_fields

    # New uniform half-levels spanning the original column extent (TOA → BOA)
    p_half = collect(range(profile.p_half[1], profile.p_half[end], length=n+1))
    p_full = (p_half[1:end-1] .+ p_half[2:end]) ./ FT(2)

    # Linear interpolation on the profile's full-level pressure grid
    old_p = profile.p_full
    function _interp(data::AbstractVector)
        grid = collect(range(minimum(old_p), maximum(old_p), length=length(data)))
        itp  = LinearInterpolation(grid, data)
        return FT.(itp.(p_full))
    end

    T       = _interp(profile.T)
    q       = _interp(profile.q)
    vmr_h2o = _interp(profile.vmr_h2o)

    # Recompute VCDs and Δz from the new layers (consistent with compute_atmos_profile_fields)
    vcd_dry = zeros(FT, n)
    vcd_h2o = zeros(FT, n)
    Δz_     = zeros(FT, n)
    for i = 1:n
        Δp      = p_half[i+1] - p_half[i]
        vmr_dry = FT(1) - vmr_h2o[i]
        M       = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd     = Nₐ * Δp / (M * g₀ * FT(100)^2) * FT(100)
        vcd_dry[i] = vmr_dry    * vcd
        vcd_h2o[i] = vmr_h2o[i] * vcd
        Δz_[i]  = (log(p_half[i+1]) - log(p_half[i])) / (g₀ * M / (R * T[i]))
    end

    # Interpolate per-species VMR profiles (fallback to scalar passthrough)
    new_vmr = Dict{String, Union{Real, Vector}}()
    for molec_i in keys(vmr)
        if profile.vmr[molec_i] isa AbstractArray
            new_vmr[molec_i] = _interp(profile.vmr[molec_i])
        else
            new_vmr[molec_i] = profile.vmr[molec_i]
        end
    end

    return AtmosphericProfile(T, p_full, q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz_)
end

"Legacy bin-averaging profile reduction (opt-in via `reduce_profile(n, p; binavg=true)`)."
function reduce_profile_binavg(n::Int, profile::AtmosphericProfile{FT}) where {FT}

    # Can only reduce the profile, not expand it
    @assert n < length(profile.T)

    # Unpack the profile vmr
    (; vmr, Δz) = profile

    # New rough half levels (boundary points)
    a = range(0, maximum(profile.p_half), length=n+1)

    # Matrices to hold new values
    T = zeros(FT, n);
    q = zeros(FT, n);
    Δz_ = zeros(FT, n);
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
        # Set the pressure levels accordingly
        p_half[i]   = a[i]
        p_half[i+1] = a[i+1]

        # Re-average the other parameters to produce new layers
        p_full[i] = mean(profile.p_full[ind])
        T[i] = mean(profile.T[ind])
        q[i] = mean(profile.q[ind])
        Δz_[i] = sum(Δz[ind])
        vcd_dry[i] = sum(profile.vcd_dry[ind])
        vcd_h2o[i] = sum(profile.vcd_h2o[ind])
        vmr_h2o[i] = vcd_h2o[i]/vcd_dry[i]
    end

    new_vmr = Dict{String, Union{Real, Vector}}()

    # TODO: This needs a VCD_dry weighted average!
    for molec_i in keys(vmr)
        if profile.vmr[molec_i] isa AbstractArray
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
    # TODO: reduce_profile and getRayleighLayerOptProp both deserve refactoring
    # beyond this merge. The current defaults are physics-forward (Bodhaine
    # 1999 Eq. 30, interpolated profile) but the implementation layout mixes
    # legacy and modern code paths via keyword args and name suffixes. Cleaner
    # architecture: a single configurable ProfileReduction strategy (dispatch on
    # type) and a single RayleighFormula strategy, both YAML-configurable.
    # Tracked as a future PR; not a merge blocker.
    Nz = length(vcd_dry)
    τRayl = zeros(FT, size(λ,1), Nz)
    # Bodhaine 1999 "On Rayleigh optical depth calculations" Eq. 30.
    # Has an implicit depolarization ratio ρ₀ ≈ 0.0279 (see Bodhaine Table 3);
    # we rescale to the caller-supplied depol_fct for flexibility.
    tau_scat = FT(0.002152) .* (FT(1.0455996) .- FT(341.29061) .* λ.^(-2) .- FT(0.90230850) .* λ.^2) ./
               (FT(1) .+ FT(0.0027059889) .* λ.^(-2) .- FT(85.968563) .* λ.^2)
    tau_scat = tau_scat .* (psurf / FT(1013.25))
    ρ₀ = FT(0.0279)
    tau_scat = tau_scat .* (FT(6) - FT(7)*ρ₀) * (FT(6) + FT(3)*depol_fct) /
                          ((FT(6) + FT(3)*ρ₀) * (FT(6) - FT(7)*depol_fct))
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
    $(FUNCTIONNAME)(total_τ, dist, profile)
    
Returns the aerosol optical depths per layer using a Distribution function in p
"""
function getAerosolLayerOptProp(total_τ::FT, dist::Distribution, profile::AtmosphericProfile) where FT
    (; p_half, p_full, Δz) = profile
    
    ρ = pdf.(dist,p_full) .* Δz
    τAer  =  (total_τ / sum(ρ)) * ρ
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
    @assert length(vmr) ==1 || length(vmr) == length(profile.p_full)  "Length of VMR array has to match profile size or be uniform"
    #@show grid
    @showprogress 1 for iz in 1:length(profile.p_full)

        # Pa -> hPa
        p = profile.p_full[iz]
        T = profile.T[iz]
        # Either use the current layer's vmr, or use the uniform vmr
        vmr_curr = vmr isa AbstractArray ? vmr[iz] : vmr
        #@show vmr_curr
        τ_abs[:,iz] += collect(absorption_cross_section(absorption_model, grid, p, T)) * profile.vcd_dry[iz] * vmr_curr
    end
    
end
