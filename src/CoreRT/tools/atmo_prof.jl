#=

This file contains functions that are related to atmospheric profile calculations

=#

"""
    compute_atmos_profile_fields(T::AbstractArray{FT,1}, p_half::AbstractArray{FT,1}, q, vmr; g₀=9.807) -> Tuple

Computes atmospheric profile fields, including volume mixing ratios (VMR) of H2O, dry and wet volume column densities (VCDs), and layer thicknesses (Δz).

# Arguments
- `T::AbstractArray{FT,1}`: Temperature profile in Kelvin (K).
- `p_half::AbstractArray{FT,1}`: Pressure at half-levels in hectopascals (hPa).
- `q`: Specific humidity as a mass fraction (kg/kg).
- `vmr`: Dictionary containing volume mixing ratios of various trace gases.
- `g₀=9.807`: Gravitational acceleration (m/s²), default is 9.807 m/s².

# Returns
- `p_full`: Pressure at full levels (hPa).
- `p_half`: Pressure at half levels (hPa), same as input.
- `vmr_h2o`: H2O/dry-air molar ratio (unitless).
- `vcd_dry`: Dry volume column density (molec/cm²).
- `vcd_h2o`: Wet volume column density (molec/cm²).
- `new_vmr`: Interpolated volume mixing ratios of trace gases (Dictionary).
- `Δz`: Layer thicknesses (m).

# Description
This function calculates various atmospheric profile fields given temperature, pressure, specific humidity, and initial volume mixing ratios of trace gases. It computes:
1. Pressure at full levels.
2. H2O/dry-air molar ratio from specific humidity.
3. Dry and wet volume column densities (VCDs).
4. Layer thicknesses (Δz).
5. Interpolated volume mixing ratios for other trace gases.
"""
function compute_atmos_profile_fields(T::AbstractArray{FT,1}, p_half::AbstractArray{FT,1}, q, vmr; g₀=9.8032465) where FT
    # g₀ default aligned with sanghavi (Bodhaine 1999 Eq.30 implementation
    # uses gravitational acceleration ≈ 9.8032 m/s² for column-air mass — see
    # sanghavi atmo_prof.jl). q is treated directly as specific humidity
    # (water-vapor mass fraction in kg/kg), matching the sanghavi convention.
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

    # Also get the H2O/dry-air molar-ratio vector.
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )
    Δz      = zeros(FT, n_layers)
    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        # `vmr_h2o` is the dry-air molar ratio r = N_H2O/N_dry. Convert to
        # moist-air mole fractions before computing the mixture molar mass and
        # splitting the hydrostatic total column. Treating r itself as the
        # moist-air mole fraction biases H2O continua at percent-level wetness.
        vmr_h2o[i] = q[i]/(1-q[i]) * ratio
        x_dry = 1 / (1 + vmr_h2o[i])
        x_h2o = vmr_h2o[i] * x_dry
        M  = x_dry * dry_mass + x_h2o * wet_mass
        vcd = Nₐ * Δp / (M  * g₀ * 100^2) * 100
        vcd_dry[i] = x_dry * vcd   # includes m2->cm2
        vcd_h2o[i] = x_h2o * vcd
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
                interp_linear = linear_interpolation(pressure_grid, vmr[molec_i])
                new_vmr[molec_i] = [interp_linear(x) for x in p_full]
            end
        else
            new_vmr[molec_i] = vmr[molec_i]
        end
    end

    return p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz

end

## IO reader methods moved to src/IO/AtmosProfile.jl to decouple CoreRT from IO.

"Geometric half-level altitudes in km above the model BOA (TOA to BOA)."
function half_level_altitudes(profile::AtmosphericProfile{FT}) where {FT}
    n = length(profile.T)
    z_half = zeros(FT, n + 1)
    @inbounds for i in n:-1:1
        z_half[i] = z_half[i + 1] + profile.Δz[i] / FT(1000)
    end
    return z_half
end

"Normalize a vertical input supplied at either layer centers or interfaces."
function _layer_centered_input(name::AbstractString,
                               values::AbstractVector,
                               p_half::AbstractVector{FT}) where {FT}
    n_layers = length(p_half) - 1
    if length(values) == n_layers
        return FT.(values)
    elseif length(values) == length(p_half)
        p_full = (p_half[1:end-1] .+ p_half[2:end]) ./ FT(2)
        interpolation = linear_interpolation(log.(p_half), FT.(values))
        return FT.(interpolation.(log.(p_full)))
    end
    throw(ArgumentError(
        "$name has $(length(values)) entries; expected $n_layers layer values " *
        "or $(length(p_half)) interface values"))
end

"Normalize trace-gas vectors that use the atmospheric interface grid."
function _layer_centered_vmr(vmr, p_half::AbstractVector{FT}) where {FT}
    n_layers = length(p_half) - 1
    normalized = Dict{String, Union{Real, Vector}}()
    for (name, values) in vmr
        normalized[string(name)] = if values isa AbstractVector &&
                                      length(values) in (n_layers, length(p_half))
            _layer_centered_input("VMR $name", values, p_half)
        elseif values isa AbstractVector
            FT.(values)  # retain the existing arbitrary-grid fallback downstream
        else
            FT(values)
        end
    end
    return normalized
end

"Interpolate an intensive layer field to a new full-level pressure grid."
function _interp_layer_field(old_p::AbstractVector{FT}, data::AbstractVector,
                             new_p::AbstractVector{FT}) where {FT}
    length(data) == length(old_p) || throw(ArgumentError(
        "vertical field has $(length(data)) entries but the profile has $(length(old_p)) layers"))
    itp = linear_interpolation(log.(old_p), FT.(data); extrapolation_bc=Flat())
    return FT.(itp.(log.(new_p)))
end

"Reframe all layer quantities onto `p_half_new`, retaining supplied geometry."
function reframe_profile(profile::AtmosphericProfile{FT},
                         p_half_new::AbstractVector{FT};
                         z_half_new::Union{Nothing, AbstractVector{FT}}=nothing) where {FT}
    p_new = collect(p_half_new)
    all(isfinite, p_new) || throw(ArgumentError("pressure interfaces must be finite"))
    all(>(zero(FT)), p_new) || throw(ArgumentError("pressure interfaces must be positive"))
    issorted(p_new) || throw(ArgumentError("pressure interfaces must be ordered TOA to BOA"))
    all(diff(p_new) .> zero(FT)) || throw(ArgumentError("pressure interfaces must be unique"))

    p_full_new = (p_new[1:end-1] .+ p_new[2:end]) ./ FT(2)
    T_new = _interp_layer_field(profile.p_full, profile.T, p_full_new)
    q_new = _interp_layer_field(profile.p_full, profile.q, p_full_new)

    vmr_new = Dict{String, Union{Real, Vector}}()
    for (name, values) in profile.vmr
        vmr_new[name] = values isa AbstractVector ?
            _interp_layer_field(profile.p_full, values, p_full_new) : values
    end

    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, vmr_out, Δz =
        compute_atmos_profile_fields(T_new, p_new, q_new, vmr_new)

    if z_half_new !== nothing
        length(z_half_new) == length(p_new) || throw(ArgumentError(
            "z_half_new must have one altitude per pressure interface"))
        Δz = FT(1000) .* (z_half_new[1:end-1] .- z_half_new[2:end])
        all(Δz .> zero(FT)) || throw(ArgumentError(
            "geometric interfaces must be strictly ordered TOA to BOA"))
    end

    return AtmosphericProfile(T_new, p_full, q_new, p_half, vmr_h2o,
                              vcd_dry, vcd_h2o, vmr_out, Δz)
end

"""
    reframe_vertical_sources(sources, old_p_full, new_p_full)

Interpolate every vertically resolved source input from the user profile's
full-level pressure grid to the final model grid. Sources without a vertical
layer field are returned unchanged. A `ThermalEmission.B_layer` may be supplied
on either the input grid or the already-final grid; input-grid data are
interpolated column-by-column in log pressure using the same flat-extrapolated
scheme as temperature, humidity, and VMR profiles. If the two grids have the
same number of layers, the input-grid interpretation takes precedence because
array shape alone cannot distinguish the two coordinates.
"""
reframe_vertical_sources(sources::AbstractSource, _old_p_full, _new_p_full) = sources

function reframe_vertical_sources(source::ThermalEmission,
                                  old_p_full::AbstractVector{FT},
                                  new_p_full::AbstractVector{FT}) where {FT}
    B = source.B_layer
    B === nothing && return source

    n_input = length(old_p_full)
    n_final = length(new_p_full)
    n_rows = size(B, 1)
    if n_rows == n_input
        B_new = Matrix{FT}(undef, n_final, size(B, 2))
        for iν in axes(B, 2)
            # Source inputs may already reside on a device; profile framing is
            # a host-side construction step, so materialise this small column.
            values = collect(view(B, :, iν))
            B_new[:, iν] .= _interp_layer_field(old_p_full, values, new_p_full)
        end
        return ThermalEmission(B_layer=B_new)
    elseif n_rows == n_final
        return source
    end

    throw(ArgumentError(
        "ThermalEmission.B_layer has $n_rows layers; expected either " *
        "$n_input input-profile layers or $n_final final model layers"))
end

function reframe_vertical_sources(sources::SourceSet, old_p_full, new_p_full)
    return SourceSet(map(source ->
        reframe_vertical_sources(source, old_p_full, new_p_full), sources.sources))
end

"Copy reframed vertical source fields into an existing model source tree."
_copy_vertical_sources!(::AbstractSource, ::AbstractSource) = nothing

function _copy_vertical_sources!(destination::ThermalEmission,
                                 source::ThermalEmission)
    source.B_layer === nothing && return nothing
    destination.B_layer === nothing && throw(ArgumentError(
        "cannot add a per-layer ThermalEmission to an existing zero-placeholder source"))
    size(destination.B_layer) == size(source.B_layer) || throw(ArgumentError(
        "updated ThermalEmission grid has size $(size(source.B_layer)); " *
        "BatchContext was built with size $(size(destination.B_layer))"))
    destination.B_layer .= source.B_layer
    return nothing
end

function _copy_vertical_sources!(destination::SourceSet, source::SourceSet)
    length(destination) == length(source) || throw(ArgumentError(
        "updated source composition differs from the BatchContext source composition"))
    for (dest, src) in zip(destination.sources, source.sources)
        _copy_vertical_sources!(dest, src)
    end
    return nothing
end

"Allocate exactly `n` reduced layers among intervals bounded by required anchors."
function _anchored_pressure_grid(profile::AtmosphericProfile{FT}, n::Int,
                                 required_p::AbstractVector{FT}) where {FT}
    anchors = sort!(unique!(vcat(FT[profile.p_half[1]], collect(required_p),
                                 FT[profile.p_half[end]])))
    n_segments = length(anchors) - 1
    n >= n_segments || throw(ArgumentError(
        "$n layers cannot contain $(length(required_p)) required interior interfaces"))

    allocation = ones(Int, n_segments)
    remaining = n - n_segments
    if remaining > 0
        spans = diff(anchors)
        raw = remaining .* spans ./ sum(spans)
        extra = floor.(Int, raw)
        allocation .+= extra
        left = remaining - sum(extra)
        if left > 0
            order = sortperm(raw .- extra; rev=true)
            allocation[order[1:left]] .+= 1
        end
    end

    p_half = FT[]
    for i in eachindex(allocation)
        segment = collect(range(anchors[i], anchors[i + 1], length=allocation[i] + 1))
        append!(p_half, i == firstindex(allocation) ? segment : segment[2:end])
    end
    return p_half
end

"Resolve scalar/vector observer-height semantics against a concrete atmosphere."
function _resolve_observer_request(obs_alt::Union{FT, Vector{FT}},
                                   toa_altitude::FT) where {FT}
    raw = obs_alt isa Real ? FT[obs_alt] : copy(obs_alt)
    isempty(raw) && throw(ArgumentError("obs_alt vector must not be empty"))
    all(isfinite, raw) || throw(ArgumentError("obs_alt values must be finite"))
    all(>=(zero(FT)), raw) || throw(ArgumentError(
        "obs_alt values must be nonnegative heights in km above BOA"))

    # A height obtained from the model and then round-tripped through user code
    # can differ from the resolved TOA by a few ulps. Treat only
    # machine-equivalent values as the endpoint; materially lower heights remain
    # strict-interior requests.
    near_toa(h) = isapprox(h, toa_altitude;
                           atol=zero(FT), rtol=FT(4) * eps(FT))
    is_vector = obs_alt isa Vector
    has_zero = any(iszero, raw)
    include_boa = is_vector ? has_zero : iszero(only(raw))
    include_toa = is_vector ?
                  (has_zero || any(h -> h > toa_altitude || near_toa(h), raw)) :
                  (!iszero(only(raw)) &&
                   (only(raw) > toa_altitude || near_toa(only(raw))))
    interior = sort!(unique(filter(
        h -> zero(FT) < h < toa_altitude && !near_toa(h), raw)); rev=true)

    return (; include_toa, include_boa, interior_altitudes=interior,
            requested_as_vector=is_vector)
end

"Index of an atmospheric interface that is machine-equivalent to `height`."
function _matching_height_index(z_half::AbstractVector{FT}, height::FT) where {FT}
    return findfirst(z -> isapprox(z, height;
                                   atol=zero(FT), rtol=FT(4) * eps(FT)),
                     z_half)
end

"Snap machine-equivalent observer heights to the exact existing interfaces."
function _snap_observer_altitudes(altitudes::AbstractVector{FT},
                                  z_half::AbstractVector{FT}) where {FT}
    snapped = FT[]
    for height in altitudes
        index = _matching_height_index(z_half, height)
        push!(snapped, index === nothing ? height : z_half[index])
    end
    return sort!(unique!(snapped); rev=true)
end

"Pressure at each requested interior geometric altitude."
function _pressures_at_altitudes(profile::AtmosphericProfile{FT},
                                 z_half::AbstractVector{FT},
                                 altitudes::AbstractVector{FT}) where {FT}
    isempty(altitudes) && return FT[]
    # `z_half` descends while interpolation grids must ascend.
    logp_of_z = linear_interpolation(reverse(z_half), reverse(log.(profile.p_half)))
    pressures = FT[]
    for height in altitudes
        index = _matching_height_index(z_half, height)
        push!(pressures, index === nothing ?
              exp(logp_of_z(height)) : profile.p_half[index])
    end
    return pressures
end

"Map exact interior altitudes to interface indices after reframing."
function _sensor_interface_indices(profile::AtmosphericProfile{FT},
                                   altitudes::AbstractVector{FT}) where {FT}
    isempty(altitudes) && return Int[]
    z_half = half_level_altitudes(profile)
    levels = Int[]
    for h in altitudes
        # Requested heights are written exactly into the reframed geometric
        # grid. Keep the fallback restricted to interior interfaces so a very
        # high but valid sensor cannot be mistaken for TOA.
        j = findfirst(==(h), z_half)
        if j === nothing
            interior = @view z_half[2:end-1]
            local_index = _matching_height_index(interior, h)
            j = local_index === nothing ? nothing : local_index + 1
        end
        j === nothing && throw(ArgumentError(
            "observer altitude $h km was not retained as an atmospheric interface"))
        1 < j < length(z_half) || throw(ArgumentError(
            "observer altitude $h km resolved to an endpoint, not an interior interface"))
        push!(levels, j - 1)
    end
    return levels
end

"""
    prepare_observer_profile(T, p_half, q, vmr, obs_alt, profile_reduction_n)

Build the final atmospheric grid and observer metadata. Heights are geometric
kilometres above the model BOA. Every strict-interior requested height is an
exact half-level interface; all intensive inputs are interpolated to the new
full levels and all integrated fields are recomputed.
"""
function prepare_observer_profile(T::Vector{FT}, p_half::Vector{FT}, q::Vector{FT},
                                  vmr, obs_alt::Union{FT, Vector{FT}},
                                  profile_reduction_n::Int) where {FT}
    (profile_reduction_n == -1 || profile_reduction_n > 0) ||
        throw(ArgumentError("profile_reduction must be -1 or a positive integer"))
    length(p_half) == length(T) + 1 || throw(ArgumentError(
        "pressure profile must have one more interface than the temperature profile has layers"))
    all(isfinite, T) || throw(ArgumentError("temperature profile must be finite"))
    all(isfinite, p_half) && all(>(zero(FT)), p_half) &&
        all(diff(p_half) .> zero(FT)) || throw(ArgumentError(
            "pressure interfaces must be finite, positive, and ordered TOA to BOA"))
    q_layer = _layer_centered_input("specific humidity", q, p_half)
    vmr_layer = _layer_centered_vmr(vmr, p_half)
    all(x -> isfinite(x) && zero(FT) <= x < one(FT), q_layer) ||
        throw(ArgumentError(
            "specific humidity must contain finite mass fractions in [0, 1)"))
    for (name, values) in vmr_layer
        samples = values isa AbstractVector ? values : (values,)
        all(isfinite, samples) || throw(ArgumentError(
            "VMR $name must contain only finite values"))
    end
    p_full, p_out, vmr_h2o, vcd_dry, vcd_h2o, vmr_out, Δz =
        compute_atmos_profile_fields(T, p_half, q_layer, vmr_layer)
    base = AtmosphericProfile(T, p_full, q_layer, p_out, vmr_h2o,
                              vcd_dry, vcd_h2o, vmr_out, Δz)
    z_base = half_level_altitudes(base)
    toa_altitude = first(z_base)
    selection = _resolve_observer_request(obs_alt, toa_altitude)
    snapped_altitudes = _snap_observer_altitudes(
        selection.interior_altitudes, z_base)
    selection = merge(selection, (; interior_altitudes=snapped_altitudes))
    requested_p = _pressures_at_altitudes(base, z_base, selection.interior_altitudes)

    # First build the finest grid implied by the input profile and requested
    # observer interfaces. A profile-reduction target at least this large is a
    # no-op: "reduction" must never expand an already adequate grid.
    natural_profile = if isempty(requested_p)
        base
    else
        p_new = sort!(unique!(vcat(copy(base.p_half), requested_p)))
        z_of_logp = linear_interpolation(log.(base.p_half), z_base)
        z_new = FT.(z_of_logp.(log.(p_new)))
        for (h, p_h) in zip(selection.interior_altitudes, requested_p)
            z_new[argmin(abs.(p_new .- p_h))] = h
        end
        reframe_profile(base, p_new; z_half_new=z_new)
    end

    final_profile = if profile_reduction_n == -1
        natural_profile
    else
        minimum_layers = length(requested_p) + 1
        effective_n = max(profile_reduction_n, minimum_layers)
        if effective_n != profile_reduction_n
            @warn "profile_reduction cannot preserve all requested observer heights; increasing layer count" requested_layers=profile_reduction_n requested_heights_km=selection.interior_altitudes effective_layers=effective_n
        end
        effective_n >= length(natural_profile.T) ? natural_profile : begin
            p_new = _anchored_pressure_grid(base, effective_n, requested_p)
            z_of_logp = linear_interpolation(log.(base.p_half), z_base)
            z_new = FT.(z_of_logp.(log.(p_new)))
            for (h, p_h) in zip(selection.interior_altitudes, requested_p)
                z_new[argmin(abs.(p_new .- p_h))] = h
            end
            reframe_profile(base, p_new; z_half_new=z_new)
        end
    end

    sensor_levels = _sensor_interface_indices(final_profile, selection.interior_altitudes)
    return final_profile, (; selection..., sensor_levels, toa_altitude)
end

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
        itp  = linear_interpolation(grid, data)
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
        x_dry   = inv(FT(1) + vmr_h2o[i])
        x_h2o   = vmr_h2o[i] * x_dry
        M       = x_dry * dry_mass + x_h2o * wet_mass
        vcd     = Nₐ * Δp / (M * g₀ * FT(100)^2) * FT(100)
        vcd_dry[i] = x_dry * vcd
        vcd_h2o[i] = x_h2o * vcd
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

"CDF of an altitude-form lognormal, including the z=0 lower boundary."
@inline function _altitude_lognormal_cdf(z, dist::LogNormal)
    z <= zero(z) && return zero(z + dist.μ + dist.σ)
    a = (log(z) - dist.μ) / dist.σ
    return erfc(-a / sqrt(oftype(a, 2))) / oftype(a, 2)
end

"""
Allocate an altitude-form lognormal aerosol by exact geometric-layer integrals.

`LogNormal(log(z₀), σ₀)` is defined in altitude `z` [km].  The returned layer
fractions are CDF differences between geometric interfaces, normalized over
the finite model column so their sum is exactly `total_τ`.
"""
function getAerosolLayerOptProp(total_τ::FT, dist::LogNormal,
                                profile::AtmosphericProfile) where FT
    z_half = half_level_altitudes(profile)
    ρ = [_altitude_lognormal_cdf(z_half[i], dist) -
         _altitude_lognormal_cdf(z_half[i + 1], dist)
         for i in eachindex(profile.Δz)]
    norm_ρ = sum(ρ)
    iszero(norm_ρ) && throw(ArgumentError(
        "altitude-lognormal aerosol has zero mass inside the model column"))
    return total_τ .* ρ ./ norm_ρ
end

@inline _h2o_moist_mole_fraction(r::Real) = r / (one(r) + r)

# Legacy vSmartMOM interpolation tables have their broadener abundance fixed
# when the table is built, so their ordinary cross-section call remains the
# fallback. AtmosphericAbsorption's line-by-line model supports a per-call
# moist-air mole fraction and must receive it for H₂O self broadening.
_layer_absorption_cross_section(model, grid, p, T, ::Nothing) =
    absorption_cross_section(model, grid, p, T)

_layer_absorption_cross_section(
    model::AtmosphericAbsorption.LineByLineModel, grid, p, T,
    self_broadener_vmr::Real) =
    absorption_cross_section(model, grid, p, T; vmr=self_broadener_vmr)

_layer_absorption_cross_section(model, grid, p, T, ::Real) =
    absorption_cross_section(model, grid, p, T)

# Native ABSCO tables carry a tabulated H₂O-broadener axis — route the
# per-layer broadener VMR into the LUT lookup (surface-split behavior).
_layer_absorption_cross_section(model::AtmosphericAbsorption.AbscoLUT, grid, p, T,
                                h2o_vmr::Real) =
    AtmosphericAbsorption.compute_cross_section(model, grid, p, T;
                                                vmr=h2o_vmr, interp=:linear)

# When the caller supplies no explicit broadener, only models with a
# tabulated H₂O-broadener axis (native ABSCO) draw it from the profile;
# legacy interpolation tables and line-by-line models keep their broadener
# semantics untouched (LBL self-broadening is an explicit, moist-air mole
# fraction supplied by `compute_h2o_absorption_profile!`).
_default_broadener(model, profile, iz) = nothing
_default_broadener(::AtmosphericAbsorption.AbscoLUT, profile, iz) =
    profile.vmr_h2o[iz]

"Given the CrossSectionModel, the grid, and the AtmosphericProfile, fill up the τ_abs array with the cross section at each layer
(using pressures/temperatures) from the profile. `self_broadener_vmr`, when
provided, is a moist-air mole fraction passed to AtmosphericAbsorption's
line-by-line cross-section call; when it is `nothing`, the per-layer
`profile.vmr_h2o` supplies the broadener so native ABSCO tables (which carry
a tabulated H₂O-broadener axis — dispatch below) keep their surface-split
behavior."
function compute_absorption_profile!(τ_abs::Array{FT,2}, 
                                     absorption_model, 
                                     grid,
                                     vmr,
                                     profile::AtmosphericProfile,
                                     ;
                                     self_broadener_vmr=nothing,
                                     batched::Bool=true,
                                     ) where FT 

    # The array to store the cross-sections must be same length as number of layers
    @assert size(τ_abs,2) == length(profile.p_full)
    @assert length(vmr) ==1 || length(vmr) == length(profile.p_full)  "Length of VMR array has to match profile size or be uniform"

    # BATCHED fast path for the tabulated models: the whole profile in ONE
    # kernel via `compute_cross_section_profile`, instead of a host loop whose
    # per-layer query cost ~7 broadcast launches, a grid upload, a full device
    # synchronize and a device->host copy — ~10 launches and 3 sync/transfer
    # points per layer, x72 layers x2 builds per cell-band-window. Per-element
    # arithmetic (bracket/clamp semantics, blend order, zero outside the nu
    # range, and the vcd*vmr weighting below) is IDENTICAL to the loop, so CPU
    # results are bit-exact; see AtmosphericAbsorption/test_profile_batch.jl.
    # The point axis is flat, so this same call can batch MANY cells' layers
    # at once when a caller concatenates them.
    if batched && _supports_batched_profile(absorption_model)
        n = length(profile.p_full)
        broadeners = _batched_broadener(absorption_model, profile,
                                        self_broadener_vmr, n)
        σmat = Array(AtmosphericAbsorption.compute_cross_section_profile(
            absorption_model, grid, profile.p_full, profile.T;
            vmr = broadeners, interp = :linear))
        @inbounds for iz in 1:n
            vmr_curr = vmr isa AbstractArray ? vmr[iz] : vmr
            @views τ_abs[:, iz] .+= σmat[:, iz] .* (profile.vcd_dry[iz] * vmr_curr)
        end
        return
    end

    @showprogress 1 for iz in eachindex(profile.p_full)

        # Pa -> hPa
        p = profile.p_full[iz]
        T = profile.T[iz]
        # Either use the current layer's vmr, or use the uniform vmr
        vmr_curr = vmr isa AbstractArray ? vmr[iz] : vmr
        broadener_curr = self_broadener_vmr isa AbstractArray ?
            self_broadener_vmr[iz] : self_broadener_vmr
        #@show vmr_curr
        # An explicit kwarg wins (LBL self-broadening); otherwise the model
        # type decides whether the profile supplies a broadener (ABSCO only).
        broadener_curr = broadener_curr === nothing ?
            _default_broadener(absorption_model, profile, iz) : broadener_curr
        σ = collect(_layer_absorption_cross_section(
            absorption_model, grid, p, T, broadener_curr))
        @views τ_abs[:,iz] .+= σ .* (profile.vcd_dry[iz] * vmr_curr)
    end
    
end

# The batched path exists for the tabulated models when the installed
# AtmosphericAbsorption ships `compute_cross_section_profile` (>= the 2026-08
# profile-batch addition); anything else keeps the per-layer loop. LBL models
# stay on the loop deliberately: their per-layer state (self-broadening) is
# already a single fused kernel per layer.
_supports_batched_profile(model) =
    isdefined(AtmosphericAbsorption, :compute_cross_section_profile) &&
    (model isa AtmosphericAbsorption.AbscoLUT ||
     model isa AtmosphericAbsorption.InterpolationModel)

# Per-layer broadener vector for the batched call, mirroring the loop exactly:
# an explicit kwarg wins (vector used as-is, scalar broadcast); otherwise
# ABSCO draws the profile's own H2O and everything else passes nothing (an
# InterpolationModel tabulates a fixed vmr, and its scalar path ignores any
# broadener — dispatch `_layer_absorption_cross_section(model, _, _, _, ::Real)`
# falls through to the broadener-free query).
function _batched_broadener(model, profile, self_broadener_vmr, n)
    model isa AtmosphericAbsorption.InterpolationModel && return nothing
    self_broadener_vmr isa AbstractArray && return self_broadener_vmr
    self_broadener_vmr !== nothing && return fill(self_broadener_vmr, n)
    return profile.vmr_h2o
end

"""
    compute_h2o_absorption_profile!(τ_abs, absorption_model, grid, profile)

Accumulate H₂O line optical depth. `profile.vmr_h2o` stores the dry-air
molar ratio `r = N_H2O/N_dry`, which remains the abundance multiplying the dry
column. Line broadening instead requires the moist-air H₂O mole fraction
`x_H2O = r/(1+r)`.
"""
function compute_h2o_absorption_profile!(τ_abs::Array{FT,2},
                                         absorption_model,
                                         grid,
                                         profile::AtmosphericProfile) where FT
    x_h2o = _h2o_moist_mole_fraction.(profile.vmr_h2o)
    return compute_absorption_profile!(
        τ_abs, absorption_model, grid, profile.vmr_h2o, profile;
        self_broadener_vmr=x_h2o)
end
