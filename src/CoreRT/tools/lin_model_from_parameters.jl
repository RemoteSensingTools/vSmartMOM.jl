"Set a common l_max in a given band for each aerosol type"
function set_uniform_lmax!(lmax::Vector{Int}, aerosol_optics, lin_aerosol_optics=nothing)
    n_bands = length(lmax)
    n_aer = length(aerosol_optics[1])
    for i_band = 1:n_bands
        pref_lmax = lmax[i_band]
        for i_aer = 1:n_aer
            aer_lmax = size(aerosol_optics[i_band][i_aer].greek_coefs.β, 1)
            if aer_lmax < pref_lmax
                for fname in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    old = getfield(aerosol_optics[i_band][i_aer].greek_coefs, fname)
                    new_arr = zeros(eltype(old), pref_lmax, size(old)[2:end]...)
                    view(new_arr, 1:aer_lmax, ntuple(_ -> Colon(), ndims(old)-1)...) .= old
                    setfield!(aerosol_optics[i_band][i_aer].greek_coefs, fname, new_arr)
                end
                if lin_aerosol_optics !== nothing
                    for fname in (:α̇, :β̇, :γ̇, :δ̇, :ϵ̇, :ζ̇)
                        old = getfield(lin_aerosol_optics[i_band][i_aer].lin_greek_coefs, fname)
                        new_arr = zeros(eltype(old), size(old,1), pref_lmax, size(old)[3:end]...)
                        view(new_arr, :, 1:aer_lmax,
                             ntuple(_ -> Colon(), ndims(old)-2)...) .= old
                        setfield!(lin_aerosol_optics[i_band][i_aer].lin_greek_coefs,
                                  fname, new_arr)
                    end
                end
            end
        end
    end
end

@inline function _zero_lin_greek(greek, ::Type{FT}) where {FT<:AbstractFloat}
    make_zero(field) = zeros(FT, 4, size(getfield(greek, field))...)
    return linGreekCoefs(make_zero(:α), make_zero(:β), make_zero(:γ),
                         make_zero(:δ), make_zero(:ϵ), make_zero(:ζ))
end

@inline function _zero_microphysics_field(value, ::Type{FT}) where {FT<:AbstractFloat}
    return value isa Number ? zeros(FT, 4) : zeros(FT, 4, size(value)...)
end

"Create structurally valid zero Mie tangents when microphysics is fixed."
function _zero_aerosol_microphysics_jacobian(aerosol, ::Type{FT}) where {FT<:AbstractFloat}
    phase_lin = aerosol.phase_greek === nothing ? nothing :
        [_zero_lin_greek(greek, FT) for greek in aerosol.phase_greek]
    return linAerosolOptics(
        lin_greek_coefs=_zero_lin_greek(aerosol.greek_coefs, FT),
        ω̃̇=_zero_microphysics_field(aerosol.ω̃, FT),
        k̇=_zero_microphysics_field(aerosol.k, FT),
        ḟᵗ=_zero_microphysics_field(aerosol.fᵗ, FT),
        phase_ν=aerosol.phase_ν === nothing ? nothing : FT.(aerosol.phase_ν),
        phase_lin_greek=phase_lin)
end

const _H2O_SELF_BROADENING = "H₂O self broadening of its line cross sections"

"Species label of each native gas slot in band `i_band` (slot 1 is q-driven H₂O)."
function _native_gas_species(abs_params, i_band, n_slots)
    variable = isnothing(abs_params) ? String[] : abs_params.variable_molecules[i_band]
    return [slot == 1 ? "q-driven H₂O" :
            slot - 1 <= length(variable) ? String(variable[slot - 1]) :
            "unused gas slot" for slot in 1:n_slots]
end

"Record that `term` responds to native gas `slot` but is absent from its tangent."
function _omit_abundance_term!(omitted, i_band, slot, term)
    term in omitted[i_band][slot] || push!(omitted[i_band][slot], term)
    return nothing
end

# ABSCO tables of every gas read the profile's H₂O broadener (line-by-line
# models of other gases take none; see `_default_broadener`). The q-driven H₂O
# columns receive that axis's exact response in both abundance modes; a model
# that reads the broadener without a tangent records the term as omitted.
function _profile_broadener_tangent!(τ̇_abs, omitted, i_band, absorber, mol_name,
                                     grid, vmr, profile)
    (_default_broadener(absorber, profile, 1) === nothing ||
     !_cross_section_uses_broadener(absorber)) && return nothing
    _has_broadener_tangent(absorber) || return _omit_abundance_term!(
        omitted, i_band, 1, "H₂O-broadener response of the $(mol_name) cross section")
    return _accumulate_broadener_tangent!(τ̇_abs, 1, absorber, grid, vmr, profile)
end

# H₂O's own broadener response is in its tangent only for ABSCO tables;
# line-by-line self broadening is the omitted term.
_omits_h2o_broadener(model) =
    _cross_section_uses_broadener(model) && !_has_broadener_tangent(model)

# A CIA pair responds to its collider abundances on this band grid wherever a
# selected panel tabulates a nonzero cross section, even if a collider's
# current VMR (and hence the present opacity) is zero.
_cia_has_cross_section(table) = any(σ -> isfinite(σ) && !iszero(σ), table.σ_νT)

#= Molecular-column parts of dτ_abs/dp_surf.

At fixed layer T, q, and mixing ratio, a line absorber has
τᵢ = σ(p_full,i, Tᵢ) · N_dry,i · vmrᵢ, so dτᵢ/dp_surf has a cross-section part
(accumulated by `_accumulate_absorption_pressure!`) and the column part
τᵢ · (dN_dry,i/dp_surf)/N_dry,i added here. Binary absorbers (CIA, MT_CKD)
scale with a number density as well as a column, τᵢ ∝ p_full,i · N_i, and
receive the extra factor (dp_full,i/dp_surf)/p_full,i. =#

"Add the column part of the line-absorption surface-pressure tangent."
function _add_column_psurf_tangent!(τ̇_psurf, τ_abs, ::FixedPressureGrid,
                                    tangents, profile)
    # Only the bottom layer's column depends on p_surf.
    dry_ratio_dot = tangents.vcd_dry_dot[end] / profile.vcd_dry[end]
    τ̇_psurf[:, end] .+= τ_abs[:, end] .* dry_ratio_dot
    return τ̇_psurf
end

function _add_column_psurf_tangent!(τ̇_psurf, τ_abs, ::FixedAltitudeGrid,
                                    tangents, profile)
    # Every layer column scales with p_surf: dN/dp_surf = N/p_surf.
    dry_ratio_dot = tangents.vcd_dry_dot ./ profile.vcd_dry
    τ̇_psurf .+= τ_abs .* dry_ratio_dot'
    return τ̇_psurf
end

"Add the surface-pressure tangent of a binary (density × column) absorber."
function _add_binary_psurf_tangent!(τ̇_psurf, τ_binary, ::FixedPressureGrid,
                                    tangents, profile)
    FT = eltype(profile.p_full)
    dry_ratio_dot = tangents.vcd_dry_dot[end] / profile.vcd_dry[end]
    binary_ratio_dot = dry_ratio_dot + FT(0.5) / profile.p_full[end]
    τ̇_psurf[:, end] .+= @view(τ_binary[:, end]) .* binary_ratio_dot
    return τ̇_psurf
end

function _add_binary_psurf_tangent!(τ̇_psurf, τ_binary, ::FixedAltitudeGrid,
                                    tangents, profile)
    binary_ratio_dot = tangents.vcd_dry_dot ./ profile.vcd_dry .+
                       tangents.p_full_dot ./ profile.p_full
    τ̇_psurf .+= τ_binary .* binary_ratio_dot'
    return τ̇_psurf
end

"""
    model_from_parameters(::LinMode, params::vSmartMOM_Parameters;
                          sources = SolarBeam(), external_solar = false,
                          aerosol_anchor_bands = nothing,
                          compute_aerosol_microphysics_jacobians = true,
                          compute_h2o_jacobians = true,
                          compute_pressure_jacobians = true,
                          abundance_jacobians = :complete,
                          vertical_grid = FixedPressureGrid())

Construct both the forward `RTModel` and the linearized `RTModelLin` objects
from the input parameters, for use in linearized (Jacobian) RT computations.

This is the **linearized** counterpart of `model_from_parameters(params)`. It computes:

1. **Forward optical properties** (same as the forward-only model):
   - Rayleigh optical depth per layer and spectral point.
   - Aerosol optical properties (via Mie theory) with δ-M truncation.
   - Trace gas absorption cross-sections and optical depths.

2. **Linearized optical properties** (derivatives of the above):
   - `τ̇_aer[iB][iaer, 7, nSpec, nLayers]`: Derivatives of aerosol τ w.r.t. 7 sub-parameters
     `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]` per aerosol type.
   - `τ̇_abs[iB][NGas*Nz, nSpec, Nz]`: Layer-resolved derivatives of gas
     absorption τ w.r.t. `VMR(gas, z)`, flattened species-major.
   - `lin_aerosol_optics[iB][iaer]`: Derivatives of Mie properties (ω̃, fᵗ, greek coefficients)
     w.r.t. Mie parameters `[nᵣ, nᵢ, rₘ, σᵣ]`.

# Returns
- `model::RTModel`: Forward model (optical properties, geometry, quadrature).
- `lin_model::RTModelLin`: Linearized model (all derivative arrays).

# Notes
- Strict-interior observer heights are inserted as exact atmospheric
  interfaces during model construction and are supported by the elastic
  analytic-linearized solver.
- Aerosol Mie derivatives of the extinction cross-section, single-scattering
  albedo, truncation factor and greek coefficients with respect to
  `[nᵣ, nᵢ, rₘ, σᵣ]` come from the analytic linearized NAI-2 calculation
  (`compute_aerosol_optical_properties(LinMode(), mie_model, FT)`), which runs
  on the CPU.
- Set `compute_aerosol_microphysics_jacobians=false` when all four Mie
  parameters are fixed. The forward aerosol optics are retained and
  structurally valid zero tangents are supplied to the downstream mixer.
- Set `compute_h2o_jacobians=false` when the q-driven H2O profile is fixed.
  Its forward absorption is retained without generating H2O tangent entries.
- Without an `absorption` block the model has no gas absorption, as in the
  forward constructor (`q` still sets the dry-air column).
- Native gas columns are the partial `∂τ/∂VMR = σ·N_dry`: line cross
  sections and the CIA and MT_CKD optical depths are held fixed. The one
  exception is ABSCO tables with an H₂O-broadener axis: the q-driven H₂O
  columns include every such table's response to the H₂O abundance (exact
  for the linear interpolation in the tabulated VMR). Where the forward
  absorption responds to an abundance in another way (line-by-line H₂O self
  broadening, MT_CKD, CIA colliders), the default
  `abundance_jacobians=:complete` rejects solves or plans requesting those
  columns. Pass `abundance_jacobians=:frozen_cross_section` to accept these
  omissions explicitly; tangent values are identical in both settings.
- Set `compute_pressure_jacobians=false` for a retrieval with fixed pressure.
  Pressure tangent fields are then `nothing`, not zero sensitivities; use a
  `PlannedRTModelLin` whose layout omits pressure. Full-layout solves reject
  these models. Forward-only absorber extensions remain usable in such plans.
- `vertical_grid` selects how the vertical grid responds to surface pressure
  in the surface-pressure tangents: [`FixedPressureGrid`](@ref) (default;
  only the bottom boundary moves) or [`FixedAltitudeGrid`](@ref) (every
  interface fixed in geometric altitude, so every interface pressure scales
  with `p_surf`). The forward optics do not depend on it. A retrieval that
  varies `p_surf` must rebuild `params.p` with the same convention, via
  [`surface_pressure_interfaces`](@ref), for the Jacobian to be the derivative
  of its forward model.
- `aerosol_anchor_bands` has the same meaning as in the forward constructor:
  it fixes aerosol spectral endpoint calculations to canonical full-band
  grids while evaluating the interpolated optics on chunked or
  shoulder-expanded solve grids. Forward and tangent optics always share the
  same anchors.
"""
function model_from_parameters(lin::LinMode,
    params::vSmartMOM_Parameters;
    sources::AbstractSource = SolarBeam(),
    external_solar::Bool = false,
    aerosol_anchor_bands=nothing,
    compute_aerosol_microphysics_jacobians::Bool = true,
    compute_h2o_jacobians::Bool = true,
    compute_pressure_jacobians::Bool = true,
    abundance_jacobians::Symbol = :complete,
    vertical_grid::AbstractVerticalGridConvention = FixedPressureGrid())
    abundance_jacobians in (:complete, :frozen_cross_section) || throw(ArgumentError(
        "abundance_jacobians must be :complete or :frozen_cross_section; " *
        "got :$abundance_jacobians"))
    FT = params.float_type
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
    aerosol_bands = aerosol_anchor_bands === nothing ? params.spec_bands :
                     aerosol_anchor_bands
    length(aerosol_bands) == n_bands || throw(DimensionMismatch(
        "aerosol_anchor_bands must have one grid per spectral band"))
    all(!isempty, aerosol_bands) || throw(ArgumentError(
        "aerosol interpolation anchor grids must be nonempty"))
    _require_surface_per_band(params)
    scat = params.scattering_params
    abs_params = params.absorption_params
    if scat !== nothing && any(_has_analytic_phase_function, scat.rt_aerosols)
        throw(ArgumentError(
            "model_from_parameters(LinMode(), ...) currently supports Mie aerosols only; " *
            "analytic phase-function aerosols do not yet define Mie-parameter Jacobians."))
    end

    truncation_type = _resolved_truncation(params, params.float_type)

    vmr = isnothing(abs_params) ? Dict() : abs_params.vmr
    T_ft, p_ft, q_ft = convert(Vector{FT}, params.T), convert(Vector{FT}, params.p), convert(Vector{FT}, params.q)
    obs_alt = params.obs_alt isa Real ? FT(params.obs_alt) : convert(Vector{FT}, params.obs_alt)
    profile, observation = prepare_observer_profile(
        T_ft, p_ft, q_ft, vmr, obs_alt, params.profile_reduction_n)
    input_p_full = (p_ft[1:end-1] .+ p_ft[2:end]) ./ FT(2)
    sources = reframe_vertical_sources(sources, input_p_full, profile.p_full)
    obs_geom = ObsGeometry{FT}(
        FT(params.sza), convert(Vector{FT}, params.vza), convert(Vector{FT}, params.vaz), obs_alt,
        observation.sensor_levels, observation.interior_altitudes,
        observation.include_toa, observation.include_boa, observation.toa_altitude)
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom,
                                 params.polarization_type, array_type(params.architecture);
                                 external_solar)
    rayleigh_molecular_T = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry)

    greek_cabannes = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    greek_rayleigh = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    ϖ_Cabannes = zeros(FT, n_bands)
    τ_rayl = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands]

    # Gas-tangent storage follows a ForwardDiff Dual CO2 VMR (hybrid AD) and is
    # otherwise the model float type.
    vmr_co2_type = isnothing(abs_params) || !haskey(abs_params.vmr, "CO2") ?
                   params.float_type : eltype(abs_params.vmr["CO2"])
    FT2 = vmr_co2_type <: ForwardDiff.Dual ? vmr_co2_type : params.float_type
    τ_abs     = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    N_fix_gas = isnothing(abs_params) ? 0 : length(unique(Iterators.flatten(abs_params.fixed_molecules)))
    N_var_gas = isnothing(abs_params) ? 0 : length(unique(Iterators.flatten(abs_params.variable_molecules)))
    # Gas Jacobian columns are species-major and layer-resolved:
    # (gas 1, z=1:Nz), (gas 2, z=1:Nz), ... .  At layer iz only the
    # corresponding gas/layer column is nonzero at the optics boundary.
    N_gas_species = 1 + N_var_gas
    Nz = length(profile.p_full)
    τ̇_abs     = [zeros(FT2, N_gas_species * Nz,
                        length(params.spec_bands[i]), Nz) for i in 1:n_bands]
    τ̇_rayl_psurf = compute_pressure_jacobians ?
        [zeros(FT2, length(params.spec_bands[i]), Nz) for i in 1:n_bands] : nothing
    τ̇_abs_psurf = compute_pressure_jacobians ?
        [zeros(FT2, length(params.spec_bands[i]), Nz) for i in 1:n_bands] : nothing
    τ̇_aer_psurf = compute_pressure_jacobians ?
        [zeros(FT2, n_aer, length(params.spec_bands[i]), Nz) for i in 1:n_bands] : nothing
    psurf_tangents = compute_pressure_jacobians ?
        psurf_profile_tangents(vertical_grid, profile) : nothing
    # Forward terms that respond to a native gas abundance but are absent from
    # its tangent, per band and gas slot (see `AbundanceTangentContract`).
    h2o_tangent = compute_h2o_jacobians && !isnothing(params.q) && any(!iszero, params.q)
    gas_species = [_native_gas_species(abs_params, i_band, N_gas_species)
                   for i_band in 1:n_bands]
    omitted = [[String[] for _ in 1:N_gas_species] for _ in 1:n_bands]
    l_max = zeros(Int, n_bands)
    l_max_aer = zeros(Int, n_aer, n_bands)

    for i_band=1:n_bands
        pressure_tangent = compute_pressure_jacobians ? τ̇_abs_psurf[i_band] : nothing

        # `params` may have been parsed as Float64 and subsequently switched
        # to Float32 by a caller. Keep the wavelength grid consistent with the
        # reframed profile and Rayleigh arrays, as the forward constructor does.
        curr_band_λ = FT.(FT(1e4) ./ params.spec_bands[i_band])
        νₘ = FT(0.5)*(params.spec_bands[i_band][1]+params.spec_bands[i_band][end])
        λₘ = FT(1.e7)/νₘ
        # Per-band molecular-constant depolarizations. ϖ_Cabannes is always
        # taken from the molecular path; the depol values feed greek coefs and
        # τ_rayl only when params.depol < 0 (auto). See model_from_parameters.jl
        # for the full rule.
        ϖ_Cabannes[i_band], depol_air_Cabannes, depol_air_Rayleigh =
            _auto_molecular_depolarizations(params.rayleigh_anisotropy, λₘ,
                                            FT(rayleigh_molecular_T))

        depol_use_Cab = params.depol < 0 ? FT(depol_air_Cabannes) : FT(params.depol)
        depol_use_Ray = params.depol < 0 ? FT(depol_air_Rayleigh) : FT(params.depol)

        push!(greek_cabannes, Scattering.get_greek_rayleigh(depol_use_Cab))
        push!(greek_rayleigh, Scattering.get_greek_rayleigh(depol_use_Ray))

        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end],
                                curr_band_λ,
                                depol_use_Ray, profile.vcd_dry)
        # τ_rayl = τ_total(p_surf) * vcd_dry/sum(vcd_dry).
        if compute_pressure_jacobians
            dry = profile.vcd_dry
            drydot = psurf_tangents.vcd_dry_dot
            total = vec(sum(τ_rayl[i_band], dims=2))
            totaldot = total ./ profile.p_half[end]
            # Use the factored quotient rule: the direct form squares the total
            # dry column (~1e25), which overflows in Float32 and makes the entire
            # surface-pressure Jacobian NaN.
            frac, fracdot = _normalized_column_fraction_tangent(dry, drydot)
            τ̇_rayl_psurf[i_band] .= totaldot * frac' .+ total * fracdot'
        end

        # No absorption block means no gas absorption (H₂O included), as in
        # the forward constructor; q still enters the dry-air column.
        isnothing(abs_params) && continue

        if !isnothing(params.q) && any(!iszero, params.q)
            jac_idx = 1
            h2o_setting = isempty(abs_params.h2o_lut) ? nothing : abs_params.h2o_lut[i_band]
            if h2o_setting === :disabled
                nothing
            elseif h2o_setting !== nothing
                h2o_tangent && _omits_h2o_broadener(h2o_setting) &&
                    _omit_abundance_term!(omitted, i_band, 1, _H2O_SELF_BROADENING)
                if compute_h2o_jacobians
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], τ̇_abs[i_band], jac_idx,
                        h2o_setting, params.spec_bands[i_band], profile;
                        pressure_tangent, vertical_grid)
                else
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], h2o_setting,
                        params.spec_bands[i_band], profile;
                        pressure_tangent, vertical_grid)
                end
            else
                @timeit "Read HITRAN" lines_h2o = _hitran_lines("H2O", FT)
                @debug "Computing profile for water vapor (q-driven) in band #$(i_band)"
                absorption_model = AtmosphericAbsorption.LineByLineModel(lines_h2o;
                    profile = abs_params.broadening_function,
                    wing_cutoff = abs_params.wing_cutoff,
                    cpf = abs_params.CEF,
                    architecture = _to_aa_arch(params.architecture),
                    vmr = 0)
                h2o_tangent && _omits_h2o_broadener(absorption_model) &&
                    _omit_abundance_term!(omitted, i_band, 1, _H2O_SELF_BROADENING)
                if compute_h2o_jacobians
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], τ̇_abs[i_band], jac_idx,
                        absorption_model, params.spec_bands[i_band], profile;
                        pressure_tangent, vertical_grid)
                else
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], absorption_model,
                        params.spec_bands[i_band], profile;
                        pressure_tangent, vertical_grid)
                end
            end
        end

        if !isnothing(abs_params)
            if !isempty(abs_params.fixed_molecules[i_band])
                for molec_i in 1:length(abs_params.fixed_molecules[i_band])
                    mol_name = abs_params.fixed_molecules[i_band][molec_i]
                    if isempty(abs_params.luts)
                        @timeit "Read HITRAN" lines = _hitran_lines(mol_name, FT)

                        @debug "Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)"
                        absorption_model = AtmosphericAbsorption.LineByLineModel(lines;
                            profile = abs_params.broadening_function,
                            wing_cutoff = abs_params.wing_cutoff,
                            cpf = abs_params.CEF,
                            architecture = _to_aa_arch(params.architecture),
                            vmr = 0)
                        @timeit "Absorption Coeff"  compute_absorption_profile!(
                                τ_abs[i_band],
                                absorption_model,
                                params.spec_bands[i_band],
                                profile.vmr[mol_name],
                                profile; pressure_tangent, vertical_grid)
                    else
                        compute_absorption_profile!(τ_abs[i_band],
                            abs_params.luts[i_band][molec_i],
                            params.spec_bands[i_band],
                            profile.vmr[mol_name],
                            profile; pressure_tangent, vertical_grid)
                        h2o_tangent && _profile_broadener_tangent!(τ̇_abs[i_band],
                            omitted, i_band, abs_params.luts[i_band][molec_i], mol_name,
                            params.spec_bands[i_band], profile.vmr[mol_name], profile)
                    end
                end
            end
            if !isempty(abs_params.variable_molecules[i_band])
                # luts[i_band] is parallel to vcat(fixed_molecules[i_band],
                # variable_molecules[i_band]); offset variable indices past fixed.
                lut_offset = length(abs_params.fixed_molecules[i_band])
                for molec_i in 1:length(abs_params.variable_molecules[i_band])
                    mol_name = abs_params.variable_molecules[i_band][molec_i]
                    jac_idx = molec_i + 1
                    if isempty(abs_params.luts)
                        @timeit "Read HITRAN" lines = _hitran_lines(mol_name, FT)
                        @debug "Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)"

                        absorption_model = AtmosphericAbsorption.LineByLineModel(lines;
                            profile = abs_params.broadening_function,
                            wing_cutoff = abs_params.wing_cutoff,
                            cpf = abs_params.CEF,
                            architecture = _to_aa_arch(params.architecture),
                            vmr = 0)
                        @timeit "Absorption Coeff"  compute_absorption_profile!(
                                τ_abs[i_band],
                                τ̇_abs[i_band],
                                jac_idx,
                                absorption_model,
                                params.spec_bands[i_band],
                                profile.vmr[mol_name],
                                profile; pressure_tangent, vertical_grid)
                    else
                        compute_absorption_profile!(
                            τ_abs[i_band],
                            τ̇_abs[i_band],
                            jac_idx,
                            abs_params.luts[i_band][lut_offset + molec_i],
                            params.spec_bands[i_band],
                            profile.vmr[mol_name],
                            profile; pressure_tangent, vertical_grid)
                        h2o_tangent && _profile_broadener_tangent!(τ̇_abs[i_band],
                            omitted, i_band, abs_params.luts[i_band][lut_offset + molec_i],
                            mol_name, params.spec_bands[i_band], profile.vmr[mol_name],
                            profile)
                    end
                end
            end
        end

        # At fixed final-grid T, q, and VMR, dτ/dp_surf contains both the
        # cross-section response (accumulated above) and molecular-column
        # response. Which layers respond depends on `vertical_grid`: only the
        # bottom layer (FixedPressureGrid) or every layer (FixedAltitudeGrid).
        # CIA and MT_CKD receive their separate density scaling below.
        if compute_pressure_jacobians
            _add_column_psurf_tangent!(τ̇_abs_psurf[i_band], τ_abs[i_band],
                                       vertical_grid, psurf_tangents, profile)
        end

        # Collision-induced absorption (HITRAN .cia files), if any. CIA is
        # held fixed in the gas-profile Jacobian τ̇_abs; a collider abundance
        # it depends on is recorded as omitted. Its surface-pressure tangent
        # is exact at fixed T and composition: τ_CIA ∝ n_midpoint * N_column.
        if !isnothing(abs_params)
            for (cia_i, cia_path) in enumerate(abs_params.cia_files)
                @timeit "CIA $(basename(cia_path))" begin
                    cia_table = _load_configured_cia_table(
                        abs_params, cia_i, params.spec_bands[i_band], FT)
                    τ_cia = zeros(eltype(τ_abs[i_band]), size(τ_abs[i_band]))
                    Absorption.compute_τ_cia!(τ_cia, cia_table, profile,
                                               profile.vmr)
                    τ_abs[i_band] .+= τ_cia
                    # Every pair depends on H₂O through the dry-air number
                    # density; dry colliders through their own VMR.
                    if _cia_has_cross_section(cia_table)
                        term = "$(cia_table.pair) collision-induced absorption"
                        h2o_tangent && _omit_abundance_term!(omitted, i_band, 1, term)
                        for (i, mol) in enumerate(abs_params.variable_molecules[i_band])
                            mol in (cia_table.species_a, cia_table.species_b) &&
                                _omit_abundance_term!(omitted, i_band, 1 + i, term)
                        end
                    end
                    if compute_pressure_jacobians
                        _add_binary_psurf_tangent!(τ̇_abs_psurf[i_band], τ_cia,
                            vertical_grid, psurf_tangents, profile)
                    end
                end
            end

            # MT_CKD H₂O continuum, if a reference table is configured.
            # It is held fixed in the H₂O-profile Jacobian (recorded as
            # omitted), but τ_cont ∝ p_midpoint * N_H2O gives the same exact
            # pressure factor as CIA when q and T are held fixed.
            if !isempty(abs_params.mtckd_file)
                @timeit "MT_CKD H2O continuum" begin
                    mtckd_table = Absorption.load_mtckd(abs_params.mtckd_file)
                    τ_continuum = zeros(
                        eltype(τ_abs[i_band]), size(τ_abs[i_band]))
                    Absorption.compute_τ_h2o_continuum!(
                        τ_continuum, mtckd_table, params.spec_bands[i_band],
                        profile, profile.vmr_h2o)
                    τ_abs[i_band] .+= τ_continuum
                    # With q nonzero, a zero continuum means no coefficients
                    # on this band grid.
                    h2o_tangent && any(!iszero, τ_continuum) &&
                        _omit_abundance_term!(omitted, i_band, 1,
                            "MT_CKD H₂O self and foreign continuum")
                    if compute_pressure_jacobians
                        _add_binary_psurf_tangent!(τ̇_abs_psurf[i_band], τ_continuum,
                            vertical_grid, psurf_tangents, profile)
                    end
                end
            end
        end
    end

    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands]
    lin_aerosol_optics = [Array{linAerosolOptics}(undef, (n_aer)) for i=1:n_bands]

    # Aerosol optics use the model float type; `FT2` (the gas-tangent storage
    # type, possibly a Dual) is kept intact for the coercion at the end.
    FT_aer = params.float_type

    τ_aer = [zeros(FT_aer, n_aer, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands]
    τ̇_aer = [zeros(FT_aer, n_aer, 7, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands]

    for i_aer=1:n_aer

        c_aero = scat.rt_aerosols[i_aer]
        τ_ref = c_aero.τ_ref
        curr_aerosol = c_aero.aerosol

        size_distribution = curr_aerosol.size_distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)

        # Linearized Mie has NO GPU variant: ForwardDiff/analytic-derivative Greek
        # coefficients are only implemented on the CPU analytic kernel. Force CPU
        # architecture here regardless of params.architecture so a GPU run still
        # builds Jacobians correctly via compute_aerosol_optical_properties(lin, …).
        _mie(λ) = make_mie_model(scat.decomp_type, mie_aerosol, λ,
            params.polarization_type, truncation_type, scat.r_max, scat.nquad_radius;
            r_min = scat.r_min, architecture = Architectures.CPU())

        # Use the same AOD convention as the forward constructor:
        # q(λ) = k(λ; n, size) / k_ref(λ_ref; n_ref, size).
        # n_ref is fixed configuration (even when the parser initialized it
        # from aerosol 1). Thus ∂k_ref/∂n = 0, while its size derivatives
        # remain in the quotient rule. See the aerosols schema for the contract.
        ref_aerosol = Aerosol(size_distribution, real(scat.n_ref), -imag(scat.n_ref))
        mie_model_ref = make_mie_model(
            scat.decomp_type, ref_aerosol, scat.λ_ref,
            params.polarization_type, truncation_type,
            scat.r_max, scat.nquad_radius; r_min=scat.r_min,
            architecture=Architectures.CPU())
        if compute_aerosol_microphysics_jacobians
            k_ref, k̇_ref = compute_ref_aerosol_extinction(
                lin, mie_model_ref, params.float_type)
            k̇_ref[1:2] .= 0  # The reference index is not a retrieved aerosol index.
        else
            k_ref = compute_ref_aerosol_extinction(mie_model_ref, params.float_type)
            k̇_ref = zeros(FT_aer, 4)
        end

        function aerosol_pair(λ)
            mie_model = _mie(λ)
            if compute_aerosol_microphysics_jacobians
                return compute_aerosol_optical_properties(lin, mie_model, FT_aer)
            end
            aerosol = compute_aerosol_optical_properties(mie_model)
            return aerosol, _zero_aerosol_microphysics_jacobian(aerosol, FT_aer)
        end

        for i_band=1:n_bands

            curr_band_λ = params.float_type(1e4) ./ params.spec_bands[i_band]
            anchor_ν = FT_aer.(aerosol_bands[i_band])
            anchor_band_λ = FT_aer(1e4) ./ anchor_ν
            already_truncated = false

            if length(curr_band_λ) == 1 && length(anchor_band_λ) == 1
                @timeit "Mie calc"  aerosol_optics_raw, lin_aerosol_optics_raw =
                                aerosol_pair((maximum(curr_band_λ)+minimum(curr_band_λ))/2)

            else
                # Multi-spectral band: compute Mie properties at band edges (λ[1], λ[end]) and
                # linearly interpolate k_ext, k_sca and their derivatives in wavenumber across
                # the band. Greek coefficients are averaged. This avoids a full Mie calculation
                # at every spectral point (Mie varies smoothly with λ).
                n_spec = length(curr_band_λ)
                length(anchor_band_λ) > 1 || throw(ArgumentError(
                    "a multi-wavelength solve grid requires at least two aerosol anchors"))
                @timeit "Mie calc"  aerosol_optics_raw_0, lin_aerosol_optics_raw_0 =
                                aerosol_pair(anchor_band_λ[1])

                @timeit "Mie calc"  aerosol_optics_raw_1, lin_aerosol_optics_raw_1 =
                                aerosol_pair(anchor_band_λ[end])

                function truncate_endpoint(a, la)
                    if truncation_type isa Scattering.δBGE &&
                       length(a.greek_coefs.β) > truncation_type.l_max
                        return Scattering.truncate_phase(
                            truncation_type, a, la; reportFit=false)
                    end
                    return Scattering.truncate_phase(
                        Scattering.NoTruncation(), a, la)
                end
                endpoint₀, lin_endpoint₀ = truncate_endpoint(
                    aerosol_optics_raw_0, lin_aerosol_optics_raw_0)
                endpoint₁, lin_endpoint₁ = truncate_endpoint(
                    aerosol_optics_raw_1, lin_aerosol_optics_raw_1)
                ν_spec = FT_aer.(params.spec_bands[i_band])
                ν_ref_phase = FT_aer(1e4) / scat.λ_ref
                endpoint_ref = lin_endpoint_ref = nothing
                if _interior_reference_node(anchor_ν, ν_ref_phase)
                    @timeit "Mie calc" aerosol_optics_raw_ref, lin_aerosol_optics_raw_ref =
                        aerosol_pair(scat.λ_ref)
                    endpoint_ref, lin_endpoint_ref = truncate_endpoint(
                        aerosol_optics_raw_ref, lin_aerosol_optics_raw_ref)
                end
                aerosol_optics_raw, lin_aerosol_optics_raw =
                    _spectralize_truncated_endpoints(
                        endpoint₀, lin_endpoint₀, endpoint₁, lin_endpoint₁, ν_spec;
                        reference=endpoint_ref, lin_reference=lin_endpoint_ref,
                        ν_ref=ν_ref_phase,
                        ν_endpoints=(first(anchor_ν), last(anchor_ν)))
                already_truncated = true
            end

            # Always go through `truncate_phase`. For NoTruncation this
            # is the cheap passthrough that resets the raw Mie sentinel
            # `fᵗ = 1` to 0 (and `ḟᵗ` to zero) so downstream
            # `delta_m_truncation_lin` doesn't silently kill the
            # aerosol SSA. For δBGE the existing least-squares fit runs
            # only when the Greek expansion is long enough to truncate.
            β_len = size(aerosol_optics_raw.greek_coefs.β, 1)
            if already_truncated
                aerosol_optics[i_band][i_aer] = aerosol_optics_raw
                lin_aerosol_optics[i_band][i_aer] = lin_aerosol_optics_raw
                l_max_aer[i_aer, i_band] = β_len
            elseif truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                aerosol_optics[i_band][i_aer], lin_aerosol_optics[i_band][i_aer] =
                    Scattering.truncate_phase(truncation_type,
                                aerosol_optics_raw, lin_aerosol_optics_raw; reportFit=false)
                l_max_aer[i_aer, i_band] = truncation_type.l_max
            else
                aerosol_optics[i_band][i_aer], lin_aerosol_optics[i_band][i_aer] =
                    Scattering.truncate_phase(Scattering.NoTruncation(),
                                aerosol_optics_raw, lin_aerosol_optics_raw)
                l_max_aer[i_aer, i_band] = β_len
            end

            k_band = aerosol_optics[i_band][i_aer].k
            k̇_band = lin_aerosol_optics[i_band][i_aer].k̇
            ν_spec = FT_aer.(params.spec_bands[i_band])
            ν_ref = FT_aer(1e4) / scat.λ_ref
            aod_scale = _aod_spectral_scale(ν_spec, k_band, k_ref, ν_ref)
            aod_scale_dot = zeros(FT_aer, 4, length(ν_spec))
            # Interpolation is already part of k_band and k̇_band. Apply the
            # same quotient rule at every wavelength, including λ_ref.
            for ctr in 1:4
                aod_scale_dot[ctr, :] .= k̇_band[ctr, :] ./ k_ref .-
                    k_band .* k̇_ref[ctr] ./ k_ref^2
            end

            # Match the production forward profile discretization exactly.
            # Normal profiles use pressure; altitude-form LogNormal profiles
            # use exact CDF differences at geometric layer interfaces.
            τₚ, dτₚdp₀, dτₚdσp =
                getAerosolLayerOptProp(lin, 1, c_aero.profile, profile)

            # ────────────────────────────────────────────────────────────────
            # Aerosol optical depth per layer:
            #   τ_aer(λ,z) = (τ_ref / k_ref) · k(λ) · τₚ(z)
            #
            # Jacobian parameters [1..7]:
            #   1: τ_ref    → ∂τ_aer/∂τ_ref = (k/k_ref) · τₚ
            #   2-5: nᵣ,nᵢ,rₘ,σᵣ → quotient rule on k(λ)/k_ref:
            #        ∂τ_aer/∂p = τ_ref·(dk/dp/k_ref − k·dk_ref/dp/k_ref²)·τₚ
            #   6-7: p₀,σₚ  → only τₚ(z) depends on these:
            #        ∂τ_aer/∂p = (τ_ref/k_ref)·k(λ)·∂τₚ/∂p
            # ────────────────────────────────────────────────────────────────
            τ_aer[i_band][i_aer,:,:] = τ_ref .* aod_scale .* τₚ'
            if compute_pressure_jacobians
                dτₚdpsurf = aerosol_profile_psurf_tangent(
                    vertical_grid, c_aero.profile, profile, psurf_tangents)
                τ̇_aer_psurf[i_band][i_aer, :, :] .=
                    τ_ref .* aod_scale .* dτₚdpsurf'
            end

            τ̇_aer[i_band][i_aer,1,:,:] .=
                aod_scale .* τₚ'

            for ctr=1:4
                τ̇_aer[i_band][i_aer,ctr+1,:,:] =
                    τ_ref .* aod_scale_dot[ctr, :] .* τₚ'
            end
            for ctr=5:6
                τ̇_aer[i_band][i_aer,ctr+1,:,:] =
                    τ_ref .* aod_scale .*
                    (ctr==5 ? dτₚdp₀' : dτₚdσp')
            end

        end
    end
    # Mirror forward path: with no aerosols, fall back to the legacy l_trunc.
    for i_band = 1:n_bands
        l_max[i_band] = n_aer > 0 ? maximum(l_max_aer[:,i_band]) : params.l_trunc
    end
    if n_aer > 0
        set_uniform_lmax!(l_max, aerosol_optics, lin_aerosol_optics)
    end

    # Per-band Fourier loop bound (order). Phase C: trait-based aggregator
    # via `component_m_max(c, ctx)`. Same helper as the forward path so
    # forward and lin can never silently disagree.
    components_per_band = [_band_components(params, aerosol_optics, sources, i_band;
                                rayleigh_active=_rayleigh_active(τ_rayl[i_band]))
                            for i_band in 1:n_bands]
    m_max_bands = _derive_m_max_bands_via_traits(l_max, params.max_m,
                                                  components_per_band,
                                                  quad_points.Nstreams;
                                                  greek_beta_cutoff=params.greek_beta_cutoff)
    n_fourier_moments_bands = m_max_bands .+ 1

    # Build the hierarchical RTModel
    solver = SolverConfig{FT, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        m_max_bands,
        n_fourier_moments_bands,
        l_max,
        params.l_trunc,
        FT(params.Δ_angle),
        FT(params.depol),
        true,   # use_component_traits — flipped on in Phase C
    )
    spec_bands_ft = [convert(Vector{FT}, b) for b in params.spec_bands]
    atm = Atmosphere(profile, spec_bands_ft)
    rayleigh_s = RayleighScattering(greek_rayleigh, greek_cabannes, FT.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    # τ_abs may have been allocated at FT2 (= eltype of the VMR for AD)
    # which can differ from FT (= params.float_type, e.g. Float32 with
    # Float64 VMR). Optics requires both τ-vectors to share the same
    # element type, so coerce τ_abs back to FT here. Matches the forward
    # path's behavior of always allocating τ_abs at params.float_type.
    τ_abs_ft = eltype(eltype(τ_abs)) === FT ? τ_abs : [convert(Matrix{FT}, m) for m in τ_abs]
    optics = Optics(rayleigh_s, aerosols_s, τ_abs_ft, τ_rayl)
    numerics = _convert_numerics(params.numerics, FT)
    model = RTModel(params.architecture, solver, numerics, obs_geom, quad_points, atm, optics, params.brdf, sources)
    # Tangent storage type: FT2 exists so ForwardDiff Duals in the VMR can
    # flow through the tangents — preserve it in that case. When FT2 is a
    # plain float that merely differs from the model float type (Float64
    # YAML VMR in a Float32 model), coerce the stored tangents to FT so a
    # Float32/GPU model never carries Float64 tangent arrays into kernels.
    _to_model_ft(x) = (FT2 <: AbstractFloat && FT2 !== FT) ?
        map(a -> convert(Array{FT, ndims(a)}, a), x) : x
    return model, RTModelLin(_to_model_ft(τ̇_abs), _to_model_ft(τ̇_aer), lin_aerosol_optics,
                             compute_pressure_jacobians ? _to_model_ft(τ̇_rayl_psurf) : nothing,
                             compute_pressure_jacobians ? _to_model_ft(τ̇_aer_psurf) : nothing,
                             compute_pressure_jacobians ? _to_model_ft(τ̇_abs_psurf) : nothing,
                             AbundanceTangentContract(abundance_jacobians, gas_species, omitted))
end

"""
    model_from_parameters_lin(params; kwargs...)

Convenience alias for `model_from_parameters(LinMode(), params; kwargs...)`.
Returns `(model, lin_model)`.
"""
model_from_parameters_lin(params::vSmartMOM_Parameters; kwargs...) =
    model_from_parameters(LinMode(), params; kwargs...)
