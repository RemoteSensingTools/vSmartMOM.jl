#=

This file contains the `model_from_parameters` function, which computes all derived information
like optical thicknesses, from the input parameters. Produces an RTModel object.

=#

"Generate default set of parameters for Radiative Transfer calculations (from ModelParameters/)"
default_parameters() = vSmartMOM.IO.parameters_from_yaml(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))

"Load one CIA table with the reference and negative-value policy parsed for it."
function _load_configured_cia_table(ap::AbsorptionParameters, cia_i::Integer,
                                    ν_grid::AbstractVector,
                                    ::Type{FT}) where {FT}
    return Absorption.load_cia_table(
        ap.cia_files[cia_i], ν_grid;
        FT,
        reference_codes=ap.cia_reference_codes[cia_i],
        negative_policy=ap.cia_negative_policies[cia_i])
end

"""
    _resolved_truncation(params, FT)

Return `params.truncation` typed for `FT`. If a downstream caller has
mutated the legacy `l_trunc` / `Δ_angle` fields after construction (an
older idiom — e.g. `test_float32.jl`), they must also mutate
`params.truncation` to keep the two in sync; this function does NOT
silently rebuild from the legacy fields, because that would discard a
user-set explicit truncation (e.g. `NoTruncation()` or a custom δBGE).
"""
@inline _resolved_truncation(params, ::Type{FT}) where {FT} =
    _resolve_auto_truncation(params.truncation, params, FT)

"""
    _resolve_auto_truncation(t, params, FT) -> AbstractTruncationType

Resolve `Scattering.AutoTruncation()` (the `truncation: auto` YAML
knob) into a concrete truncation type at model-build time:

- **No aerosols** → `NoTruncation()`. Rayleigh has only β₀, β₁, β₂
  so the phase function fits any `stream_l_cap ≥ 5`.
- **With aerosols** → `δBGE(stream_l_cap, Δ_angle)`. Mie phase
  functions can have hundreds of Greek coefficients; the δ-BGE
  forward-peak fit (Sanghavi & Stephens 2015) keeps the projection
  tractable. The Mie call sites in `model_from_parameters` guard
  the `truncate_phase(::δBGE, …)` call with a length check
  (`length(greek.β) > l_max`) so short phase functions degrade
  to `NoTruncation()` automatically — no crash on small-particle
  Mie that produces a series shorter than `stream_l_cap`.

Mirrors VLIDORT's `DO_DELTAM_SCALING` philosophy. Logs the chosen
branch via `@info`. Non-`AutoTruncation` inputs pass through
unchanged so users keep full control.
"""
function _resolve_auto_truncation(t::Scattering.AbstractTruncationType, params, ::Type{FT}) where {FT}
    t isa Scattering.AutoTruncation || return t

    has_aerosols = params.scattering_params !== nothing &&
                   !isempty(params.scattering_params.rt_aerosols)
    if has_aerosols
        l_cap = params.stream_l_cap > 0 ? params.stream_l_cap : params.l_trunc
        chosen = Scattering.δBGE{FT}(l_cap, FT(params.Δ_angle))
        @info "truncation: auto → δBGE(stream_l_cap=$(l_cap), Δ_angle=$(FT(params.Δ_angle))) (per-band greek length checked at Mie call site; falls back to NoTruncation if phase function already fits)"
        return chosen
    else
        @info "truncation: auto → NoTruncation() (no aerosols; Rayleigh-only scene)"
        return Scattering.NoTruncation()
    end
end

_has_analytic_phase_function(c_aero::RT_Aerosol) =
    c_aero.phase_function !== nothing

"""
    _derive_m_max_bands(l_max::Vector{Int}, legacy_max_m_count::Int) -> Vector{Int}

Per-band Fourier loop bound (**order semantics**) for both the forward
and linearized RT paths. Identical aggregator across both paths so the
loops can never silently diverge — see Phase B of the v0.7
fourier-stream-resolution refactor.

The inner-paren division `Int(ceil((l+1)/2))` is the **forward path's**
historical formula. The lin path used to write
`Int(ceil(l+1)/2)` (outer-paren division), which evaluates `Int` of a
non-integer Float at even `l_max` and either throws `InexactError` or
silently rounds down by 1 depending on Julia version. Phase B unifies
both paths through this helper, fixing the latent precedence bug.

Returns a vector of orders `m_max_bands[i] = ceil((l_max[i]+1)/2) - 1`,
clamped from above by the user-supplied `legacy_max_m_count - 1`.
"""
function _derive_m_max_bands(l_max::AbstractVector{Int}, legacy_max_m_count::Int)
    m_max_bands = similar(l_max, Int)
    legacy_order_cap = max(legacy_max_m_count - 1, 0)
    for i in eachindex(l_max)
        m_max_bands[i] = Int(ceil((l_max[i] + 1) / 2)) - 1
        m_max_bands[i] = min(m_max_bands[i], legacy_order_cap)
    end
    return m_max_bands
end

"""
    _derive_m_max_bands_via_traits(l_max, legacy_max_m_count, components_per_band) -> Vector{Int}

Phase C aggregator (active when `SolverConfig.use_component_traits == true`).
Computes per-band `m_max_bands` (order) by taking
`maximum(component_m_max(c, ctx))` across the band's active components and
clamping to `min(stream_l_cap, legacy_l_cap)` where:

- `stream_l_cap = 2·l_max - 1` (the projection cap implied by the chosen
  Legendre truncation; conservative since `l_max` is itself capped by
  the user's `l_trunc`)
- `legacy_l_cap = legacy_max_m_count - 1` (user's `max_m` knob, in order
  semantics)

For aerosol components, `component_m_max` optionally applies
`greek_beta_cutoff` after the size distribution has been integrated into Greek
coefficients: it returns the last degree whose `abs(β_l)` reaches the cutoff at
any wavelength in the band. Thus the maximum-size relationship remains the
allocation ceiling, while β supplies a second effective-support filter. No
other Greek coefficient family participates.

`components_per_band[i]` is an iterable of components participating in
band `i`. The trait dispatch already lives in `component_m_max.jl`.
"""
function _derive_m_max_bands_via_traits(l_max::AbstractVector{Int},
                                        legacy_max_m_count::Int,
                                        components_per_band,
                                        nstreams::Int;
                                        greek_beta_cutoff=nothing)
    @assert length(l_max) == length(components_per_band)
    m_max_bands = similar(l_max, Int)
    # Public contract: stream_l_cap = 2·Nstreams - 1 (see Phase A; same
    # formula for both Gauss and Radau — Radau's internal subinterval
    # node count is bookkeeping). Codex review of Phase C (P2) flagged
    # that using `2·l_max - 1` here lets the trait aggregator hand
    # back loop bounds the quadrature can't actually resolve.
    quadrature_l_cap = max(2 * nstreams - 1, 0)
    legacy_order_cap = max(legacy_max_m_count - 1, 0)
    for i in eachindex(l_max)
        # `user_l_cap` is the most restrictive of: the stream resolving
        # power, the legacy `max_m` knob, and the per-band Legendre
        # ceiling derived from aerosol greek length.
        user_l_cap = min(quadrature_l_cap, legacy_order_cap, l_max[i])
        # Aerosol traits may tighten this further using the already-integrated
        # beta tail; the cutoff does not alter quadrature construction here.
        ctx = (; user_l_cap, stream_l_cap = quadrature_l_cap,
                m_max_override = nothing, truncation = nothing,
                greek_beta_cutoff)
        m_max_bands[i] = _aggregate_m_max(components_per_band[i], ctx)
    end
    return m_max_bands
end

"""
    _band_components(params, aerosol_optics, sources, i_band) -> Vector

Build the per-band component list consumed by the trait aggregator.
Always includes Rayleigh (m_max=2 floor) plus the band's surface BRDF,
any truncated aerosol optics, and the active source(s) — typically a
`SolarBeam` (m_max=0, neutral) but may include `SurfaceSIF` or a
future thermal source whose Fourier support is non-trivial.

The Rayleigh contribution is added as the type itself
(`RayleighScattering`) and `component_m_max(::Type{...})` is dispatched
to `2` — this avoids constructing a dummy `RayleighScattering`
instance at trait-resolution time.

Codex review of Phase C (P2) flagged that omitting `sources` here
would silently drop source-driven Fourier moments for any future
source whose trait isn't `0`.
"""
function _band_components(params, aerosol_optics, sources, i_band)
    comps = Any[RayleighScattering, params.brdf[i_band], sources]
    if !isempty(aerosol_optics) && i_band <= length(aerosol_optics)
        for ao in aerosol_optics[i_band]
            push!(comps, ao)
        end
    end
    return comps
end

function _finite_truncation_lmax(params, truncation_type)
    fallback = max(1, Int(params.l_trunc), 2 * Int(params.max_m) - 1)
    hasproperty(truncation_type, :l_max) || return fallback
    trunc_lmax = Int(getproperty(truncation_type, :l_max))
    return trunc_lmax == typemax(Int) ? fallback : max(1, trunc_lmax)
end

function _analytic_phase_lmax(params, truncation_type)
    trunc_lmax = _finite_truncation_lmax(params, truncation_type)
    base_lmax = max(1, trunc_lmax, 2 * Int(params.max_m) - 1)
    if truncation_type isa Scattering.δBGE
        return max(base_lmax + 1, 2 * trunc_lmax)
    end
    return base_lmax
end

function _analytic_aerosol_optics(c_aero::RT_Aerosol, params, truncation_type,
                                  ::Type{FT}) where {FT}
    lmax = _analytic_phase_lmax(params, truncation_type)
    return Scattering.analytic_aerosol_optics(
        c_aero.phase_function;
        single_scattering_albedo = convert(FT, c_aero.ϖ),
        extinction_cross_section = one(FT),
        l_max = lmax,
        nquad = max(2lmax + 1, 64),
        FT_out = FT)   # emit greek in the model FT (type stability; mirrors the Mie path)
end

"""
    model_from_parameters(params::vSmartMOM_Parameters;
                          sources=SolarBeam(), external_solar=false) -> RTModel

Construct an [`RTModel`](@ref) from user-supplied [`vSmartMOM_Parameters`](@ref).

Computes all derived quantities needed by the RT solver:
- Observation geometry and quadrature points
- Rayleigh scattering coefficients (Greek/Cabannes)
- Aerosol optics via Mie theory with δ-M truncation
- Gas absorption cross-sections and optical depth profiles
- Per-band Fourier truncation limits

# Arguments
- `params::vSmartMOM_Parameters`: Configuration struct (typically from `parameters_from_yaml`).
- `external_solar::Bool=false`: When `true`, keep the direct solar direction
  outside the diffuse operator (rectangular `Z₀` source columns instead of an
  embedded zero-weight μ₀ node). This is an opt-in, TOA-only fast path — it
  requires Gauss quadrature and is consumed through [`rt_run_toa`](@ref).
  The default embedded-μ₀ mode supports every entry point: full `R, T`
  output, interior/BOA observer heights, the atmosphere/surface split, and
  Radau quadrature.

# Returns
- `model::RTModel{ARCH, FT}`: Hierarchical model ready for `rt_run(model)`.

# See also
- `model_from_parameters(LinMode(), params)` for the linearized (Jacobian) variant.
- `parameters_from_yaml(path)` to load parameters from a YAML file.
"""
function model_from_parameters(params::vSmartMOM_Parameters;
                               sources::AbstractSource = SolarBeam(),
                               external_solar::Bool = false)
    FT = params.float_type
    #@show FT
    # Number of total bands and aerosols (for convenience)
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    # Truncation method (typed; NoTruncation, δBGE, ...). The legacy
    # `params.Δ_angle` is only consulted via the default δBGE built in
    # `parameters_from_dict` when the user has not set `truncation`
    # explicitly. If a downstream caller mutates `l_trunc` / `Δ_angle`
    # after construction, they must also mutate `params.truncation` —
    # see `_resolved_truncation`.
    truncation_type = _resolved_truncation(params, FT)
    # Get AtmosphericProfile from parameters.
    # Convert T, p, q to the requested FT so that profile arrays (vcd_dry,
    # p_half, etc.) are consistent with the model float type.  Without this,
    # mutating params.float_type = Float32 on a Float64-constructed
    # vSmartMOM_Parameters causes getRayleighLayerOptProp to receive mixed
    # Float32/Float64 arguments and raises a MethodError.
    vmr = isnothing(params.absorption_params) ? Dict() : params.absorption_params.vmr
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

    # Rayleigh optical depth per spectral point per layer (uses reduced profile size).
    #
    # Depolarization sourcing rule (2026-04-24):
    #   params.depol < 0  → "auto": derive both Rayleigh (greek_rayleigh, τ_rayl)
    #                       and Cabannes (greek_cabannes) depolarizations per
    #                       band from the N₂/O₂ molecular constants (sanghavi
    #                       convention; appropriate for Earth atmospheres).
    #   params.depol ≥ 0  → "explicit": use params.depol uniformly for
    #                       greek_rayleigh, greek_cabannes, AND τ_rayl. Required
    #                       for idealizations such as Natraj (1980) which fix
    #                       depol = 0 — no molecular path can reproduce that.
    # ϖ_Cabannes is always taken from the molecular path (it relates Cabannes
    # single-scattering albedo to N₂/O₂ line ratios, not depol).
    τ_rayl = [zeros(FT,length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];

    # Per-band Cabannes / Rayleigh greek coefs (depolarizations from molecular constants)
    greek_rayleigh_arr = Vector{Scattering.GreekCoefs{FT}}()
    greek_cabannes = Vector{Scattering.GreekCoefs{FT}}()
    ϖ_Cabannes = zeros(FT, n_bands)
    
    τ_abs     = [zeros(FT, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    
    # Track per-band l_max from aerosol greek coef lengths
    l_max_aer = zeros(Int, max(n_aer, 1), n_bands)
    
    # Loop over all bands:
    for i_band=1:n_bands

        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])
        
        # Per-band molecular-constant depolarizations (always needed for ϖ_Cab;
        # used as the depol values when params.depol < 0).
        νₘ = FT(0.5) * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ = FT(1.0e7) / νₘ
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(FT(1.0e7) / λₘ, FT(rayleigh_molecular_T))
        ϖ_Cab = InelasticScattering.compute_ϖ_Cabannes(λₘ, _n2, _o2)
        γ_air_Cab, _ = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
        γ_air_Ray, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ, _n2, _o2)
        ϖ_Cabannes[i_band] = FT(ϖ_Cab)
        depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
        depol_air_Ray = 2γ_air_Ray / (1 + γ_air_Ray)

        # Apply the auto vs. explicit depol rule (see comment above the loop).
        depol_use_Cab = params.depol < 0 ? FT(depol_air_Cab) : FT(params.depol)
        depol_use_Ray = params.depol < 0 ? FT(depol_air_Ray) : FT(params.depol)

        push!(greek_cabannes,     Scattering.get_greek_rayleigh(depol_use_Cab))
        push!(greek_rayleigh_arr, Scattering.get_greek_rayleigh(depol_use_Ray))

        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end],
                                curr_band_λ,
                                depol_use_Ray, profile.vcd_dry);
        #@show τ_rayl[i_band]
        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        ap = params.absorption_params

        # Loop over fixed_molecules ∪ variable_molecules in this band; H2O is
        # handled separately below (driven by q, not by these lists).
        all_species = vcat(ap.fixed_molecules[i_band], ap.variable_molecules[i_band])
        for (molec_i, mol_name) in enumerate(all_species)
            if isempty(ap.luts)
                @timeit "Read HITRAN" lines = AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(artifact(mol_name)); FT)
                @debug "Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)"
                absorption_model = AtmosphericAbsorption.LineByLineModel(lines;
                    profile = ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    cpf = ap.CEF,
                    architecture = _to_aa_arch(params.architecture),
                    vmr = 0)
                @timeit "Absorption Coeff" compute_absorption_profile!(τ_abs[i_band], absorption_model, params.spec_bands[i_band], profile.vmr[mol_name], profile)
            else
                # TIMED: this branch used to be the only unmeasured stage in the
                # model build, which is how a 19x regression in it went unnoticed
                # (a verbose run reported "Tot / % measured: 11.3%", with the whole
                # LUT absorption cost sitting in the unmeasured remainder).
                @timeit "Absorption Coeff LUT" compute_absorption_profile!(τ_abs[i_band], ap.luts[i_band][molec_i], params.spec_bands[i_band], profile.vmr[mol_name], profile)
            end
        end

        # H₂O line absorption (driven by q). Use the band's H2O LUT if the
        # parser found one inside LUTfiles; otherwise fall back to artifact.
        if any(!iszero, params.q)
            if ap.h2o_lut[i_band] !== nothing
                @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                    τ_abs[i_band], ap.h2o_lut[i_band],
                    params.spec_bands[i_band], profile)
            else
                @timeit "Read HITRAN H2O" lines_h2o = AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(artifact("H2O")); FT)
                @debug "Computing profile for H2O (q-driven) for band #$(i_band)"
                h2o_model = AtmosphericAbsorption.LineByLineModel(lines_h2o;
                    profile = ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    cpf = ap.CEF,
                    architecture = _to_aa_arch(params.architecture),
                    vmr = 0)
                @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                    τ_abs[i_band], h2o_model, params.spec_bands[i_band], profile)
            end
        end

        # Collision-induced absorption (HITRAN .cia files), if any.
        for (cia_i, cia_path) in enumerate(ap.cia_files)
            @timeit "CIA $(basename(cia_path))" begin
                cia_table = _load_configured_cia_table(
                    ap, cia_i, params.spec_bands[i_band], FT)
                Absorption.compute_τ_cia!(τ_abs[i_band], cia_table, profile,
                                           profile.vmr)
            end
        end

        # MT_CKD H₂O continuum (self + foreign), if a reference table is configured.
        if !isempty(ap.mtckd_file)
            @timeit "MT_CKD H2O continuum" begin
                mtckd_table = Absorption.load_mtckd(ap.mtckd_file)
                Absorption.compute_τ_h2o_continuum!(τ_abs[i_band], mtckd_table,
                                                     params.spec_bands[i_band], profile,
                                                     profile.vmr_h2o)
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];

    # τ_aer[iBand][iAer, nSpec, iZ] — 3-D so the forward path matches the lin path:
    # per-spectral aerosol optical depth enables wavelength-dependent k(λ) scaling
    # within wide spectral bands. Single-λ bands remain bit-identical to the old 2-D
    # layout (nSpec=1 is a degenerate case of the 3-D array).
    τ_aer = [zeros(FT, n_aer, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]

        if _has_analytic_phase_function(c_aero)
            τ_profile = getAerosolLayerOptProp(one(FT), c_aero.profile, profile)
            for i_band=1:n_bands
                aerosol_optics_raw =
                    _analytic_aerosol_optics(c_aero, params, truncation_type, FT)
                β_len = length(aerosol_optics_raw.greek_coefs.β)
                aerosol_optics[i_band][i_aer] =
                    if truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                        Scattering.truncate_phase_lowconf(truncation_type,
                                                          aerosol_optics_raw; reportFit=false)
                    else
                        Scattering.truncate_phase(Scattering.NoTruncation(),
                                                  aerosol_optics_raw)
                    end
                l_max_aer[i_aer, i_band] =
                    length(aerosol_optics[i_band][i_aer].greek_coefs.β)
                # Analytic aerosols have no wavelength-dependent k(λ): τ is constant across λ.
                # Write (nSpec × nLayers) slice; τ_profile is an nLayers vector.
                nSpec_band = length(params.spec_bands[i_band])
                τ_aer[i_band][i_aer, :, :] =
                    c_aero.τ_ref .* ones(FT, nSpec_band) .* τ_profile'
                @debug "AOD at band $i_band : $(sum(τ_aer[i_band][i_aer,:,:])), analytic phase function = $(typeof(c_aero.phase_function)), truncation factor = $(aerosol_optics[i_band][i_aer].fᵗ)"
            end
            continue
        end

        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        #@show typeof(curr_aerosol.nᵣ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # k_ref uses the reference refractive index n_ref (a separate parameter that
        # defines the normalisation wavelength/index convention), NOT the aerosol's
        # own nᵣ/nᵢ. We build a separate aerosol object so that mie_aerosol (which
        # carries the aerosol's own refractive index) is never mutated — the in-band
        # Mie closures below must use the aerosol's true nᵣ/nᵢ, not n_ref.
        ref_aerosol    = Aerosol(size_distribution,
                                 real(params.scattering_params.n_ref),
                                 -imag(params.scattering_params.n_ref))
        mie_model_ref  = make_mie_model(params.scattering_params.decomp_type,
                                        ref_aerosol,
                                        params.scattering_params.λ_ref,
                                        params.polarization_type,
                                        truncation_type,
                                        params.scattering_params.r_max,
                                        params.scattering_params.nquad_radius)
        # k for reference wavelength (uses n_ref index, not the aerosol's own index)
        k_ref          = compute_ref_aerosol_extinction(mie_model_ref, params.float_type)

        # Closure for building a Mie model at a given wavelength λ (μm).
        # Uses CPU() to ensure the Greek-coef interpolation below stays on host.
        _mie_fwd(λ) = make_mie_model(params.scattering_params.decomp_type,
                                      mie_aerosol, λ,
                                      params.polarization_type, truncation_type,
                                      params.scattering_params.r_max,
                                      params.scattering_params.nquad_radius;
                                      architecture = params.architecture)

        # Loop over bands
        for i_band=1:n_bands

            # i'th spectral band wavelengths (convert from cm⁻¹ to μm)
            curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])
            n_spec = length(curr_band_λ)

            if n_spec == 1
                # ── Single-wavelength band: Mie at band center (bit-identical to old behavior) ──
                mie_model = _mie_fwd((maximum(curr_band_λ)+minimum(curr_band_λ))/2)
                @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model)

                β_len = length(aerosol_optics_raw.greek_coefs.β)
                aerosol_optics[i_band][i_aer] =
                    if truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                        Scattering.truncate_phase(truncation_type,
                                                  aerosol_optics_raw; reportFit=false)
                    else
                        Scattering.truncate_phase(Scattering.NoTruncation(),
                                                  aerosol_optics_raw)
                    end

                l_max_aer[i_aer, i_band] =
                    truncation_type isa Scattering.δBGE ?
                        min(length(aerosol_optics[i_band][i_aer].greek_coefs.β), truncation_type.l_max) :
                        length(aerosol_optics[i_band][i_aer].greek_coefs.β)

                τ_profile = getAerosolLayerOptProp(1, c_aero.profile, profile)
                τ_aer[i_band][i_aer, 1, :] = vec(_aerosol_τ_slice(
                    params.scattering_params.rt_aerosols[i_aer].τ_ref,
                    aerosol_optics[i_band][i_aer].k, k_ref, τ_profile,
                    FT.(params.spec_bands[i_band]),
                    FT(1e4) / params.scattering_params.λ_ref))
                @debug "AOD at band $i_band (single-λ): $(sum(τ_aer[i_band][i_aer,:,:])), truncation factor = $(aerosol_optics[i_band][i_aer].fᵗ)"

            else
                # ── Multi-wavelength band: Mie at band endpoints, linearly interpolate k(λ) ──
                # Linear endpoint interpolation is an approximation.
                # A future improvement could also compute the band-center point and use
                # quadratic (3-point) interpolation. With smaller spectral bands,
                # in-band aerosol dispersion is negligible and linear interpolation
                # introduces at most ~(Δk/k)²/8 relative error in τ.
                mie_model_0 = _mie_fwd(curr_band_λ[1])
                @timeit "Mie calc"  aerosol_optics_raw_0 = compute_aerosol_optical_properties(mie_model_0)

                mie_model_1 = _mie_fwd(curr_band_λ[end])
                @timeit "Mie calc"  aerosol_optics_raw_1 = compute_aerosol_optical_properties(mie_model_1)

                # Truncate each Mie endpoint independently, then interpolate the
                # resulting physical optics in wavenumber. Because Z is linear in
                # the Greek coefficients, this is exactly endpoint-Z interpolation
                # while remaining independent of the eventual quadrature grid.
                truncate_endpoint(a) =
                    truncation_type isa Scattering.δBGE &&
                    length(a.greek_coefs.β) > truncation_type.l_max ?
                        Scattering.truncate_phase(truncation_type, a; reportFit=false) :
                        Scattering.truncate_phase(Scattering.NoTruncation(), a)
                endpoint₀ = truncate_endpoint(aerosol_optics_raw_0)
                endpoint₁ = truncate_endpoint(aerosol_optics_raw_1)
                ν_spec = FT.(params.spec_bands[i_band])
                ν_ref_phase = FT(1e4) / params.scattering_params.λ_ref
                endpoint_ref = nothing
                if first(extrema(ν_spec)) < ν_ref_phase < last(extrema(ν_spec))
                    @timeit "Mie calc" aerosol_optics_raw_ref =
                        compute_aerosol_optical_properties(_mie_fwd(params.scattering_params.λ_ref))
                    endpoint_ref = truncate_endpoint(aerosol_optics_raw_ref)
                end
                aerosol_optics[i_band][i_aer] =
                    _spectralize_truncated_endpoints(endpoint₀, endpoint₁, ν_spec;
                        reference=endpoint_ref, ν_ref=ν_ref_phase)
                l_max_aer[i_aer, i_band] =
                    size(aerosol_optics[i_band][i_aer].greek_coefs.β, 1)

                τ_profile = getAerosolLayerOptProp(1, c_aero.profile, profile)
                τ_aer[i_band][i_aer, :, :] = _aerosol_τ_slice(
                    params.scattering_params.rt_aerosols[i_aer].τ_ref,
                    aerosol_optics[i_band][i_aer].k, k_ref, τ_profile,
                    ν_spec, FT(1e4) / params.scattering_params.λ_ref)
                @debug "AOD at band $i_band (multi-λ): $(sum(τ_aer[i_band][i_aer,:,:])), truncation factor = $(aerosol_optics[i_band][i_aer].fᵗ)"
            end
        end
    end

    # Compute per-band l_max from aerosol greek coefficient lengths.
    l_max = zeros(Int, n_bands)
    for i_band = 1:n_bands
        if n_aer > 0
            l_max[i_band] = maximum(l_max_aer[:, i_band])
        else
            l_max[i_band] = params.l_trunc
        end
    end
    # Per-band Fourier loop bound (order). Phase C: per-component traits via
    # `component_m_max(c, ctx)` (see src/CoreRT/component_m_max.jl). Each band's
    # component list contains:
    #   - RayleighScattering (always present, contributes m_max=2)
    #   - the band's truncated AerosolOptics list (contributes length(β)-1)
    #   - the band's surface BRDF (contributes 0 for Lambertian, user_l_cap for
    #     Cox-Munk / RossLi / RPV / canopy)
    # Codex review of Phase B (P1) flagged that the previous count-only
    # aggregator silently half-truncated Cox-Munk forward — traits restore the
    # full surface-driven Fourier resolution.
    components_per_band = [_band_components(params, aerosol_optics, sources, i_band)
                            for i_band in 1:n_bands]
    m_max_bands = _derive_m_max_bands_via_traits(l_max, params.max_m,
                                                  components_per_band,
                                                  quad_points.Nstreams;
                                                  greek_beta_cutoff=params.greek_beta_cutoff)
    n_fourier_moments_bands = m_max_bands .+ 1

    # Build the hierarchical RTModel
    FT_ = FT
    solver = SolverConfig{FT_, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        m_max_bands,
        n_fourier_moments_bands,
        l_max,
        params.l_trunc,
        FT_(params.Δ_angle),
        FT_(params.depol),
        true,   # use_component_traits — flipped on in Phase C
    )
    spec_bands_ft = [convert(Vector{FT_}, b) for b in params.spec_bands]
    atm = Atmosphere(profile, spec_bands_ft)
    rayleigh = RayleighScattering(greek_rayleigh_arr, greek_cabannes, FT_.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    optics = Optics(rayleigh, aerosols_s, τ_abs, τ_rayl)
    numerics = _convert_numerics(params.numerics, FT_)
    return RTModel(params.architecture, solver, numerics, obs_geom, quad_points, atm, optics, params.brdf, sources)
end

"""
    remake_geometry(model::RTModel, params::vSmartMOM_Parameters;
                   sza = model.geometry.sza,
                   vza = model.geometry.vza,
                   vaz = model.geometry.vaz) -> RTModel

Rebuild an [`RTModel`](@ref) at a different observation geometry WITHOUT
re-running the atmosphere build (`model_from_parameters`'s HITRAN
absorption, Mie aerosol optics, and per-band Fourier-support derivation) —
the expensive part of `model_from_parameters` for a spectral-band scene.

Only `geometry::ObsGeometry` and `quad_points::QuadPoints` depend on
`(sza, vza, vaz)`: `quad_points` is built by [`rt_set_streams`](@ref),
which appends SZA/VZA cosines as (mostly) zero-weight output nodes onto a
fixed set of weighted quadrature nodes. Everything else — `atmosphere`
(profile, τ_abs, τ_rayl), `optics` (aerosol Mie/greek coefs), `solver`
(m_max_bands, truncation) and `surfaces` — is geometry-independent, so
this function shares them by reference with `model` and only constructs
fresh `ObsGeometry`/`QuadPoints` objects, THE SAME WAY
`model_from_parameters` does (identical `rt_set_streams` call).

`params` must be the SAME parameters object (or an unmutated `deepcopy`)
used to build `model` in the first place — it supplies the quadrature
scheme, `l_trunc`, polarization type and architecture that `model.solver`
was derived from (`params.quadrature_type`, `params.l_trunc`,
`params.polarization_type`, `params.architecture`); this function does
NOT re-derive `solver.m_max_bands` from these, it reuses `model.solver`
verbatim (see the bit-exactness note below for why that's safe). Mutate
`sza`/`vza`/`vaz` via the keyword arguments here, not by mutating
`params` — a `params` with anything OTHER than geometry changed (e.g. a
different `l_trunc`) would silently desync `model.solver` from the
quadrature this function builds.

# Bit-exactness

`rt_run(remake_geometry(model, params; sza=θ))` is bit-exact with
`rt_run(model_from_parameters(p′))` for `p′ = deepcopy(params)` with
`p′.sza = θ` (and same `vza`/`vaz`) — see
`test/test_scenario_sweep.jl`. This relies on `quad_points.Nstreams`
(hence `solver.m_max_bands`, which is NOT recomputed here) being
independent of `sza`/`vza`: `rt_set_streams` builds the weighted
Gauss-Legendre/Radau nodes first (always nonzero weight) and only
afterwards appends the SZA/VZA cosines as additional nodes, padded with
zero weight; `Nstreams = count(!iszero, wt_μ)` therefore only changes if
a SZA/VZA cosine happens to land exactly on an existing weighted node in
floating point (it would then be dropped by `unique` instead of
appended) — a pre-existing knife-edge property of `rt_set_streams`
itself (shared with the monolithic `model_from_parameters` path), not
something this function introduces. `quad_points.Nquad` (the total node
count, including zero-weight SZA/VZA duplicates) MAY legitimately differ
between geometries — that's expected and harmless, since `rt_run` sizes
its per-call `NquadN` from `model.quad_points.Nquad` fresh on every
call, not from anything cached at model-build time.

See `docs/dev_notes/proposals/surface_split_albedo_sweep.md` §6/§7 (PR 3) for the
design context — this is the seam `run_sweep` uses to amortise the
per-SZA HITRAN/Mie rebuild across a scenario sweep.
"""
function remake_geometry(model::RTModel, params::vSmartMOM_Parameters;
                         sza::Real = model.geometry.sza,
                         vza::AbstractVector = model.geometry.vza,
                         vaz::AbstractVector = model.geometry.vaz)
    FT = float_type(model)
    obs_geom = ObsGeometry{FT}(FT(sza), FT.(collect(vza)), FT.(collect(vaz)), model.geometry.obs_alt)
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom,
                                 params.polarization_type, array_type(params.architecture))
    return RTModel(model.architecture, model.solver, model.numerics, obs_geom, quad_points,
                   model.atmosphere, model.optics, model.surfaces, model.sources)
end

"Re-type the user-supplied `RTNumericalParameters` to the resolved
RTModel float type. No-op when types already match."
@inline _convert_numerics(n::RTNumericalParameters{FT}, ::Type{FT}) where {FT} = n
@inline function _convert_numerics(n::RTNumericalParameters, ::Type{FT}) where {FT}
    RTNumericalParameters{FT}(
        dτ_max_threshold = FT(n.dτ_max_threshold),
        dτ_min_floor     = FT(n.dτ_min_floor),
        blas_threads     = n.blas_threads,
        verbose          = n.verbose,
        fourier_convergence = n.fourier_convergence,
        ss_correction    = n.ss_correction,
    )
end


#=

Modified version for vibrational Raman scattering

=#


"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into an RTModel"
function model_from_parameters(RS_type::Union{VS_0to1_plus, VS_1to0_plus},
                    λ₀,
                    params::vSmartMOM_Parameters;
                    sources::AbstractSource = SolarBeam())
    # Number of total bands and aerosols (for convenience)
    n_bands = 3 #length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    # Truncation method (typed; NoTruncation, δBGE, ...). Matches the elastic
    # `model_from_parameters` site: respects user-set `truncation` field and
    # honours `truncation: auto` via `_resolved_truncation`.
    truncation_type = _resolved_truncation(params, params.float_type)

    # Get AtmosphericProfile from parameters.
    # Convert T, p, q to params.float_type so profile arrays are type-consistent.
    FT_vrs_early = params.float_type
    vmr = isnothing(params.absorption_params) ? Dict() : params.absorption_params.vmr
    T_ft, p_ft, q_ft = convert(Vector{FT_vrs_early}, params.T), convert(Vector{FT_vrs_early}, params.p), convert(Vector{FT_vrs_early}, params.q)
    obs_alt = params.obs_alt isa Real ? FT_vrs_early(params.obs_alt) :
              convert(Vector{FT_vrs_early}, params.obs_alt)
    profile, observation = prepare_observer_profile(
        T_ft, p_ft, q_ft, vmr, obs_alt, params.profile_reduction_n)
    input_p_full = (p_ft[1:end-1] .+ p_ft[2:end]) ./ FT_vrs_early(2)
    sources = reframe_vertical_sources(sources, input_p_full, profile.p_full)
    obs_geom = ObsGeometry{FT_vrs_early}(
        FT_vrs_early(params.sza), convert(Vector{FT_vrs_early}, params.vza),
        convert(Vector{FT_vrs_early}, params.vaz), obs_alt,
        observation.sensor_levels, observation.interior_altitudes,
        observation.include_toa, observation.include_boa, observation.toa_altitude)
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom,
                                 params.polarization_type, array_type(params.architecture))

    effT = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry);
    # Define RS type
    # Compute N2 and O2
    RS_type.n2, RS_type.o2 =
        InelasticScattering.getRamanAtmoConstants(1.e7/λ₀,effT);
    # VS-plus signature: getRamanSSProp!(RS, depol, λ_inc) — the depol
    # arg sets the Rayleigh/Raman cross-section ratio per Hovenier convention.
    InelasticScattering.getRamanSSProp!(RS_type, params.depol, λ₀);
    n_bands = length(RS_type.iBand)
    params.spec_bands = RS_type.grid_in

    # Rayleigh optical properties calculation
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τ_rayl = [zeros(params.float_type,1, length(profile.p_full)) for i=1:n_bands];

    # Per-band Cabannes / Rayleigh depolarization (for inelastic scattering support)
    FT_vrs = params.float_type
    greek_cabannes = typeof(greek_rayleigh)[]
    ϖ_Cabannes = zeros(FT_vrs, n_bands)
    l_max_aer = zeros(Int, max(n_aer, 1), n_bands)

    # Pre-allocated absorption arrays
    τ_abs     = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    # Loop over all bands:
    for i_band=1:n_bands
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]

        # Compute per-band Cabannes properties using the same effective
        # temperature as the Raman object for closure between elastic and
        # inelastic fractions.
        νₘ = 0.5 * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ = 1.0e7 / νₘ
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(1.0e7 / λₘ, effT)
        ϖ_Cab = InelasticScattering.compute_ϖ_Cabannes(λₘ, _n2, _o2)
        γ_air_Cab, _ = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
        ϖ_Cabannes[i_band] = FT_vrs(ϖ_Cab)
        depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
        push!(greek_cabannes, Scattering.get_greek_rayleigh(FT_vrs(depol_air_Cab)))

        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                            (maximum(curr_band_λ) + minimum(curr_band_λ))/2, 
                            params.depol, profile.vcd_dry);

        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        ap = params.absorption_params

        # Loop over fixed_molecules ∪ variable_molecules in this band; H2O is
        # handled separately below (driven by q, not by these lists).
        all_species = vcat(ap.fixed_molecules[i_band], ap.variable_molecules[i_band])
        for (molec_i, mol_name) in enumerate(all_species)
            if isempty(ap.luts)
                @timeit "Read HITRAN" lines = AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(artifact(mol_name)); FT)
                @debug "Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)"
                absorption_model = AtmosphericAbsorption.LineByLineModel(lines;
                    profile = ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    cpf = ap.CEF,
                    architecture = _to_aa_arch(params.architecture),
                    vmr = 0)
                @timeit "Absorption Coeff" compute_absorption_profile!(τ_abs[i_band], absorption_model, params.spec_bands[i_band], profile.vmr[mol_name], profile)
            else
                # TIMED: this branch used to be the only unmeasured stage in the
                # model build, which is how a 19x regression in it went unnoticed
                # (a verbose run reported "Tot / % measured: 11.3%", with the whole
                # LUT absorption cost sitting in the unmeasured remainder).
                @timeit "Absorption Coeff LUT" compute_absorption_profile!(τ_abs[i_band], ap.luts[i_band][molec_i], params.spec_bands[i_band], profile.vmr[mol_name], profile)
            end
        end

        # H₂O line absorption (driven by q). Use the band's H2O LUT if the
        # parser found one inside LUTfiles; otherwise fall back to artifact.
        if any(!iszero, params.q)
            if ap.h2o_lut[i_band] !== nothing
                @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                    τ_abs[i_band], ap.h2o_lut[i_band],
                    params.spec_bands[i_band], profile)
            else
                @timeit "Read HITRAN H2O" lines_h2o = AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(artifact("H2O")); FT)
                @debug "Computing profile for H2O (q-driven) for band #$(i_band)"
                h2o_model = AtmosphericAbsorption.LineByLineModel(lines_h2o;
                    profile = ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    cpf = ap.CEF,
                    architecture = _to_aa_arch(params.architecture),
                    vmr = 0)
                @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                    τ_abs[i_band], h2o_model, params.spec_bands[i_band], profile)
            end
        end

        # Collision-induced absorption (HITRAN .cia files), if any.
        for (cia_i, cia_path) in enumerate(ap.cia_files)
            @timeit "CIA $(basename(cia_path))" begin
                cia_table = _load_configured_cia_table(
                    ap, cia_i, params.spec_bands[i_band], FT_vrs)
                Absorption.compute_τ_cia!(τ_abs[i_band], cia_table, profile,
                                           profile.vmr)
            end
        end

        # MT_CKD H₂O continuum (self + foreign), if a reference table is configured.
        if !isempty(ap.mtckd_file)
            @timeit "MT_CKD H2O continuum" begin
                mtckd_table = Absorption.load_mtckd(ap.mtckd_file)
                Absorption.compute_τ_h2o_continuum!(τ_abs[i_band], mtckd_table,
                                                     params.spec_bands[i_band], profile,
                                                     profile.vmr_h2o)
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];

    # τ_aer[iBand][iAer, nSpec, iZ] — 3-D matching the elastic forward path
    τ_aer = [zeros(FT_vrs, n_aer, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # k_ref uses the reference refractive index n_ref (normalisation convention),
        # NOT the aerosol's own nᵣ/nᵢ. Build a separate ref aerosol so that
        # mie_aerosol is never mutated — the in-band Mie below uses the aerosol's
        # true nᵣ/nᵢ, not n_ref.
        ref_aerosol_vs = Aerosol(size_distribution,
                                 real(params.scattering_params.n_ref),
                                 -imag(params.scattering_params.n_ref))
        mie_model_ref_vs = make_mie_model(params.scattering_params.decomp_type,
                                          ref_aerosol_vs,
                                          params.scattering_params.λ_ref,
                                          params.polarization_type,
                                          truncation_type,
                                          params.scattering_params.r_max,
                                          params.scattering_params.nquad_radius)
        k_ref          = compute_ref_aerosol_extinction(mie_model_ref_vs, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols (architecture from params selects CPU vs GPU Mie):
            mie_model      = make_mie_model(params.scattering_params.decomp_type,
                                            mie_aerosol,
                                            (maximum(curr_band_λ)+minimum(curr_band_λ))/2,
                                            params.polarization_type,
                                            truncation_type,
                                            params.scattering_params.r_max,
                                            params.scattering_params.nquad_radius;
                                            architecture = params.architecture)

            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually).
            # Single-verb call: dispatches CPU/GPU off mie_model.architecture.
            # Output type defaults to the MieModel's FT parameter (== params.float_type == FT_vrs),
            # so a Float32 model produces Float32 greek coefficients with no explicit override.
            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model);

            # Compute truncated aerosol optical properties (phase function and fᵗ).
            # Safety guard: only run δBGE forward-peak truncation when the raw
            # Greek series exceeds the projection cap. See the matching guard in
            # the forward-model `model_from_parameters` aerosol loop.
            β_len = length(aerosol_optics_raw.greek_coefs.β)
            aerosol_optics[i_band][i_aer] =
                if truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                    Scattering.truncate_phase(truncation_type,
                                              aerosol_optics_raw; reportFit=false)
                else
                    Scattering.truncate_phase(Scattering.NoTruncation(),
                                              aerosol_optics_raw)
                end

            # Track greek coef length for l_max computation (was previously
            # left at zero, which made `l_max[i_band]` always 0 in the
            # presence of aerosols and silently dropped Mie Fourier moments).
            l_max_aer[i_aer, i_band] =
                truncation_type isa Scattering.δBGE ?
                    min(length(aerosol_optics[i_band][i_aer].greek_coefs.β), truncation_type.l_max) :
                    length(aerosol_optics[i_band][i_aer].greek_coefs.β)

            # Compute nAer aerosol optical thickness profiles. RT_Aerosol stores
            # the vertical pressure distribution as `c_aero.profile` (a
            # `Distributions.Distribution`); use the 3-arg `getAerosolLayerOptProp`
            # that consumes it directly, matching the elastic site.
            # τ_aer is 3-D [iAer, nSpec, iLayer]; k is a scalar at band-center.
            τ_profile_vs = CoreRT.getAerosolLayerOptProp(1, c_aero.profile, profile)
            τ_aer[i_band][i_aer, :, :] =
                params.scattering_params.rt_aerosols[i_aer].τ_ref *
                (aerosol_optics[i_band][i_aer].k/k_ref) .*
                ones(FT_vrs, size(τ_aer[i_band], 2)) .* τ_profile_vs'
        end
    end

    # Compute per-band l_max from aerosol greek coefficient lengths.
    l_max = zeros(Int, n_bands)
    for i_band = 1:n_bands
        if n_aer > 0
            l_max[i_band] = maximum(l_max_aer[:, i_band])
        else
            l_max[i_band] = params.l_trunc
        end
    end
    # Per-band Fourier loop bound (order). See Phase C note in the matching
    # `model_from_parameters` site above.
    # VS-plus expands the user's single spectral band into `n_bands` sub-bands
    # (incident + vibrational ±). The user's YAML normally has one `surface:`
    # entry; replicate it so per-band lookups in `_band_components` and the
    # RTModel constructor have one BRDF per sub-band.
    if length(params.brdf) == 1 && n_bands > 1
        params.brdf = [params.brdf[1] for _ in 1:n_bands]
    end
    components_per_band = [_band_components(params, aerosol_optics, sources, i_band)
                            for i_band in 1:n_bands]
    m_max_bands = _derive_m_max_bands_via_traits(l_max, params.max_m,
                                                  components_per_band,
                                                  quad_points.Nstreams;
                                                  greek_beta_cutoff=params.greek_beta_cutoff)
    n_fourier_moments_bands = m_max_bands .+ 1

    # Build the hierarchical RTModel
    FT_vrs2 = params.float_type
    solver = SolverConfig{FT_vrs2, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        m_max_bands,
        n_fourier_moments_bands,
        l_max,
        params.l_trunc,
        FT_vrs2(params.Δ_angle),
        FT_vrs2(params.depol),
        true,   # use_component_traits — flipped on in Phase C
    )
    spec_bands_ft2 = [convert(Vector{FT_vrs2}, b) for b in params.spec_bands]
    atm = Atmosphere(profile, spec_bands_ft2)
    rayleigh_s = RayleighScattering(greek_rayleigh, greek_cabannes, FT_vrs2.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    optics = Optics(rayleigh_s, aerosols_s, τ_abs, τ_rayl)
    numerics = _convert_numerics(params.numerics, FT_vrs2)
    return RTModel(params.architecture, solver, numerics, obs_geom, quad_points, atm, optics, params.brdf, sources)
end

function loadAbsco(file; scale=(1.0))
    absco = Dataset(file)
    mol = absco["Gas_Index"][1]
    
    cs_name = "Gas_"* mol * "_Absorption"
    # Loading cross sections:
    σ = Float32(scale)*absco[cs_name][:]
    # Temperature
    T = absco["Temperature"][:]
    p = absco["Pressure"][:]/100
    ν = absco["Wavenumber"][:]
    return Absorption.AbscoTable(parse(Int,mol), -1, ν, σ, p, T )
end
"Linear interpolation weights on the native wavenumber coordinate."
@inline _spectral_weight(ν, ν₀, ν₁) = (ν - ν₀) / (ν₁ - ν₀)

"Interpolate equal-shaped node arrays along a new final spectral dimension."
function _interpolate_phase_nodes(ν_spec, ν_nodes, values)
    p = sortperm(ν_nodes)
    x = collect(ν_nodes[p])
    y = values[p]
    length(x) in (2, 3) || throw(ArgumentError("phase interpolation requires 2 or 3 nodes"))
    samples = if length(x) == 2
        [begin
            w = _spectral_weight(ν, x[1], x[2])
            (one(w) - w) .* y[1] .+ w .* y[2]
         end for ν in ν_spec]
    else
        [_natural_cubic_three((ν,), x, y)[1] for ν in ν_spec]
    end
    # Write columns into a preallocated dense array instead of
    # `cat(samples...; dims=...)`: splatting 12301 arrays into one cat call is
    # quadratic-ish in dispatch and shape recursion (and routes through
    # SparseArrays' _cat), and profiled as ~80% of an entire model build —
    # the dominant term of the 19x post-merge build regression. The per-point
    # sample values above are untouched, so results are bit-identical.
    proto = first(samples)
    out = Array{eltype(proto)}(undef, size(proto)..., length(samples))
    for (i, sample) in enumerate(samples)
        selectdim(out, ndims(out), i) .= sample
    end
    return out
end

function _pad_moment(v::AbstractVector, n::Int)
    out = zeros(eltype(v), n)
    copyto!(out, 1, v, 1, length(v))
    return out
end

function _pad_moment(v::AbstractMatrix, n::Int)
    out = zeros(eltype(v), size(v, 1), n)
    @views out[:, 1:size(v, 2)] .= v
    return out
end

"""Interpolate independently truncated endpoint optics onto every ν sample."""
function _spectralize_truncated_endpoints(a₀::AerosolOptics,
                                          a₁::AerosolOptics,
                                          ν_spec;
                                          reference::Union{Nothing,AerosolOptics}=nothing,
                                          ν_ref=nothing)
    FT = promote_type(eltype(ν_spec), typeof(a₀.k), typeof(a₁.k))
    ν₀, ν₁ = FT(first(ν_spec)), FT(last(ν_spec))
    w = FT[_spectral_weight(ν, ν₀, ν₁) for ν in ν_spec]
    n_spec = length(w)
    n_l = maximum(length(x.greek_coefs.β) for x in
                  (reference === nothing ? (a₀, a₁) : (a₀, reference, a₁)))
    fields = Dict{Symbol,Matrix{FT}}()
    for fn in (:α, :β, :γ, :δ, :ϵ, :ζ)
        g₀ = FT.(_pad_moment(getfield(a₀.greek_coefs, fn), n_l))
        g₁ = FT.(_pad_moment(getfield(a₁.greek_coefs, fn), n_l))
        fields[fn] = reference === nothing ?
            g₀ .* reshape(one(FT) .- w, 1, n_spec) .+
            g₁ .* reshape(w, 1, n_spec) :
            _interpolate_phase_nodes(ν_spec, [ν₀, ν_ref, ν₁],
                [g₀, FT.(_pad_moment(getfield(reference.greek_coefs, fn), n_l)), g₁])
    end
    gc = GreekCoefs(fields[:α], fields[:β], fields[:γ],
                    fields[:δ], fields[:ϵ], fields[:ζ])
    k = reference === nothing ? FT.(a₀.k .* (one(FT) .- w) .+ a₁.k .* w) :
        vec(_interpolate_phase_nodes(ν_spec, [ν₀, ν_ref, ν₁],
                                     [fill(FT(a₀.k),1), fill(FT(reference.k),1), fill(FT(a₁.k),1)]))
    ksca₀, ksca₁ = a₀.k * a₀.ω̃, a₁.k * a₁.ω̃
    ksca = reference === nothing ? FT.(ksca₀ .* (one(FT) .- w) .+ ksca₁ .* w) :
        vec(_interpolate_phase_nodes(ν_spec, [ν₀, ν_ref, ν₁],
            [fill(FT(ksca₀),1), fill(FT(reference.k * reference.ω̃),1), fill(FT(ksca₁),1)]))
    ω̃ = ksca ./ k
    fᵗ = reference === nothing ? FT.(a₀.fᵗ .* (one(FT) .- w) .+ a₁.fᵗ .* w) :
        vec(_interpolate_phase_nodes(ν_spec, [ν₀, ν_ref, ν₁],
                                     [fill(FT(a₀.fᵗ),1), fill(FT(reference.fᵗ),1), fill(FT(a₁.fᵗ),1)]))
    node_ν = reference === nothing ? FT[ν₀, ν₁] : FT[ν₀, ν_ref, ν₁]
    node_greek = reference === nothing ? [a₀.greek_coefs, a₁.greek_coefs] :
                                           [a₀.greek_coefs, reference.greek_coefs, a₁.greek_coefs]
    # Untruncated node Greeks (for the TMS exact-SS evaluation), carried
    # through from truncate_phase; nothing when any endpoint was not
    # actually truncated (exact == truncated there).
    node_endpoints = reference === nothing ? (a₀, a₁) : (a₀, reference, a₁)
    _gc_ft(g) = GreekCoefs(FT.(g.α), FT.(g.β), FT.(g.γ),
                           FT.(g.δ), FT.(g.ϵ), FT.(g.ζ))
    node_raw = all(e -> e.phase_greek_raw !== nothing, node_endpoints) ?
        [_gc_ft(e.phase_greek_raw[1]) for e in node_endpoints] : nothing
    return AerosolOptics(greek_coefs=gc, ω̃=ω̃, k=k, fᵗ=fᵗ,
                          phase_ν=node_ν, phase_greek=node_greek,
                          phase_greek_raw=node_raw)
end

function _spectralize_truncated_endpoints(a₀::AerosolOptics,
                                          l₀::linAerosolOptics,
                                          a₁::AerosolOptics,
                                          l₁::linAerosolOptics,
                                          ν_spec;
                                          reference::Union{Nothing,AerosolOptics}=nothing,
                                          lin_reference::Union{Nothing,linAerosolOptics}=nothing,
                                          ν_ref=nothing)
    a = _spectralize_truncated_endpoints(a₀, a₁, ν_spec;
                                         reference, ν_ref)
    FT = eltype(a.k)
    w = FT[_spectral_weight(ν, first(ν_spec), last(ν_spec)) for ν in ν_spec]
    n_spec = length(w)
    n_param = size(l₀.lin_greek_coefs.β̇, 1)
    n_l = size(a.greek_coefs.β, 1)
    fields = Dict{Symbol,Array{FT,3}}()
    for lfn in (:α̇, :β̇, :γ̇, :δ̇, :ϵ̇, :ζ̇)
        g₀ = FT.(_pad_moment(getfield(l₀.lin_greek_coefs, lfn), n_l))
        g₁ = FT.(_pad_moment(getfield(l₁.lin_greek_coefs, lfn), n_l))
        fields[lfn] = lin_reference === nothing ?
            reshape(g₀, n_param, n_l, 1) .* reshape(one(FT) .- w, 1, 1, n_spec) .+
            reshape(g₁, n_param, n_l, 1) .* reshape(w, 1, 1, n_spec) :
            _interpolate_phase_nodes(ν_spec, a.phase_ν,
                [g₀, FT.(_pad_moment(getfield(lin_reference.lin_greek_coefs, lfn), n_l)), g₁])
    end
    lgc = linGreekCoefs(fields[:α̇], fields[:β̇], fields[:γ̇],
                        fields[:δ̇], fields[:ϵ̇], fields[:ζ̇])
    # `_interpolate_phase_nodes` keeps the (n_param, 1) node shape and returns
    # (n_param, 1, n_spec); collapse the singleton node axis so both branches
    # hand downstream code the same (n_param, n_spec) layout.
    _drop_node_axis(x) = dropdims(x; dims=2)
    k̇₀ = reshape(l₀.k̇, n_param, 1)
    k̇₁ = reshape(l₁.k̇, n_param, 1)
    ω̃̇₀ = reshape(l₀.ω̃̇, n_param, 1)
    ω̃̇₁ = reshape(l₁.ω̃̇, n_param, 1)
    k̇ = lin_reference === nothing ?
        k̇₀ .* reshape(one(FT) .- w, 1, n_spec) .+ k̇₁ .* reshape(w, 1, n_spec) :
        _drop_node_axis(_interpolate_phase_nodes(ν_spec, a.phase_ν,
            [k̇₀, reshape(lin_reference.k̇, n_param, 1), k̇₁]))
    kscȧ₀ = k̇₀ .* a₀.ω̃ .+ a₀.k .* ω̃̇₀
    kscȧ₁ = k̇₁ .* a₁.ω̃ .+ a₁.k .* ω̃̇₁
    kscȧ = lin_reference === nothing ?
        kscȧ₀ .* reshape(one(FT) .- w, 1, n_spec) .+ kscȧ₁ .* reshape(w, 1, n_spec) :
        _drop_node_axis(_interpolate_phase_nodes(ν_spec, a.phase_ν,
            [kscȧ₀,
             reshape(lin_reference.k̇ .* reference.ω̃ .+
                     reference.k .* lin_reference.ω̃̇, n_param, 1), kscȧ₁]))
    ω̃̇ = (kscȧ .- reshape(a.ω̃, 1, n_spec) .* k̇) ./
          reshape(a.k, 1, n_spec)
    ḟᵗ = lin_reference === nothing ?
        reshape(l₀.ḟᵗ, n_param, 1) .* reshape(one(FT) .- w, 1, n_spec) .+
        reshape(l₁.ḟᵗ, n_param, 1) .* reshape(w, 1, n_spec) :
        _drop_node_axis(_interpolate_phase_nodes(ν_spec, a.phase_ν,
            [reshape(l₀.ḟᵗ,n_param,1), reshape(lin_reference.ḟᵗ,n_param,1),
             reshape(l₁.ḟᵗ,n_param,1)]))
    node_lin = lin_reference === nothing ?
        [l₀.lin_greek_coefs, l₁.lin_greek_coefs] :
        [l₀.lin_greek_coefs, lin_reference.lin_greek_coefs, l₁.lin_greek_coefs]
    return a, linAerosolOptics(lin_greek_coefs=lgc, ω̃̇=ω̃̇, k̇=k̇, ḟᵗ=ḟᵗ,
                                phase_ν=copy(a.phase_ν), phase_lin_greek=node_lin)
end

"Natural cubic spline through exactly three ordered knots."
function _natural_cubic_three(x, knots, values)
    p = sortperm(knots)
    x₀, x₁, x₂ = knots[p]
    y₀, y₁, y₂ = values[p]
    h₀, h₁ = x₁ - x₀, x₂ - x₁
    M₁ = 3 * ((y₂ - y₁) / h₁ - (y₁ - y₀) / h₀) / (h₀ + h₁)
    function segment(xv, xa, xb, ya, yb, Ma, Mb)
        h = xb - xa
        return Ma * (xb - xv)^3 / (6h) + Mb * (xv - xa)^3 / (6h) +
               (ya - Ma * h^2 / 6) * (xb - xv) / h +
               (yb - Mb * h^2 / 6) * (xv - xa) / h
    end
    return [xv <= x₁ ? segment(xv, x₀, x₁, y₀, y₁, zero(M₁), M₁) :
                       segment(xv, x₁, x₂, y₁, y₂, M₁, zero(M₁)) for xv in x]
end

"""
    _aerosol_τ_slice(τ_ref, k, k_ref, τ_profile, ν_spec, ν_ref)

The (nSpec × nLayers) aerosol optical-depth slice
`τ_ref · aod_scale(λ) · τ_profile'`, shared by the fresh model build and the
`update_model!`/`update_aerosol_loading!` fast paths. A single definition is
what makes scene updates bit-exact against a freshly built model — do not
inline this product elsewhere. `k` may be a scalar (single-λ band) or an
nSpec vector; the spectral scale reduces to `k/k_ref` unless λ_ref falls
inside the band (three-node natural cubic, see `_aod_spectral_scale`).
"""
_aerosol_τ_slice(τ_ref, k, k_ref, τ_profile, ν_spec, ν_ref) =
    τ_ref .* _aod_spectral_scale(ν_spec, k, k_ref, ν_ref) .* τ_profile'

"AOD spectral scale, with exact unity at an in-band reference wavenumber."
function _aod_spectral_scale(ν_spec, k_spec, k_ref, ν_ref)
    νlo, νhi = extrema(ν_spec)
    if νlo < ν_ref < νhi
        return _natural_cubic_three(ν_spec,
            [first(ν_spec), ν_ref, last(ν_spec)],
            [first(k_spec) / k_ref, one(eltype(k_spec)), last(k_spec) / k_ref])
    end
    return k_spec ./ k_ref
end
