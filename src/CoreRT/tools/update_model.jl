#=

BatchContext and update_model! — batch-processing facility for vSmartMOM.

Phase 1: Only p/T profiles, specific humidity q, and trace-gas VMRs change
per scene; spectral bands, quadrature, geometry, surfaces, aerosols, and
polarization remain fixed. BatchContext caches everything that is expensive
and scene-invariant so that per-scene updates skip HITRAN re-parsing and
Mie re-computation entirely.

Phase 2 (future): aerosol loading updates will also benefit from k_ref
caching; it is stored here now.

=#

"""
    BatchContext

Pre-built context for efficient batch radiative-transfer over many atmospheric
scenes that share the same spectral configuration, geometry, surface BRDF,
aerosol microphysics, and quadrature setup but differ in temperature/pressure
profiles, humidity, and/or trace-gas VMRs.

Construct once from a `vSmartMOM_Parameters` object; then call
[`update_model!`](@ref) before each `rt_run` invocation.

# Example usage

```julia
params = parameters_from_yaml("my_config.yaml")
ctx    = BatchContext(params)

for scene in scenes
    update_model!(ctx;
        T      = scene.T,        # length-Nz temperature profile [K]
        p_half = scene.p_half,   # length-(Nz+1) half-level pressures [hPa]
        q      = scene.q,        # length-Nz specific humidity [kg/kg]
        vmr    = scene.vmr)      # Dict{String,Any} — same keys as params.absorption_params.vmr
    R, T = rt_run(ctx.model)
    # ... process R, T ...
end
```

# What may change per scene

- Temperature profile `T` [K]
- Pressure half-levels `p_half` [hPa]
- Specific humidity `q` [kg/kg]
- Trace-gas VMRs via the `vmr` keyword dict (same keys as the original
  `params.absorption_params.vmr`; adding new species requires a full rebuild)

# What requires a full rebuild (new `BatchContext`)

- Spectral bands (`params.spec_bands`)
- Geometry (SZA, VZA, VAZ, observer altitude)
- Surface BRDF
- Aerosol microphysics (size distribution, refractive index, τ_ref)
- Polarization type, quadrature type, truncation

# Thread safety

One `BatchContext` per thread. `update_model!` and `rt_run(ctx.model)` mutate
shared arrays and must **not** be called concurrently on the same `ctx`. If you
run parallel scenes, create one `BatchContext` per worker thread (each
`BatchContext` owns its own `RTModel` with independent mutable arrays).

# Fields

- `model`: the [`RTModel`](@ref) that is updated in place by `update_model!`
- `params`: the original `vSmartMOM_Parameters` used for construction (kept for
  meta-information such as molecule lists, YAML knobs, float type)
- `absorption_models`: cached `AtmosphericAbsorption.LineByLineModel` objects,
  one per `(band, species)` pair — eliminates HITRAN re-parsing per scene
- `h2o_models`: cached H₂O `LineByLineModel` per band (`nothing` when a LUT is
  used or when `q` is all-zero)
- `k_ref`: reference-wavelength aerosol extinction coefficients, one per aerosol
  species (needed for Phase 2 aerosol-loading updates)
- `n_bands`, `n_aerosols`, `Nz`: scene-invariant dimension bookmarks
- `profile_reduction_n`: reduction target passed to `reduce_profile`; `-1` means
  no reduction
"""
mutable struct BatchContext
    "The RTModel updated in place by update_model!"
    model::RTModel
    "Original parameters used to build this context"
    params::vSmartMOM_Parameters
    "Cached LineByLineModel per (band, species index in all_species list)"
    absorption_models::Vector{Vector{Any}}      # [i_band][molec_i]
    "Cached H₂O LineByLineModel per band (nothing when LUT / q≡0)"
    h2o_models::Vector{Any}                     # [i_band]
    "Reference-wavelength extinction coefficient per aerosol species (Phase 2)"
    k_ref::Vector{Float64}
    "Number of spectral bands (scene-invariant)"
    n_bands::Int
    "Number of aerosol species (scene-invariant)"
    n_aerosols::Int
    "Number of atmospheric layers after profile reduction (scene-invariant)"
    Nz::Int
    "Profile reduction target (-1 = no reduction)"
    profile_reduction_n::Int
end

"""
    BatchContext(params::vSmartMOM_Parameters; sources::AbstractSource = SolarBeam()) -> BatchContext

Build a `BatchContext` from the given parameters.

This is the **expensive** constructor: it calls `model_from_parameters(params)`
(which runs Mie, builds the aerosol optics, and reads HITRAN files), then
caches the `AtmosphericAbsorption.LineByLineModel` objects so that subsequent
[`update_model!`](@ref) calls can skip HITRAN re-parsing entirely.

Keyword argument `sources` is forwarded to `model_from_parameters`.
"""
function BatchContext(params::vSmartMOM_Parameters;
                      sources::AbstractSource = SolarBeam())
    FT = params.float_type

    # ── 1. Build the full RTModel (expensive; done once) ───────────────────
    model = model_from_parameters(params; sources)

    n_bands    = length(params.spec_bands)
    n_aerosols = isnothing(params.scattering_params) ? 0 :
                 length(params.scattering_params.rt_aerosols)
    Nz         = length(model.atmosphere.profile.p_full)
    profile_reduction_n = params.profile_reduction_n

    # ── 2. Cache per-(band, species) LineByLineModel objects ───────────────
    # Mirror exactly what model_from_parameters does in the absorption block,
    # but store the model objects instead of discarding them after use.
    absorption_models = Vector{Vector{Any}}(undef, n_bands)
    h2o_models        = Vector{Any}(undef, n_bands)

    for i_band in 1:n_bands
        absorption_models[i_band] = Any[]
        h2o_models[i_band]        = nothing

        isnothing(params.absorption_params) && continue
        ap = params.absorption_params

        all_species = vcat(ap.fixed_molecules[i_band], ap.variable_molecules[i_band])

        for (molec_i, mol_name) in enumerate(all_species)
            if isempty(ap.luts)
                lines = AtmosphericAbsorption.load_lines(
                    AtmosphericAbsorption.HitranPort(artifact(mol_name)); FT)
                absorption_model = AtmosphericAbsorption.LineByLineModel(lines;
                    profile      = ap.broadening_function,
                    wing_cutoff  = ap.wing_cutoff,
                    cpf          = ap.CEF,
                    architecture = _to_aa_arch(params.architecture),
                    vmr          = 0)
                push!(absorption_models[i_band], absorption_model)
            else
                # LUT — store the LUT object itself so _compute_band_absorption!
                # can dispatch through the same compute_absorption_profile! path.
                push!(absorption_models[i_band], ap.luts[i_band][molec_i])
            end
        end

        # H₂O model (if needed)
        if any(!iszero, params.q)
            if ap.h2o_lut[i_band] !== nothing
                h2o_models[i_band] = ap.h2o_lut[i_band]
            else
                lines_h2o = AtmosphericAbsorption.load_lines(
                    AtmosphericAbsorption.HitranPort(artifact("H2O")); FT)
                h2o_models[i_band] = AtmosphericAbsorption.LineByLineModel(lines_h2o;
                    profile      = ap.broadening_function,
                    wing_cutoff  = ap.wing_cutoff,
                    cpf          = ap.CEF,
                    architecture = _to_aa_arch(params.architecture),
                    vmr          = 0)
            end
        end
    end

    # ── 3. Cache k_ref per aerosol (Mie reference extinction) ─────────────
    # model_from_parameters discards k_ref after building τ_aer; we keep it
    # for Phase 2 aerosol-loading updates (changing τ_ref or profile without
    # rebuilding the full Mie optics).
    k_ref_vec = Float64[]
    if !isnothing(params.scattering_params)
        truncation_type = _resolved_truncation(params, FT)
        for i_aer in 1:n_aerosols
            c_aero = params.scattering_params.rt_aerosols[i_aer]
            if _has_analytic_phase_function(c_aero)
                # Analytic phase function: k is per-aerosol constant = 1 because
                # extinction_cross_section is set to 1 in _analytic_aerosol_optics.
                push!(k_ref_vec, 1.0)
            else
                curr_aerosol   = c_aero.aerosol
                mie_aerosol    = Aerosol(curr_aerosol.size_distribution,
                                         curr_aerosol.nᵣ, curr_aerosol.nᵢ)
                mie_model_ref  = make_mie_model(params.scattering_params.decomp_type,
                                                mie_aerosol,
                                                params.scattering_params.λ_ref,
                                                params.polarization_type,
                                                truncation_type,
                                                params.scattering_params.r_max,
                                                params.scattering_params.nquad_radius)
                mie_model_ref.aerosol.nᵣ = real(params.scattering_params.n_ref)
                mie_model_ref.aerosol.nᵢ = -imag(params.scattering_params.n_ref)
                push!(k_ref_vec, Float64(compute_ref_aerosol_extinction(mie_model_ref, FT)))
            end
        end
    end

    return BatchContext(model, params, absorption_models, h2o_models,
                        k_ref_vec, n_bands, n_aerosols, Nz, profile_reduction_n)
end

# ============================================================================
# Internal shared accumulation helper
# ============================================================================

"""
    _compute_band_absorption!(τ_abs_band, absorption_models_band, h2o_model,
                              spec_band, all_species, profile, q_nonzero)

Accumulate gas absorption optical depths into `τ_abs_band` for one spectral band.

**IMPORTANT**: `τ_abs_band` must be zeroed by the caller before this function
is called. `compute_absorption_profile!` uses `+=`; zeroing first is essential
to prevent stale values from previous scenes from contaminating the new result.

This helper is shared between `update_model!` and (optionally) any future
wrapper that mirrors the `model_from_parameters` accumulation logic without
duplicating it. Using this helper keeps the gas-accumulation path in one place
so it cannot silently diverge across call sites.

# Arguments
- `τ_abs_band`: pre-zeroed `(nSpec × Nz)` matrix to accumulate into
- `absorption_models_band`: vector of cached model objects (one per non-H₂O
  species, in the same order as `all_species`)
- `h2o_model`: cached H₂O model (or `nothing` when not needed)
- `spec_band`: wavenumber grid for this band [cm⁻¹]
- `all_species`: `Vector{String}` — ordered species names matching
  `absorption_models_band`
- `profile`: `AtmosphericProfile` (already reduced if needed)
- `q_nonzero`: `Bool` — whether the humidity profile has any non-zero values
"""
function _compute_band_absorption!(τ_abs_band, absorption_models_band, h2o_model,
                                   spec_band, all_species, profile, q_nonzero)
    for (molec_i, mol_name) in enumerate(all_species)
        model_obj = absorption_models_band[molec_i]
        vmr_curr  = profile.vmr[mol_name]
        compute_absorption_profile!(τ_abs_band, model_obj, spec_band, vmr_curr, profile)
    end

    # H₂O line absorption (driven by q / vmr_h2o)
    if q_nonzero && h2o_model !== nothing
        compute_absorption_profile!(τ_abs_band, h2o_model, spec_band,
                                    profile.vmr_h2o, profile)
    end
end

# ============================================================================
# Public API: update_model!
# ============================================================================

"""
    update_model!(ctx::BatchContext;
                  T      = nothing,
                  p_half = nothing,
                  q      = nothing,
                  vmr    = nothing)

Update `ctx.model` in place for a new atmospheric scene.

Only the arguments you supply are updated; `nothing` means "keep the current
value from the model". All updated arguments are type-converted to the model's
float type (`params.float_type`) before use.

After this call `rt_run(ctx.model)` will produce radiances for the new scene.

# Keyword arguments

- `T::AbstractVector`: Temperature profile [K], length `Nz` (full levels,
  TOA to BOA). If `profile_reduction_n != -1`, pass the **unreduced** profile
  length (same as the original `params.T`).
- `p_half::AbstractVector`: Half-level pressures [hPa], length `Nz + 1` (or
  `N_orig + 1` before reduction). Must be strictly positive and monotonically
  increasing (TOA to BOA, i.e. `p_half[end]` is surface pressure).
- `q::AbstractVector`: Specific humidity [kg/kg], length `Nz` (or `N_orig`).
- `vmr::Dict{String,Any}`: Trace-gas volume mixing ratios. **Keys must be a
  subset of the configured molecules** — supplying an unknown species key raises
  an error (add it to `params.absorption_params.vmr` and rebuild the context
  instead). Values may be scalars or length-`Nz` vectors.

# Thread safety

`update_model!` and `rt_run(ctx.model)` share mutable arrays. Do not call them
concurrently on the same `ctx`. Use one `BatchContext` per thread for parallelism.

# Example

```julia
params = parameters_from_yaml("config.yaml")
ctx    = BatchContext(params)

for scene in scenes
    update_model!(ctx;
        T      = scene.T,
        p_half = scene.p_half,
        q      = scene.q,
        vmr    = Dict("O2" => 0.21, "CH4" => scene.ch4_vmr))
    R, T_trans = rt_run(ctx.model)
end
```
"""
function update_model!(ctx::BatchContext;
                       T      = nothing,
                       p_half = nothing,
                       q      = nothing,
                       vmr    = nothing)
    model  = ctx.model
    params = ctx.params
    FT     = params.float_type
    ap     = params.absorption_params   # may be nothing
    # profile_cur is the currently-stored (already-reduced) profile; used only
    # for in-place copying after the new profile is computed.
    profile_cur = model.atmosphere.profile

    # ── 0. Validate and merge inputs ────────────────────────────────────────
    # The public API always accepts **unreduced** profile inputs (i.e. the
    # same length as params.T / params.p / params.q). This is necessary for
    # partial updates: if the user supplies only T (no p_half), we must fall
    # back to the original params.p (unreduced) so that the validate-then-
    # reduce path is consistent. Falling back to the stored (already-reduced)
    # profile_cur.p_half would produce a length mismatch when T was passed
    # as an unreduced vector.
    #
    # Rule: if a field is `nothing`, use the original params value (unreduced).
    T_new      = T      === nothing ? convert(Vector{FT}, params.T) : convert(Vector{FT}, T)
    q_new      = q      === nothing ? convert(Vector{FT}, params.q) : convert(Vector{FT}, q)
    p_half_new = p_half === nothing ? convert(Vector{FT}, params.p) : convert(Vector{FT}, p_half)

    # Validate p_half length: it must be Nz+1 (post-reduction) or N_orig
    # (params.p is already the half-level array; it has length Nz_orig + 1,
    #  which equals length(params.T) + 1 — the same thing as length(params.p)).
    n_half_ok = (length(p_half_new) == ctx.Nz + 1) ||
                (ctx.profile_reduction_n != -1 && length(p_half_new) == length(params.p))
    n_half_ok || error(
        "update_model!: p_half length $(length(p_half_new)) does not match " *
        "expected Nz+1 = $(ctx.Nz + 1)" *
        (ctx.profile_reduction_n != -1 ?
            " or unreduced length $(length(params.p))" : ""))

    # T and q must match p_half (they are full-level quantities: length(p_half) - 1)
    n_T_expect = length(p_half_new) - 1
    length(T_new) == n_T_expect || error(
        "update_model!: T length $(length(T_new)) must equal length(p_half)-1 = $n_T_expect")
    length(q_new) == n_T_expect || error(
        "update_model!: q length $(length(q_new)) must equal length(p_half)-1 = $n_T_expect")

    # Check all T and p values are finite and positive
    all(isfinite, T_new)      || error("update_model!: T contains non-finite values")
    all(>(0), T_new)          || error("update_model!: T must be strictly positive")
    all(isfinite, p_half_new) || error("update_model!: p_half contains non-finite values")
    all(>(0), p_half_new)     || error("update_model!: p_half must be strictly positive")
    all(isfinite, q_new)      || error("update_model!: q contains non-finite values")

    # VMR validation: new keys must be subset of configured molecules
    if vmr !== nothing && !isnothing(ap)
        configured_keys = keys(ap.vmr)
        for k in keys(vmr)
            k in configured_keys || error(
                "update_model!: VMR key \"$k\" was not in the original " *
                "params.absorption_params.vmr (configured: $(collect(configured_keys))). " *
                "To add a new species, rebuild the BatchContext.")
        end
    end

    # ── 1. Recompute AtmosphericProfile fields ──────────────────────────────
    # Build a merged vmr dict: start from the current profile's vmr, apply overrides.
    vmr_merged = if vmr === nothing || isnothing(ap)
        ap === nothing ? Dict{String, Any}() : ap.vmr
    else
        # Overlay the supplied keys over the configured defaults
        d = Dict{String, Any}(ap.vmr)
        for (k, v) in vmr
            d[k] = v
        end
        d
    end

    # compute_atmos_profile_fields: same call as in model_from_parameters
    p_full_new, p_half_out, vmr_h2o_new, vcd_dry_new, vcd_h2o_new, new_vmr, Δz_new =
        compute_atmos_profile_fields(T_new, p_half_new, q_new, vmr_merged)

    new_profile_unreduced = AtmosphericProfile(
        T_new, p_full_new, q_new, p_half_out, vmr_h2o_new,
        vcd_dry_new, vcd_h2o_new, new_vmr, Δz_new)

    # Apply profile reduction if needed
    new_profile = if ctx.profile_reduction_n != -1
        reduce_profile(ctx.profile_reduction_n, new_profile_unreduced)
    else
        new_profile_unreduced
    end

    Nz_new = length(new_profile.p_full)
    Nz_new == ctx.Nz || error(
        "update_model!: profile reduction produced Nz=$Nz_new layers but " *
        "BatchContext was built with Nz=$(ctx.Nz). This should not happen; " *
        "check that p_half spans the same column.")

    # Copy new profile into the existing profile arrays in-place.
    # RTModel/Atmosphere structs are immutable but every leaf Array is mutable.
    cur = profile_cur
    cur.T       .= new_profile.T
    cur.p_full  .= new_profile.p_full
    cur.q       .= new_profile.q
    cur.p_half  .= new_profile.p_half
    cur.vmr_h2o .= new_profile.vmr_h2o
    cur.vcd_dry .= new_profile.vcd_dry
    cur.vcd_h2o .= new_profile.vcd_h2o
    cur.Δz      .= new_profile.Δz
    # VMR dict: scalar values or vectors are updated key-by-key
    for (k, v) in new_profile.vmr
        if haskey(cur.vmr, k)
            if cur.vmr[k] isa AbstractArray && v isa AbstractArray
                cur.vmr[k] .= v
            else
                # scalar or mixed: just replace (Dict is mutable)
                cur.vmr[k] = v
            end
        else
            cur.vmr[k] = v
        end
    end

    # ── 2. Recompute τ_rayl in place per band ──────────────────────────────
    # Use the same depol logic as model_from_parameters.
    for i_band in 1:ctx.n_bands
        curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])

        νₘ  = FT(0.5) * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ  = FT(1.0e7) / νₘ
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(FT(1.0e7) / λₘ, FT(300))
        _, _ = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
        γ_air_Ray, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ, _n2, _o2)
        depol_air_Ray = 2γ_air_Ray / (1 + γ_air_Ray)
        depol_use_Ray = params.depol < 0 ? FT(depol_air_Ray) : FT(params.depol)

        model.optics.τ_rayl[i_band] .= getRayleighLayerOptProp(
            new_profile.p_half[end],
            curr_band_λ,
            depol_use_Ray,
            new_profile.vcd_dry)
    end

    # ── 3. Recompute τ_aer in place per band per aerosol ───────────────────
    # Microphysics (aerosol_optics, k) are fixed; only the vertical
    # redistribution changes with the new p_full/Δz.
    if !isnothing(params.scattering_params)
        for i_aer in 1:ctx.n_aerosols
            c_aero    = params.scattering_params.rt_aerosols[i_aer]
            k_ref_aer = ctx.k_ref[i_aer]

            if _has_analytic_phase_function(c_aero)
                τ_profile = getAerosolLayerOptProp(one(FT), c_aero.profile, new_profile)
                for i_band in 1:ctx.n_bands
                    model.optics.aerosols.τ_aer[i_band][i_aer, :] .= c_aero.τ_ref * τ_profile
                end
            else
                for i_band in 1:ctx.n_bands
                    k_aer = model.optics.aerosols.aerosol_optics[i_band][i_aer].k
                    τ_profile = getAerosolLayerOptProp(one(FT), c_aero.profile, new_profile)
                    model.optics.aerosols.τ_aer[i_band][i_aer, :] .=
                        c_aero.τ_ref * (k_aer / k_ref_aer) * τ_profile
                end
            end
        end
    end

    # ── 4. Re-accumulate τ_abs per band ────────────────────────────────────
    # IMPORTANT: zero first, then accumulate with +=.
    # compute_absorption_profile! (in atmo_prof.jl) uses +=; zeroing before
    # accumulation is essential to prevent stale optical depths from a previous
    # scene from contaminating the current result.
    if !isnothing(ap)
        q_nonzero = any(!iszero, q_new)
        for i_band in 1:ctx.n_bands
            # Zero entire band τ_abs (gas lines + CIA + continuum all live here)
            fill!(model.optics.τ_abs[i_band], zero(FT))

            all_species = vcat(ap.fixed_molecules[i_band], ap.variable_molecules[i_band])

            # Accumulate line absorption via the shared helper.
            # The helper iterates absorption_models_band (cached LineByLine/LUT
            # objects) in the same order as all_species.
            _compute_band_absorption!(
                model.optics.τ_abs[i_band],
                ctx.absorption_models[i_band],
                ctx.h2o_models[i_band],
                params.spec_bands[i_band],
                all_species,
                new_profile,
                q_nonzero)

            # Collision-induced absorption (CIA tables — re-read from disk like
            # model_from_parameters; CIA tables are small and rarely configured).
            for cia_path in ap.cia_files
                cia_table = Absorption.load_cia_table(
                    cia_path, params.spec_bands[i_band]; FT)
                Absorption.compute_τ_cia!(
                    model.optics.τ_abs[i_band], cia_table, new_profile, ap.vmr)
            end

            # MT_CKD H₂O continuum
            if !isempty(ap.mtckd_file)
                mtckd_table = Absorption.load_mtckd(ap.mtckd_file)
                Absorption.compute_τ_h2o_continuum!(
                    model.optics.τ_abs[i_band], mtckd_table,
                    params.spec_bands[i_band], new_profile, new_profile.vmr_h2o)
            end
        end
    end

    return nothing
end

# Export hint (actual export list lives in CoreRT.jl)
# export BatchContext, update_model!
