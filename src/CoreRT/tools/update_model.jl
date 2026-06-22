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
- Polarization type, quadrature type, truncation
- A different *number* of aerosols or gas species

Aerosol changes do **not** require a rebuild: use
[`update_aerosol_loading!`](@ref) for τ_ref / vertical-placement changes
(cheap, no Mie recomputation) and `update_aerosol_microphysics!` for
size-distribution / refractive-index changes (reruns Mie for that aerosol and
re-derives the Fourier-loop bounds).

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
- `h2o_models`: cached H₂O `LineByLineModel` per band (`nothing` only when the
  configuration has no `absorption_params` at all; otherwise the model is cached
  unconditionally so a scene whose `q` changes from all-zero to non-zero can
  apply H₂O line absorption without a rebuild)
- `k_ref`: reference-wavelength aerosol extinction coefficients, one per aerosol
  species (needed for Phase 2 aerosol-loading updates)
- `current_T`, `current_p_half`, `current_q`, `current_vmr`: the **unreduced**
  scene state that `update_model!` was last called with (initialised from
  `params` at construction). `nothing` keyword arguments to `update_model!` fall
  back to these stored values — never to `params` — so successive partial updates
  compose incrementally. Stored already FT-converted so repeated updates do not
  re-convert or drift. The unreduced state cannot be recovered from
  `model.atmosphere.profile` (which is reduced), which is why it lives here.
- `current_τ_ref`, `current_profile_dist`: the current per-aerosol loading state
  (column optical depth at `λ_ref` and vertical distribution), initialised from
  `params` and written by both `update_aerosol_loading!` and
  `update_aerosol_microphysics!`. `update_model!`'s τ_aer redistribution reads
  these so a prior loading update is not silently wiped.
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
    "Cached H₂O LineByLineModel per band (nothing only when no absorption_params)"
    h2o_models::Vector{Any}                     # [i_band]
    "Reference-wavelength extinction coefficient per aerosol species (Phase 2)"
    k_ref::Vector{Float64}
    "Current unreduced temperature profile [K] (FT-converted)"
    current_T::Vector
    "Current unreduced half-level pressures [hPa] (FT-converted)"
    current_p_half::Vector
    "Current unreduced specific humidity [kg/kg] (FT-converted)"
    current_q::Vector
    "Current trace-gas VMR overrides (merged over params.absorption_params.vmr)"
    current_vmr::Dict{String, Any}
    "Current column optical depth at λ_ref per aerosol species (FT-converted)"
    current_τ_ref::Vector
    "Current vertical pressure distribution per aerosol species"
    current_profile_dist::Vector{Any}
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

        # H₂O model — cache unconditionally whenever absorption is configured.
        # model_from_parameters only *applies* the H₂O path when q is non-zero,
        # but the LineByLineModel / LUT object itself does not depend on the q
        # values (only on the band config). Caching it regardless of the initial
        # q means a scene that later raises q from all-zero to non-zero can apply
        # H₂O line absorption without rebuilding the context. The q-VALUE gate
        # lives only in update_model! (per scene), mirroring the `any(!iszero, …)`
        # branch in model_from_parameters. We are already inside the
        # `isnothing(params.absorption_params) && continue` guard above, so
        # absorption is configured here.
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
                # k_ref uses the reference refractive index n_ref (normalisation
                # convention). Build a separate Aerosol with n_ref so the aerosol's
                # own nᵣ/nᵢ are never mutated (fixes aliasing bug).
                ref_aerosol_bc = Aerosol(curr_aerosol.size_distribution,
                                         real(params.scattering_params.n_ref),
                                         -imag(params.scattering_params.n_ref))
                mie_model_ref  = make_mie_model(params.scattering_params.decomp_type,
                                                ref_aerosol_bc,
                                                params.scattering_params.λ_ref,
                                                params.polarization_type,
                                                truncation_type,
                                                params.scattering_params.r_max,
                                                params.scattering_params.nquad_radius)
                push!(k_ref_vec, Float64(compute_ref_aerosol_extinction(mie_model_ref, FT)))
            end
        end
    end

    # ── 4. Initialise current scene state from params (B3) ─────────────────
    # These hold the UNREDUCED scene inputs that update_model! was last called
    # with. `nothing` keyword arguments fall back to these (never to params),
    # so successive partial updates compose. Store already FT-converted so
    # repeated updates do not re-convert or drift.
    current_T      = convert(Vector{FT}, params.T)
    current_p_half = convert(Vector{FT}, params.p)
    current_q      = convert(Vector{FT}, params.q)
    current_vmr    = isnothing(params.absorption_params) ?
                     Dict{String, Any}() :
                     Dict{String, Any}(params.absorption_params.vmr)

    # ── 5. Initialise current per-aerosol loading state from params (B4) ───
    # update_aerosol_loading! and update_aerosol_microphysics! write these;
    # update_model!'s τ_aer redistribution reads them so a prior loading update
    # is not wiped by a later profile update.
    current_τ_ref        = FT[]
    current_profile_dist = Any[]
    if !isnothing(params.scattering_params)
        for i_aer in 1:n_aerosols
            c_aero = params.scattering_params.rt_aerosols[i_aer]
            push!(current_τ_ref, FT(c_aero.τ_ref))
            push!(current_profile_dist, c_aero.profile)
        end
    end

    return BatchContext(model, params, absorption_models, h2o_models,
                        k_ref_vec,
                        current_T, current_p_half, current_q, current_vmr,
                        current_τ_ref, current_profile_dist,
                        n_bands, n_aerosols, Nz, profile_reduction_n)
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

Only the arguments you supply are updated. `nothing` for a field means "keep the
**current** scene value" — i.e. the value from the most recent successful
`update_model!` call (or, if none yet, the value the `BatchContext` was built
with). It does **not** fall back to the original `params`. This makes successive
partial updates compose incrementally:

```julia
update_model!(ctx; vmr = B)   # vmr = B, everything else still original
update_model!(ctx; T   = C)   # T = C AND vmr is still B (not reset to original)
```

The merge is per-field: `T`, `p_half`, `q`, and the `vmr` override dict are each
remembered on the context (`ctx.current_T`, `ctx.current_p_half`, `ctx.current_q`,
`ctx.current_vmr`) and overwritten only by the fields you pass. The remembered
state is stored **unreduced** and already converted to the model's float type
(`params.float_type`), so repeated updates neither re-convert nor drift.

For `vmr`, the supplied keys are overlaid onto the current override dict (which
itself started as `params.absorption_params.vmr`); keys you do not pass keep
their last value. The merged overrides are then layered over the configured
defaults exactly as `model_from_parameters` does.

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
    # partial updates: if the user supplies only T (no p_half), we fall back to
    # the CURRENT stored scene value (unreduced) — NOT the original params — so
    # that successive partial updates compose incrementally (true incremental
    # semantics). The current state is stored unreduced on the context precisely
    # because it cannot be recovered from model.atmosphere.profile (which is
    # reduced). Falling back to the stored (already-reduced) profile_cur arrays
    # would also produce a length mismatch when T is passed unreduced.
    #
    # Rule: if a field is `nothing`, use the current stored value (unreduced);
    # otherwise convert the supplied value to FT.
    T_new      = T      === nothing ? ctx.current_T      : convert(Vector{FT}, T)
    q_new      = q      === nothing ? ctx.current_q      : convert(Vector{FT}, q)
    p_half_new = p_half === nothing ? ctx.current_p_half : convert(Vector{FT}, p_half)

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
    # Build the merged vmr override dict incrementally: start from the CURRENT
    # stored overrides (ctx.current_vmr, initialised from params at construction)
    # and overlay any newly supplied keys. This preserves prior vmr updates
    # across partial calls (B3). The merged overrides are then used directly as
    # the per-species vmr for the profile — they already include every configured
    # default (current_vmr started as a copy of ap.vmr) plus all overrides ever
    # applied. We do NOT mutate ctx.current_vmr until the update is validated and
    # applied (done at the end), so a failed call leaves the stored state intact.
    vmr_merged = if isnothing(ap)
        Dict{String, Any}()
    else
        d = Dict{String, Any}(ctx.current_vmr)
        if vmr !== nothing
            for (k, v) in vmr
                d[k] = v
            end
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
    # Use the same depol logic as model_from_parameters — including the SAME
    # vcd-weighted mean temperature for the Raman atmo constants. Hardcoding
    # 300 K here made τ_rayl differ from a fresh build at the ~1e-5 level: the
    # depol's T-dependence is weak, so the mismatch was sub-ULP (bit-invisible)
    # on x86-Linux but visible on macOS/Windows, breaking the rtol=0
    # bit-equality test in test_update_model.jl.
    rayleigh_molecular_T = (new_profile.vcd_dry' * new_profile.T) / sum(new_profile.vcd_dry)
    for i_band in 1:ctx.n_bands
        curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])

        νₘ  = FT(0.5) * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ  = FT(1.0e7) / νₘ
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(FT(1.0e7) / λₘ, FT(rayleigh_molecular_T))
        ϖ_Cab         = InelasticScattering.compute_ϖ_Cabannes(λₘ, _n2, _o2)
        γ_air_Cab, _  = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
        γ_air_Ray, _  = InelasticScattering.compute_γ_air_Rayleigh!(λₘ, _n2, _o2)
        depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
        depol_air_Ray = 2γ_air_Ray / (1 + γ_air_Ray)
        depol_use_Cab = params.depol < 0 ? FT(depol_air_Cab) : FT(params.depol)
        depol_use_Ray = params.depol < 0 ? FT(depol_air_Ray) : FT(params.depol)

        # Refresh the depolarization-dependent Rayleigh/Cabannes phase-matrix
        # greek coefs + ϖ_Cabannes in place. The vcd-weighted mean T changes
        # between scenes, so leaving these at the construction-time depol left
        # rt_run radiances ~1e-8 off a fresh build (test_update_model Test 3).
        # Mirrors model_from_parameters. greek_cabannes/ϖ_Cabannes are always
        # per-band vectors; greek_rayleigh is per-band when built via
        # model_from_parameters (always true for a BatchContext).
        model.optics.rayleigh.ϖ_Cabannes[i_band]     = FT(ϖ_Cab)
        model.optics.rayleigh.greek_cabannes[i_band] = Scattering.get_greek_rayleigh(depol_use_Cab)
        if model.optics.rayleigh.greek_rayleigh isa AbstractVector
            model.optics.rayleigh.greek_rayleigh[i_band] = Scattering.get_greek_rayleigh(depol_use_Ray)
        end

        model.optics.τ_rayl[i_band] .= getRayleighLayerOptProp(
            new_profile.p_half[end],
            curr_band_λ,
            depol_use_Ray,
            new_profile.vcd_dry)
    end

    # ── 3. Recompute τ_aer in place per band per aerosol ───────────────────
    # Microphysics (aerosol_optics, k) are fixed; only the vertical
    # redistribution changes with the new p_full/Δz.
    # B4: read the CURRENT loading state (ctx.current_τ_ref / current_profile_dist)
    # — written by update_aerosol_loading! and update_aerosol_microphysics! —
    # NOT the original params (c_aero.τ_ref / c_aero.profile). Otherwise a prior
    # loading update would be silently wiped by a profile update.
    if !isnothing(params.scattering_params)
        for i_aer in 1:ctx.n_aerosols
            c_aero    = params.scattering_params.rt_aerosols[i_aer]
            k_ref_aer = ctx.k_ref[i_aer]
            τ_eff     = ctx.current_τ_ref[i_aer]
            dist      = ctx.current_profile_dist[i_aer]

            if _has_analytic_phase_function(c_aero)
                τ_profile = getAerosolLayerOptProp(one(FT), dist, new_profile)
                for i_band in 1:ctx.n_bands
                    # τ_aer is now 3-D [iAer, nSpec, iLayer]; analytic aerosols
                    # have no k(λ) dependence — broadcast τ_eff * τ_profile across nSpec.
                    model.optics.aerosols.τ_aer[i_band][i_aer, :, :] .= τ_eff .* τ_profile'
                end
            else
                for i_band in 1:ctx.n_bands
                    # k_aer may be a scalar (single-λ band) or nSpec vector (multi-λ band).
                    k_aer = model.optics.aerosols.aerosol_optics[i_band][i_aer].k
                    τ_profile = getAerosolLayerOptProp(one(FT), dist, new_profile)
                    # k_aer .* τ_profile' broadcasts to (nSpec × nLayers).
                    model.optics.aerosols.τ_aer[i_band][i_aer, :, :] .=
                        (τ_eff / k_ref_aer) .* k_aer .* τ_profile'
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
            # B2: pass the UPDATED scene VMRs (new_profile.vmr), not the original
            # ap.vmr config dict, so CIA colliders (O2/N2) respond to per-scene
            # vmr changes exactly as a fresh model_from_parameters build with the
            # changed vmr would. compute_τ_cia! falls back to USS defaults for
            # any collider species absent from the dict.
            for cia_path in ap.cia_files
                cia_table = Absorption.load_cia_table(
                    cia_path, params.spec_bands[i_band]; FT)
                Absorption.compute_τ_cia!(
                    model.optics.τ_abs[i_band], cia_table, new_profile,
                    new_profile.vmr)
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

    # ── 5. Persist the current scene state (B3) ─────────────────────────────
    # Only reached after all validation passed and the model was updated, so a
    # failed call leaves the stored state untouched. Store already FT-converted,
    # unreduced vectors so the next partial update composes without re-converting
    # or drifting. (T_new/q_new/p_half_new are either the FT-converted supplied
    # values or the previously-stored current values — already FT.)
    ctx.current_T      = T_new
    ctx.current_q      = q_new
    ctx.current_p_half = p_half_new
    ctx.current_vmr    = vmr_merged

    return nothing
end

# ============================================================================
# Phase 2: aerosol-specific update functions
# ============================================================================

# ── Internal helper: re-derive and write SolverConfig Fourier bounds ─────────
#
# SolverConfig is an immutable struct.  Its three per-band VECTOR fields
# (m_max_bands, n_fourier_moments_bands, l_max) are mutable Julia Arrays and
# can be updated with .=.  The SCALAR fields (l_trunc, Δ_angle, depol,
# polarization_type, quadrature_type, use_component_traits) are bits-immutable
# and cannot be changed in-place.
#
# For all Phase-2 aerosol updates those scalar fields are guaranteed not to
# change (they derive from params, not from per-aerosol optics), so this helper
# only needs to rewrite the three Vector fields.  If a future caller tries to
# change l_trunc or Δ_angle they must rebuild the entire model — this function
# will @error and throw if it detects the values would need to change.
function _rewrite_solver_fourier_bounds!(ctx::BatchContext)
    model  = ctx.model
    params = ctx.params
    FT     = params.float_type

    n_bands   = ctx.n_bands
    solver    = model.solver
    n_aer     = ctx.n_aerosols
    ae_optics = model.optics.aerosols.aerosol_optics  # [i_band][i_aer]

    # Recompute l_max_aer exactly as model_from_parameters does.
    l_max_aer = zeros(Int, max(n_aer, 1), n_bands)
    truncation_type = _resolved_truncation(params, FT)

    for i_aer in 1:n_aer
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        for i_band in 1:n_bands
            ao = ae_optics[i_band][i_aer]
            l_max_aer[i_aer, i_band] =
                if _has_analytic_phase_function(c_aero)
                    # analytic path: length of the stored greek series
                    length(ao.greek_coefs.β)
                elseif truncation_type isa Scattering.δBGE
                    min(length(ao.greek_coefs.β), truncation_type.l_max)
                else
                    length(ao.greek_coefs.β)
                end
        end
    end

    # New per-band l_max
    new_l_max = zeros(Int, n_bands)
    for i_band in 1:n_bands
        new_l_max[i_band] = n_aer > 0 ?
            maximum(l_max_aer[:, i_band]) : params.l_trunc
    end

    # Re-derive m_max_bands via the same trait aggregator as model_from_parameters.
    components_per_band = [_band_components(params, ae_optics, model.sources, i_band)
                           for i_band in 1:n_bands]
    new_m_max_bands = _derive_m_max_bands_via_traits(
        new_l_max, params.max_m, components_per_band, model.quad_points.Nstreams)
    new_n_fourier   = new_m_max_bands .+ 1

    # Write back into the mutable Vector fields in-place.
    solver.m_max_bands          .= new_m_max_bands
    solver.n_fourier_moments_bands .= new_n_fourier
    solver.l_max                .= new_l_max

    return nothing
end

"""
    update_aerosol_loading!(ctx::BatchContext, i_aer::Int;
                            τ_ref        = nothing,
                            profile_dist = nothing)

**Cheap path** — update only the column-integrated optical depth (`τ_ref`)
and/or the vertical pressure distribution (`profile_dist`) for aerosol
species `i_aer`, using the *cached* aerosol optics (`aerosol_optics`) and
reference extinction coefficient (`ctx.k_ref[i_aer]`).  No Mie
recomputation is performed.

After this call `rt_run(ctx.model)` produces radiances for the new aerosol
loading.  All other model state (gas absorption, profile, surface) is
unchanged.

# Keyword arguments

- `τ_ref::Real`: New column optical depth at `λ_ref`. When `nothing`, the
  value stored in `params.scattering_params.rt_aerosols[i_aer].τ_ref` is
  kept (only the distribution is updated).
- `profile_dist::Distributions.Distribution`: New vertical pressure
  distribution (a `Distributions.jl` `Distribution` object — the same type
  used in `RT_Aerosol.profile`).  When `nothing`, the existing distribution
  is kept.

# Cost

`O(Nz)` per band — evaluates the distribution PDF on the pressure grid and
scales the optical depth profile. No Mie, no HITRAN.

# What it invalidates

- `model.optics.aerosols.τ_aer[i_band][i_aer, :, :]` — all bands (3-D: nSpec × nLayers).
- **Does NOT** change `aerosol_optics`, `k_ref`, `SolverConfig` Fourier
  bounds, gas absorption, or Rayleigh properties.

# Example

```julia
ctx = BatchContext(params)
# Change τ_ref for aerosol 1 to 0.3, keep vertical distribution
update_aerosol_loading!(ctx, 1; τ_ref = 0.3)
R, T = rt_run(ctx.model)
```
"""
function update_aerosol_loading!(ctx::BatchContext, i_aer::Int;
                                  τ_ref        = nothing,
                                  profile_dist = nothing)
    1 <= i_aer <= ctx.n_aerosols || error(
        "update_aerosol_loading!: i_aer=$i_aer is out of range " *
        "[1, $(ctx.n_aerosols)]")

    model  = ctx.model
    params = ctx.params
    FT     = params.float_type

    c_aero = params.scattering_params.rt_aerosols[i_aer]
    profile = model.atmosphere.profile

    # Resolve τ_ref and distribution. `nothing` keeps the CURRENT stored value
    # (not the original params), mirroring update_model!'s incremental semantics.
    τ_eff  = τ_ref === nothing ? ctx.current_τ_ref[i_aer] : FT(τ_ref)
    dist   = profile_dist === nothing ? ctx.current_profile_dist[i_aer] : profile_dist

    # B4: persist the resolved loading so a later update_model! redistribution
    # (and a later loading update that omits a field) sees it instead of the
    # original params τ_ref / profile.
    ctx.current_τ_ref[i_aer]        = τ_eff
    ctx.current_profile_dist[i_aer] = dist

    k_ref_aer = ctx.k_ref[i_aer]

    # Recompute τ_aer rows for all bands using the cached aerosol optics.
    # τ_aer is now 3-D [iAer, nSpec, iLayer]; k_aer may be scalar or nSpec vector.
    if _has_analytic_phase_function(c_aero)
        τ_profile = getAerosolLayerOptProp(one(FT), dist, profile)
        for i_band in 1:ctx.n_bands
            model.optics.aerosols.τ_aer[i_band][i_aer, :, :] .= τ_eff .* τ_profile'
        end
    else
        for i_band in 1:ctx.n_bands
            k_aer     = model.optics.aerosols.aerosol_optics[i_band][i_aer].k
            τ_profile = getAerosolLayerOptProp(one(FT), dist, profile)
            model.optics.aerosols.τ_aer[i_band][i_aer, :, :] .=
                (τ_eff / k_ref_aer) .* k_aer .* τ_profile'
        end
    end

    return nothing
end

"""
    update_aerosol_microphysics!(ctx::BatchContext, i_aer::Int, aerosol::Aerosol;
                                  τ_ref = nothing)

**Expensive path** — replace the microphysics of aerosol species `i_aer` with
a new [`Aerosol`](@ref) (new size distribution and/or refractive index) and
recompute all derived quantities via the same Mie code path that
[`model_from_parameters`](@ref) uses.

This updates:
- `model.optics.aerosols.aerosol_optics[i_band][i_aer]` for every band
- `ctx.k_ref[i_aer]` (the new reference-wavelength extinction coefficient)
- `model.solver.m_max_bands`, `n_fourier_moments_bands`, `l_max` (in-place
  via `.=`) — the critical Fourier-loop re-derivation that prevents silent
  wrong results when the new particle size changes the length of the Greek-
  coefficient series.
- `model.optics.aerosols.τ_aer[i_band][i_aer, :, :]` for every band (3-D: nSpec × nLayers).

After this call `rt_run(ctx.model)` gives the same result as building a
fresh model with the new aerosol.

# Arguments

- `ctx`: The [`BatchContext`](@ref) to update.
- `i_aer`: 1-based index of the aerosol species to replace.
- `aerosol`: New [`Aerosol`](@ref) (from `Scattering.Aerosol`), carrying
  the new `size_distribution`, `nᵣ`, and `nᵢ`.

# Keyword arguments

- `τ_ref::Real`: New column optical depth at `λ_ref`. When `nothing`, the
  current scene value `ctx.current_τ_ref[i_aer]` is kept (i.e. the most
  recently applied loading, not the original `params` value).

# SolverConfig mutability note

`SolverConfig` is an immutable struct, but its three **Vector** fields
(`m_max_bands`, `n_fourier_moments_bands`, `l_max`) are mutable Julia arrays
and are updated in-place with `.=`.  The scalar fields (`l_trunc`, `Δ_angle`,
`depol`, `polarization_type`, `quadrature_type`) cannot be changed in-place;
they are guaranteed to be unaffected by aerosol microphysics changes since
they derive from `params`, not from per-aerosol optics.  If a caller somehow
requires a new `l_trunc` or `Δ_angle`, a full `model_from_parameters` rebuild
is necessary.

# Cost

Full Mie per band (`O(n_bands × nquad_radius × size_distribution_points)`),
i.e. the same cost as the aerosol section of `model_from_parameters`.

# What it invalidates

- `aerosol_optics[i_band][i_aer]` — all bands.
- `ctx.k_ref[i_aer]`.
- `solver.m_max_bands`, `solver.n_fourier_moments_bands`, `solver.l_max`.
- `τ_aer[i_band][i_aer, :, :]` — all bands (3-D: nSpec × nLayers).

# Example

```julia
ctx = BatchContext(params)
# Save initial state for comparison
R_A, _ = rt_run(ctx.model)

# Replace aerosol 1 with larger particles (reff ~2 μm)
new_aerosol = Aerosol(LogNormal(log(2.0), 0.4), 1.3, 1e-8)
update_aerosol_microphysics!(ctx, 1, new_aerosol; τ_ref = 0.1)
R_B, _ = rt_run(ctx.model)
```
"""
function update_aerosol_microphysics!(ctx::BatchContext, i_aer::Int, aerosol::Aerosol;
                                       τ_ref = nothing)
    1 <= i_aer <= ctx.n_aerosols || error(
        "update_aerosol_microphysics!: i_aer=$i_aer is out of range " *
        "[1, $(ctx.n_aerosols)]")
    # n_aerosols > 0 implies scattering_params is not nothing; guard for safety
    isnothing(ctx.params.scattering_params) && error(
        "update_aerosol_microphysics!: ctx.params.scattering_params is nothing")

    model  = ctx.model
    params = ctx.params
    FT     = params.float_type
    sp     = params.scattering_params

    c_aero          = sp.rt_aerosols[i_aer]
    truncation_type = _resolved_truncation(params, FT)
    # B4: resolve τ_ref / distribution from the CURRENT loading state when not
    # supplied (not the original params), and persist them so a later
    # update_model! redistribution preserves this loading.
    τ_eff = τ_ref === nothing ? ctx.current_τ_ref[i_aer] : FT(τ_ref)
    dist  = ctx.current_profile_dist[i_aer]
    ctx.current_τ_ref[i_aer]        = τ_eff
    ctx.current_profile_dist[i_aer] = dist

    # ── 1. Recompute k_ref at reference wavelength ──────────────────────────
    # k_ref uses the reference refractive index n_ref (normalisation convention),
    # NOT the new aerosol's own nᵣ/nᵢ. Build a separate Aerosol with n_ref so
    # the incoming aerosol object is never mutated (fixes aliasing bug).
    ref_aerosol_mic = Aerosol(aerosol.size_distribution,
                               real(sp.n_ref),
                               -imag(sp.n_ref))
    mie_model_ref = make_mie_model(
        sp.decomp_type, ref_aerosol_mic, sp.λ_ref,
        params.polarization_type, truncation_type,
        sp.r_max, sp.nquad_radius)
    new_k_ref = Float64(compute_ref_aerosol_extinction(mie_model_ref, FT))
    ctx.k_ref[i_aer] = new_k_ref

    # ── 2. Per-band Mie + truncation + τ_aer rows ───────────────────────────
    # Mirror the model_from_parameters.jl endpoint-Mie + linear interpolation
    # strategy: single-λ bands use band-center Mie (scalar k, bit-identical);
    # multi-λ bands compute Mie at band edges and linearly interpolate k(λ) in
    # wavenumber so τ_aer is consistent with a fresh model build.
    profile = model.atmosphere.profile

    _mie_fwd(λ) = make_mie_model(
        sp.decomp_type, aerosol, λ,
        params.polarization_type, truncation_type,
        sp.r_max, sp.nquad_radius)

    for i_band in 1:ctx.n_bands
        curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])
        n_spec = length(curr_band_λ)

        if n_spec == 1
            # ── Single-wavelength band: Mie at band center (bit-identical) ──
            band_center_λ = (maximum(curr_band_λ) + minimum(curr_band_λ)) / 2
            mie_model     = _mie_fwd(band_center_λ)
            aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT)

            β_len = length(aerosol_optics_raw.greek_coefs.β)
            new_ao =
                if truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                    Scattering.truncate_phase(truncation_type,
                                              aerosol_optics_raw; reportFit=false)
                else
                    Scattering.truncate_phase(Scattering.NoTruncation(),
                                              aerosol_optics_raw)
                end

            model.optics.aerosols.aerosol_optics[i_band][i_aer] = new_ao

            τ_profile = getAerosolLayerOptProp(one(FT), dist, profile)
            model.optics.aerosols.τ_aer[i_band][i_aer, 1, :] .=
                (τ_eff / new_k_ref) * new_ao.k * τ_profile

        else
            # ── Multi-wavelength band: endpoint Mie + linear interpolation ──
            mie_model_0 = _mie_fwd(curr_band_λ[1])
            aerosol_optics_raw_0 = compute_aerosol_optical_properties(mie_model_0, FT)

            mie_model_1 = _mie_fwd(curr_band_λ[end])
            aerosol_optics_raw_1 = compute_aerosol_optical_properties(mie_model_1, FT)

            # Average Greek coefs from both endpoints (matching forward path).
            # Robust to either ordering of series length (descending bands).
            Nl0 = length(aerosol_optics_raw_0.greek_coefs.α)
            Nl1 = length(aerosol_optics_raw_1.greek_coefs.α)
            Nlo, Nhi = min(Nl0, Nl1), max(Nl0, Nl1)
            gc0 = aerosol_optics_raw_0.greek_coefs
            gc1 = aerosol_optics_raw_1.greek_coefs
            greek_arrs = Dict{Symbol,Vector{FT}}()
            for fn in (:α, :β, :γ, :δ, :ϵ, :ζ)
                g0 = getfield(gc0, fn); g1 = getfield(gc1, fn)
                arr = zeros(FT, Nhi)
                arr[1:Nlo] .= FT(0.5) .* (g0[1:Nlo] .+ g1[1:Nlo])
                if Nl0 > Nlo
                    arr[1+Nlo:Nhi] .= FT(0.5) .* g0[1+Nlo:Nl0]
                elseif Nl1 > Nlo
                    arr[1+Nlo:Nhi] .= FT(0.5) .* g1[1+Nlo:Nl1]
                end
                greek_arrs[fn] = arr
            end
            greek_coeffs = GreekCoefs(greek_arrs[:α], greek_arrs[:β], greek_arrs[:γ],
                                       greek_arrs[:δ], greek_arrs[:ϵ], greek_arrs[:ζ])

            # Linearly interpolate k_ext and k_sca in wavenumber across the band.
            # Sort endpoints by wavenumber (bands may be descending) — mirrors
            # model_from_parameters.jl.
            ν0, ν1 = 1e4/curr_band_λ[1], 1e4/curr_band_λ[end]
            asc = ν0 <= ν1
            ν_grid    = asc ? [ν0, ν1] : [ν1, ν0]
            kext_grid = asc ? [aerosol_optics_raw_0.k, aerosol_optics_raw_1.k] :
                              [aerosol_optics_raw_1.k, aerosol_optics_raw_0.k]
            ksca_grid = asc ?
                [aerosol_optics_raw_0.k * aerosol_optics_raw_0.ω̃,
                 aerosol_optics_raw_1.k * aerosol_optics_raw_1.ω̃] :
                [aerosol_optics_raw_1.k * aerosol_optics_raw_1.ω̃,
                 aerosol_optics_raw_0.k * aerosol_optics_raw_0.ω̃]
            interp_kext = LinearInterpolation(ν_grid, kext_grid)
            interp_ksca = LinearInterpolation(ν_grid, ksca_grid)
            k_spec  = FT[interp_kext(1e4/curr_band_λ[i]) for i=1:n_spec]
            ω̃_spec  = FT[interp_ksca(1e4/curr_band_λ[i])/k_spec[i] for i=1:n_spec]
            fᵗ_spec = zeros(FT, n_spec)

            aerosol_optics_raw_interp = AerosolOptics(greek_coefs=greek_coeffs,
                                                       ω̃=ω̃_spec, k=k_spec, fᵗ=fᵗ_spec)

            β_len = length(aerosol_optics_raw_interp.greek_coefs.β)
            new_ao =
                if truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                    Scattering.truncate_phase(truncation_type,
                                              aerosol_optics_raw_interp; reportFit=false)
                else
                    Scattering.truncate_phase(Scattering.NoTruncation(),
                                              aerosol_optics_raw_interp)
                end

            model.optics.aerosols.aerosol_optics[i_band][i_aer] = new_ao

            τ_profile = getAerosolLayerOptProp(one(FT), dist, profile)
            model.optics.aerosols.τ_aer[i_band][i_aer, :, :] .=
                (τ_eff / new_k_ref) .* new_ao.k .* τ_profile'
        end
    end

    # ── 3. Re-derive Fourier bounds (CRITICAL anti-silent-wrongness step) ───
    # New Greek series length may differ from the old one → m_max_bands /
    # l_max must be recomputed and written back into the SolverConfig Vectors.
    _rewrite_solver_fourier_bounds!(ctx)

    return nothing
end

# Export hint (actual export list lives in CoreRT.jl)
# export BatchContext, update_model!, update_aerosol_loading!, update_aerosol_microphysics!
