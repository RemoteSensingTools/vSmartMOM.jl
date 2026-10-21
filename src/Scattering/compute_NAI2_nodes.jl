#=
=====================================================================
Caller-defined-size-distribution ("caller-node") bulk Mie API.
=====================================================================

`compute_aerosol_optical_properties_nodes` is the public entry point: given
an arbitrary set of (radius, weight) nodes — e.g. a sectional/tabulated size
distribution such as GCHPIO's TOMAS two-moment bins — plus one complex
refractive index for the whole ensemble, it returns the BULK AerosolOptics
(number-mean cross sections, ω̃, and all six Greek coefficient arrays), using
exactly the same Siewert NAI-2 physics and Hovenier/Greek conventions as the
log-normal path in `compute_NAI2.jl`.

Key differences from the log-normal (`MieModel`) entry point:
- The caller supplies radii/weights directly; weights need not be normalized
  (they are normalized internally to number-mean properties).
- The angular resolution (`n_max`/`n_mu`) is derived from the ACTUAL largest
  supplied radius — no artificial `r_max` cap.
- The GPU path launches the previously-unused device-resident reduction and
  Greek-projection kernels (`size_reduction_kernel!`, `weighted_sum_kernel!`,
  `greek_coefficients_kernel!` in `gpu_mie_kernels.jl`) so that only the bulk
  Greek arrays + ω̃/k are copied back to the host — no per-node aₙ/bₙ or
  angular arrays are cached or returned.

This is the intended seam for sectional/tabulated aerosol size distributions;
GCHPIO's TOMAS two-moment bins are the first consumer (see
`docs/dev_notes/proposals/gchp_aerosol_optics/`).
=====================================================================
=#

"""
    _slice_greek(a::AerosolOptics, l_max::Int) -> AerosolOptics

Truncate (or leave unchanged, if `l_max` is already ≥ the natural length) all
six Greek coefficient arrays of `a` to `min(l_max, length(a.greek_coefs.β))`
terms. Plain array slicing — no re-fitting (`δBGE` re-fitting is a distinct,
separate step, see [`truncate_phase`](@ref)).
"""
function _slice_greek(a::AerosolOptics, l_max::Integer)
    gc = a.greek_coefs
    L = min(l_max, length(gc.β))
    greek_coefs = GreekCoefs(gc.α[1:L], gc.β[1:L], gc.γ[1:L],
                             gc.δ[1:L], gc.ϵ[1:L], gc.ζ[1:L])
    return AerosolOptics(greek_coefs=greek_coefs, ω̃=a.ω̃, k=a.k, fᵗ=a.fᵗ, derivs=a.derivs)
end

"""
    _apply_requested_truncation(raw::AerosolOptics, truncation::AbstractTruncationType) -> AerosolOptics

Apply `truncation` the same way the RT pipeline does downstream of the raw
Mie call (see `src/CoreRT/tools/model_from_parameters.jl` /
`update_model.jl`): a `δBGE` truncation is actually fit only if it would
shorten the Greek series (`length(β) > truncation.l_max`); otherwise (and for
`NoTruncation`/`AutoTruncation`) the identity `truncate_phase(NoTruncation(), raw)`
is applied, which resets `fᵗ` from the raw "untruncated yet" sentinel (`1`) to
the true no-op value (`0`).
"""
function _apply_requested_truncation(raw::AerosolOptics, truncation::AbstractTruncationType)
    β_len = length(raw.greek_coefs.β)
    if truncation isa δBGE && β_len > truncation.l_max
        return truncate_phase(truncation, raw)
    else
        return truncate_phase(NoTruncation(), raw)
    end
end

@doc raw"""
    compute_aerosol_optical_properties_nodes(radii, weights, n_real, n_imag,
        wavelength, polarization, truncation;
        architecture = Architectures.CPU(),
        l_max = nothing,
        precision_policy = nothing) -> AerosolOptics

Bulk Mie aerosol optics for a CALLER-DEFINED set of size-distribution nodes,
using the Siewert NAI-2 formulation (same physics/conventions as
[`compute_aerosol_optical_properties`](@ref) for a log-normal `MieModel`, and
the same [`AerosolOptics`](@ref)/[`GreekCoefs`](@ref) output shape).

This is the intended seam for **sectional or tabulated size distributions**
that are not naturally expressed as a `Distributions.jl` continuous
distribution — e.g. GCHPIO's TOMAS two-moment sectional bins. No log-normal
(or any other analytic) size distribution is assumed anywhere in this path.

# Arguments
- `radii::AbstractVector`: wet particle radii, one per node, **[μm]**.
- `weights::AbstractVector`: non-negative node weights (number-mixing-ratio
  or number-concentration weights). **Need not be normalized** — they are
  divided by their sum internally, so the returned bulk quantities are
  number-mean properties and are invariant to a uniform rescaling of
  `weights` (e.g. `weights .* 7.3` gives the identical `AerosolOptics`).
- `n_real, n_imag`: the ONE complex refractive index shared by the whole
  node ensemble, using vSmartMOM's convention `n = n_real - i·n_imag`
  (`n_imag ≥ 0`; positive `n_imag` is absorption — same convention as
  [`Aerosol`](@ref) and asserted identically).
- `wavelength`: wavelength, same units as `radii` (μm internally).
- `polarization::AbstractPolarizationType`: accepted for API parity with
  [`make_mie_model`](@ref) / `MieModel`; as in the existing NAI2 path, all six
  Greek coefficient arrays are always computed regardless of polarization
  type (the type does not currently gate the computation — matching
  `compute_aerosol_optical_properties`'s existing behavior).
- `truncation::AbstractTruncationType`: applied to the raw bulk Greek
  coefficients via [`truncate_phase`](@ref), exactly mirroring how the RT
  pipeline (`model_from_parameters`/`update_model`) applies a `MieModel`'s
  `truncation_type` downstream of the raw Mie call: a `δBGE` is fit only if it
  would shorten the natural (node-driven) Greek series; otherwise (including
  `NoTruncation`) the untruncated series passes through with `fᵗ = 0`.

# Keyword arguments
- `architecture`: `Architectures.CPU()` (default) or a GPU architecture with a
  registered Mie pipeline (`Architectures.has_gpu_mie`, e.g. CUDA `GPU()`) —
  same architecture routing as `compute_aerosol_optical_properties(model)`.
  Non-CPU architectures without a GPU Mie pipeline warn once and fall back to
  CPU (mirroring the existing `MieModel` router).
- `l_max`: optional cap on the number of OUTPUT Greek terms, applied by plain
  array slicing (`min(l_max, series_length)`) AFTER `truncation` is applied —
  a `δBGE` fit always sees the full natural series, and the cap then shortens
  whatever the truncation produced. This is an efficiency knob, independent of
  `truncation`: because angular resolution here scales with the largest
  supplied node (no `r_max` cap), a very wide node ensemble can produce a long
  natural Greek series; `l_max` lets a caller who does not need the full
  series avoid paying for it downstream, without engaging `δBGE`'s
  forward-peak re-fit. `nothing` (default) keeps the full length.
- `precision_policy`: GPU-only, see [`MiePrecisionPolicy`](@ref)
  (`NativeFloat64`/`DSEmulated`/`NativeFloat32`); `nothing` (default)
  auto-selects via `Architectures.default_mie_precision_policy` exactly as the
  `MieModel` GPU path does (`NativeFloat32` on `MetalGPU()`, since Metal has no
  Float64 hardware at all). Ignored on the CPU path.

# Angular resolution
`n_max` is `get_n_max(maximum(2π .* radii ./ wavelength))` — i.e. derived from
the size parameter of the ACTUAL largest supplied radius, with **no
artificial cap**. `n_mu = 2n_max - 1`, exactly as the existing NAI2 path
computes it from `r_max` when the size distribution comes from a
log-normal `MieModel`.

# Returns
[`AerosolOptics`](@ref) with `greek_coefs`, number-mean `ω̃`, number-mean `k`,
and `fᵗ` from the requested `truncation`.
"""
function compute_aerosol_optical_properties_nodes(
    radii::AbstractVector, weights::AbstractVector,
    n_real::Real, n_imag::Real,
    wavelength::Real,
    polarization::AbstractPolarizationType,
    truncation::AbstractTruncationType;
    architecture::Architectures.AbstractArchitecture = Architectures.CPU(),
    l_max::Union{Nothing,Integer} = nothing,
    precision_policy::Union{Nothing,MiePrecisionPolicy} = nothing,
)
    _validate_node_inputs(radii, weights, n_imag)
    @assert sum(weights) > 0 "weights must not sum to zero"

    # Output float type: promoted user type of all numeric inputs (mirrors
    # _mie_output_type's role for the MieModel path).
    FT = float(promote_type(eltype(radii), eltype(weights),
                            typeof(n_real), typeof(n_imag), typeof(wavelength)))

    raw = if architecture isa Architectures.CPU
        _aerosol_optical_properties_nodes_cpu(radii, weights, n_real, n_imag, wavelength, FT)
    elseif Architectures.has_gpu_mie(architecture)
        backend = Architectures.ka_backend(architecture)
        policy = precision_policy === nothing ?
                 Architectures.default_mie_precision_policy(architecture, FT) :
                 precision_policy
        compute_aerosol_optical_properties_nodes_gpu(radii, weights, n_real, n_imag,
                                                      wavelength, backend;
                                                      precision_policy = policy, l_max = nothing)
    else
        if !_WARNED_NO_GPU_MIE_NODES[]
            @warn "no GPU Mie pipeline for $(architecture); computing caller-node Mie on CPU"
            _WARNED_NO_GPU_MIE_NODES[] = true
        end
        _aerosol_optical_properties_nodes_cpu(radii, weights, n_real, n_imag, wavelength, FT)
    end

    # Truncation FIRST (on the full raw series — a δBGE fit must see every
    # natural Greek term), output-length cap AFTER. The reverse order would
    # let a small `l_max` silently disable or distort a requested δBGE fit
    # (the fit would run on — or be skipped because of — a pre-sliced series).
    out = _apply_requested_truncation(raw, truncation)
    return l_max === nothing ? out : _slice_greek(out, l_max)
end

# One-time warning state, mirroring `_WARNED_NO_GPU_MIE` in phase_function_autodiff.jl
# but scoped to the node API (kept independent so tests for one don't clear the other).
const _WARNED_NO_GPU_MIE_NODES = Ref(false)

"""
    _validate_node_inputs(radii, weights, n_imag)

Shared input validation for BOTH exported caller-node entry points
([`compute_aerosol_optical_properties_nodes`](@ref) and
[`compute_aerosol_optical_properties_nodes_gpu`](@ref)) — one list, so the two
public seams cannot drift apart (a negative weight passed only to the `_gpu`
entry point would otherwise produce silently wrong bulk optics rather than an
error). `@assert` matches this module's established input-validation
convention (see `compute_NAI2.jl`/`compute_NAI2_gpu.jl`).
"""
@inline function _validate_node_inputs(radii, weights, n_imag)
    @assert length(radii) == length(weights) "radii and weights must have the same length"
    @assert length(radii) ≥ 1 "need at least one size-distribution node"
    @assert all(w -> w ≥ 0, weights) "weights must be ≥ 0"
    @assert n_imag ≥ 0 "Imaginary part of the refractive index must be ≥ 0 (convention n = n_real - i·n_imag)"
    return nothing
end

"""
    _aerosol_optical_properties_nodes_cpu(radii, weights, n_real, n_imag, λ, FT) -> AerosolOptics

CPU implementation backing [`compute_aerosol_optical_properties_nodes`](@ref).
Normalizes `weights`, derives `(n_max, n_mu)` from the actual node radii (no
cap), and delegates to the shared `_nai2_bulk_optics` core — the same
core the log-normal `MieModel` NAI2 path uses. Internal computation always
runs in Float64 (`IC`), narrowed to `FT` only at the end, exactly mirroring
`compute_aerosol_optical_properties`'s IC/FT_out convention.
"""
function _aerosol_optical_properties_nodes_cpu(radii, weights, n_real, n_imag, λ, FT::Type)
    IC = FT <: AbstractFloat ? Float64 : FT

    r = IC.(radii)
    w = IC.(weights)          # broadcast always materializes a fresh copy, so
    w ./= sum(w)              # in-place normalization never aliases the input

    k = IC(2π) / IC(λ)
    x_size_param = k .* r
    n_max = get_n_max(maximum(x_size_param))
    n_mu  = 2n_max - 1

    core = _nai2_bulk_optics(r, w, IC(n_real), IC(n_imag), IC(λ), n_max, n_mu, IC)

    if IC <: AbstractFloat && FT <: AbstractFloat
        greek_coefs = GreekCoefs(convert.(FT, core.α), convert.(FT, core.β),
                                 convert.(FT, core.γ), convert.(FT, core.δ),
                                 convert.(FT, core.ϵ), convert.(FT, core.ζ))
        return AerosolOptics(greek_coefs=greek_coefs,
                             ω̃=FT(core.bulk_C_sca / core.bulk_C_ext),
                             k=FT(core.bulk_C_ext), fᵗ=FT(1))
    else
        greek_coefs = GreekCoefs(core.α, core.β, core.γ, core.δ, core.ϵ, core.ζ)
        return AerosolOptics(greek_coefs=greek_coefs,
                             ω̃=(core.bulk_C_sca / core.bulk_C_ext),
                             k=(core.bulk_C_ext), fᵗ=one(eltype(core.α)))
    end
end

"""
    _prepare_node_mie_gpu(radii, weights, n_real, n_imag, wavelength, backend,
                          precision_policy) -> NamedTuple

Shared GPU preamble for BOTH [`compute_aerosol_optical_properties_nodes_gpu`](@ref)
and [`compute_aerosol_extinction_nodes_gpu`](@ref): validates inputs, derives
the output float type `FT` and the size-distribution-reduction accumulator
type `RA` (policy-dependent, see [`MiePrecisionPolicy`](@ref) /
[`_mie_reduction_type`](@ref)), builds the per-node size parameters and Mie
recursion inputs, allocates the per-node device arrays, and launches Kernel 1
(Mie coefficients aₙ, bₙ). Everything downstream of Kernel 1
(amplitude/phase-matrix/cross-section/reduction/Greek kernels) differs between
the two callers and is NOT part of this shared preamble.

`RA` is the key lever for [`NativeFloat32`](@ref): unlike
[`NativeFloat64`](@ref) (`RA = Float64` always) and [`DSEmulated`](@ref)
(`RA = Float64` when `FT === Float32`, matching the historical Float64-widened
+ Neumaier-compensated reduction), `NativeFloat32` sets `RA = Float32` — so the
caller allocates Float32 reduction/Greek/Legendre device arrays instead of
Float64 ones. This is what makes the node GPU pipeline allocate ZERO Float64
device arrays under `NativeFloat32` (required for Apple Metal, which has no
Float64 hardware support at all, and for consumer CUDA with crippled FP64
throughput).

On a Metal backend (`_is_metal_backend`, a capability trait defaulting to
`false` on `Any` and refined to `true` for `Metal.MetalBackend` by the weakly
loaded `vSmartMOMMetalExt` extension — real dispatch, not name/string
matching, the same precompile-safe pattern as `Architectures.ka_backend` /
`has_gpu_mie`), this throws `ArgumentError` if `FT === Float64` or
`RA === Float64` — i.e. Metal only supports `NativeFloat32` this round
(`DSEmulated`-on-Metal, which would need Float32-COMPENSATED reductions
instead of the Float64-widened ones used elsewhere, is not implemented yet;
see [`MiePrecisionPolicy`](@ref)).
"""
function _prepare_node_mie_gpu(radii::AbstractVector, weights::AbstractVector,
                                n_real::Real, n_imag::Real, wavelength::Real, backend,
                                precision_policy::MiePrecisionPolicy)
    _validate_node_inputs(radii, weights, n_imag)

    FT = float(promote_type(eltype(radii), eltype(weights),
                            typeof(n_real), typeof(n_imag), typeof(wavelength)))
    _check_policy_ft(precision_policy, FT)
    RA = _mie_reduction_type(precision_policy, FT)

    if _is_metal_backend(backend) && (FT === Float64 || RA === Float64)
        throw(ArgumentError(
            "Metal has no Float64 hardware support: $(nameof(typeof(precision_policy))) " *
            "with node element type $FT would need a Float64 device array (the " *
            "Mie-coefficient kernel and/or the size-distribution reduction). On " *
            "MetalGPU(), use Float32 node inputs with NativeFloat32() -- NativeFloat64() " *
            "is unsupported on Metal entirely, and DSEmulated-on-Metal (a Float32-" *
            "COMPENSATED reduction, instead of the Float64-widened one used elsewhere) " *
            "is not implemented yet (tracked as a follow-up)."))
    end

    r = FT.(radii)
    # Normalize in Float64 BEFORE narrowing to the kernel type: a Float32 sum
    # of large-but-finite weights can overflow to Inf (zeroing every weight and
    # poisoning ω̃/Greek output with NaNs), and normalized weights are O(1) so
    # the narrowed copy is safe. The CPU path gets this for free via IC=Float64.
    w64 = Float64.(weights)
    w_sum = sum(w64)
    @assert isfinite(w_sum) && w_sum > 0 "weights must sum to a positive, finite value"
    w64 ./= w_sum
    w = FT.(w64)

    k_wavenum = FT(2π / wavelength)
    x_size_param = k_wavenum .* r
    n_max_global = get_n_max(maximum(x_size_param))
    nquad_radius = length(r)

    m_ref_re = FT(n_real)
    m_ref_im = FT(-n_imag)  # convention: n = n_real - i*n_imag, so m_im is negative

    nmax_per_r = [get_n_max(x) for x in x_size_param]
    y_max = maximum(x_size_param) * sqrt(m_ref_re^2 + m_ref_im^2)
    nmx_max = round(Int, max(n_max_global, y_max) + 51)

    AT = KernelAbstractions.allocate

    # --- Per-node arrays: native kernel precision FT (speed) ---
    x_dev    = AT(backend, FT, nquad_radius)
    nmax_dev = AT(backend, Int, nquad_radius)
    an_dev   = AT(backend, Complex{FT}, (nquad_radius, n_max_global))
    bn_dev   = AT(backend, Complex{FT}, (nquad_radius, n_max_global))

    KernelAbstractions.copyto!(backend, x_dev, FT.(x_size_param))
    KernelAbstractions.copyto!(backend, nmax_dev, nmax_per_r)

    fill!(an_dev, zero(Complex{FT}))
    fill!(bn_dev, zero(Complex{FT}))

    # --- Kernel 1: Mie coefficients (unchanged, existing kernel; policy selects
    #     mie_coefficients_kernel_f64! for BOTH NativeFloat64 and NativeFloat32,
    #     mie_coefficients_kernel_ds! for DSEmulated -- see _mie_kernel1) ---
    kernel1 = _mie_kernel1(precision_policy, backend)
    kernel1(an_dev, bn_dev, x_dev, m_ref_re, m_ref_im, nmax_dev, nmx_max; ndrange=nquad_radius)
    KernelAbstractions.synchronize(backend)

    return (FT=FT, RA=RA, r=r, w=w, k_wavenum=k_wavenum, n_max_global=n_max_global,
            nquad_radius=nquad_radius, x_dev=x_dev, nmax_dev=nmax_dev, an_dev=an_dev, bn_dev=bn_dev)
end

@doc raw"""
    compute_aerosol_optical_properties_nodes_gpu(radii, weights, n_real, n_imag,
        wavelength, backend; precision_policy = NativeFloat64(), l_max = nothing)
        -> AerosolOptics

GPU/KernelAbstractions implementation backing
[`compute_aerosol_optical_properties_nodes`](@ref) — takes an explicit KA
`backend` (e.g. `KernelAbstractions.CPU()` or `CUDA.CUDABackend()`), mirroring
[`compute_aerosol_optical_properties_gpu`](@ref)'s calling convention so it can
be exercised directly in tests on a CPU backend without a real GPU device.

Shares its preamble (validation, `FT`/`RA` derivation, per-node arrays, Kernel
1) with [`compute_aerosol_extinction_nodes_gpu`](@ref) via the internal
`_prepare_node_mie_gpu` helper. Reuses the existing per-node kernels from
`gpu_mie_kernels.jl` (`mie_coefficients_kernel_f64!`/`_ds!`,
`amplitude_phase_kernel!`, `cross_sections_kernel!`) UNCHANGED. The weighted
bulk reduction across nodes and the Greek-coefficient angular projection
additionally launch the previously-unused device-resident kernels
`weighted_sum_kernel!`, `size_reduction_kernel!`, and
`greek_coefficients_kernel!` — so, unlike
[`compute_aerosol_optical_properties_gpu`](@ref) (which copies the per-node
phase-matrix/cross-section arrays back to the host and reduces there), this
path keeps the reduction device-resident and copies back only the six bulk
Greek arrays plus the two bulk scalars (`ω̃`, `k`).

**Precision**: the per-node kernels (1-3) run in the kernel's native float
type `FT = eltype(radii)` (Float32 under `DSEmulated`/`NativeFloat32`,
matching the existing GPU path's speed rationale). The reduction/Greek
kernels (4-5) accumulate in `RA = _mie_reduction_type(precision_policy, FT)`
— `Float64` for `NativeFloat64` always, `Float64` for `DSEmulated` when
`FT === Float32` (Float64-widened + Neumaier-compensated, exactly like
`compute_aerosol_optical_properties_gpu`'s host-side reduction — except here
the widened reduction runs as device arrays/kernels rather than after an
`Array()` copy-back), and — new — `Float32` (never widened) for
`NativeFloat32`, so the entire pipeline allocates zero Float64 device arrays
under that policy. See [`MiePrecisionPolicy`](@ref) for the full table.
"""
function compute_aerosol_optical_properties_nodes_gpu(
    radii::AbstractVector, weights::AbstractVector,
    n_real::Real, n_imag::Real, wavelength::Real, backend;
    precision_policy::MiePrecisionPolicy = NativeFloat64(),
    l_max::Union{Nothing,Integer} = nothing,
)
    (; FT, RA, r, w, k_wavenum, n_max_global, nquad_radius, x_dev, nmax_dev, an_dev, bn_dev) =
        _prepare_node_mie_gpu(radii, weights, n_real, n_imag, wavelength, backend, precision_policy)

    n_mu = 2n_max_global - 1
    μ, w_μ = gausslegendre(n_mu)
    leg_π, leg_τ = compute_mie_π_τ(FT.(μ), n_max_global)

    AT = KernelAbstractions.allocate

    # --- Remaining per-node arrays: native kernel precision FT (speed) ---
    f11_dev     = AT(backend, FT, (n_mu, nquad_radius))
    f33_dev     = AT(backend, FT, (n_mu, nquad_radius))
    f12_dev     = AT(backend, FT, (n_mu, nquad_radius))
    f34_dev     = AT(backend, FT, (n_mu, nquad_radius))
    C_sca_dev   = AT(backend, FT, nquad_radius)
    C_ext_dev   = AT(backend, FT, nquad_radius)
    leg_pi_dev  = AT(backend, FT, (n_mu, n_max_global))
    leg_tau_dev = AT(backend, FT, (n_mu, n_max_global))

    KernelAbstractions.copyto!(backend, leg_pi_dev, FT.(leg_π))
    KernelAbstractions.copyto!(backend, leg_tau_dev, FT.(leg_τ))

    # --- Kernel 2+3: amplitude functions + phase matrix (unchanged, existing kernel) ---
    kernel23 = amplitude_phase_kernel!(backend)
    kernel23(f11_dev, f33_dev, f12_dev, f34_dev,
             an_dev, bn_dev, leg_pi_dev, leg_tau_dev,
             x_dev, nmax_dev; ndrange=(n_mu, nquad_radius))
    KernelAbstractions.synchronize(backend)

    # --- Kernel 4a: per-node cross-sections (unchanged, existing kernel) ---
    kernel4a = cross_sections_kernel!(backend)
    kernel4a(C_sca_dev, C_ext_dev, an_dev, bn_dev, k_wavenum, nmax_dev; ndrange=nquad_radius)
    KernelAbstractions.synchronize(backend)

    # --- Device-resident weighted reduction across nodes (RA precision) ---
    # These kernels (weighted_sum_kernel!, size_reduction_kernel!) already existed
    # in gpu_mie_kernels.jl but had no call sites before this API.
    w_dev  = AT(backend, RA, nquad_radius); KernelAbstractions.copyto!(backend, w_dev, RA.(w))
    wr_dev = AT(backend, RA, nquad_radius); KernelAbstractions.copyto!(backend, wr_dev, RA.(4π .* r.^2 .* w))

    bulk_Csca_dev = AT(backend, RA, 1)
    bulk_Cext_dev = AT(backend, RA, 1)
    wsum = weighted_sum_kernel!(backend)
    wsum(bulk_Csca_dev, C_sca_dev, w_dev; ndrange=1)
    wsum(bulk_Cext_dev, C_ext_dev, w_dev; ndrange=1)
    KernelAbstractions.synchronize(backend)

    bulk_f11_dev = AT(backend, RA, n_mu)
    bulk_f33_dev = AT(backend, RA, n_mu)
    bulk_f12_dev = AT(backend, RA, n_mu)
    bulk_f34_dev = AT(backend, RA, n_mu)
    sred = size_reduction_kernel!(backend)
    sred(bulk_f11_dev, f11_dev, wr_dev; ndrange=n_mu)
    sred(bulk_f33_dev, f33_dev, wr_dev; ndrange=n_mu)
    sred(bulk_f12_dev, f12_dev, wr_dev; ndrange=n_mu)
    sred(bulk_f34_dev, f34_dev, wr_dev; ndrange=n_mu)
    KernelAbstractions.synchronize(backend)

    # Normalize the bulk phase matrix by the bulk scattering cross section.
    # A single scalar is pulled back to host to drive the broadcast (negligible
    # cost; the phase-matrix ARRAYS themselves stay device-resident).
    bulk_C_sca_host = Array(bulk_Csca_dev)[1]
    bulk_C_ext_host = Array(bulk_Cext_dev)[1]
    inv_bulk_C_sca = one(RA) / bulk_C_sca_host
    bulk_f11_dev .*= inv_bulk_C_sca
    bulk_f33_dev .*= inv_bulk_C_sca
    bulk_f12_dev .*= inv_bulk_C_sca
    bulk_f34_dev .*= inv_bulk_C_sca

    # --- Device-resident Greek-coefficient projection (RA precision) ---
    P, P², R², T² = compute_legendre_poly(RA.(μ), n_mu)
    P_dev  = AT(backend, RA, (n_mu, n_mu)); KernelAbstractions.copyto!(backend, P_dev,  RA.(P))
    P2_dev = AT(backend, RA, (n_mu, n_mu)); KernelAbstractions.copyto!(backend, P2_dev, RA.(P²))
    R2_dev = AT(backend, RA, (n_mu, n_mu)); KernelAbstractions.copyto!(backend, R2_dev, RA.(R²))
    T2_dev = AT(backend, RA, (n_mu, n_mu)); KernelAbstractions.copyto!(backend, T2_dev, RA.(T²))
    wmu_dev = AT(backend, RA, n_mu); KernelAbstractions.copyto!(backend, wmu_dev, RA.(w_μ))

    α_dev = AT(backend, RA, n_mu); β_dev = AT(backend, RA, n_mu)
    γ_dev = AT(backend, RA, n_mu); δ_dev = AT(backend, RA, n_mu)
    ϵ_dev = AT(backend, RA, n_mu); ζ_dev = AT(backend, RA, n_mu)

    greek_kernel = greek_coefficients_kernel!(backend)
    greek_kernel(α_dev, β_dev, γ_dev, δ_dev, ϵ_dev, ζ_dev,
                 bulk_f11_dev, bulk_f33_dev, bulk_f12_dev, bulk_f34_dev,
                 P_dev, P2_dev, R2_dev, T2_dev, wmu_dev, n_mu; ndrange=n_mu)
    KernelAbstractions.synchronize(backend)

    # --- Copy back ONLY bulk quantities: 6 Greek vectors + 2 scalars ---
    α = Array(α_dev); β = Array(β_dev); γ = Array(γ_dev)
    δ = Array(δ_dev); ϵ = Array(ϵ_dev); ζ = Array(ζ_dev)

    if RA !== FT
        greek_coefs = GreekCoefs(convert.(FT, α), convert.(FT, β), convert.(FT, γ),
                                 convert.(FT, δ), convert.(FT, ϵ), convert.(FT, ζ))
        raw = AerosolOptics(greek_coefs=greek_coefs,
                            ω̃=FT(bulk_C_sca_host / bulk_C_ext_host),
                            k=FT(bulk_C_ext_host), fᵗ=FT(1))
    else
        greek_coefs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
        raw = AerosolOptics(greek_coefs=greek_coefs,
                            ω̃=(bulk_C_sca_host / bulk_C_ext_host),
                            k=bulk_C_ext_host, fᵗ=one(FT))
    end

    return l_max === nothing ? raw : _slice_greek(raw, l_max)
end

#=
=====================================================================
EXTINCTION-ONLY fast path.
=====================================================================
Motivation: a downstream hybrid consumer (GCHPIO's `:parametric_exact_tau`
mode) needs only the number-mean extinction cross-section from ONE node set
(e.g. a per-layer sectional/tabulated size distribution) while taking its
phase-matrix SHAPE (Greek coefficients) from a DIFFERENT node set (e.g. a
representative parametric ensemble). Before this existed, that consumer had
to call the full `compute_aerosol_optical_properties_nodes` API and
immediately discard the amplitude matrices (S₁/S₂) and all six Greek series —
work that dominates the cost at large angular orders (`n_mu = 2·n_max - 1`
scales with the size parameter of the largest node, and the Greek projection
is `O(n_mu × l_max)`), while `k` itself only needs the per-node Mie
coefficients (aₙ, bₙ) and a single weighted sum — no angular quadrature at
all.
=====================================================================
=#

# Sibling warn-once state for `compute_aerosol_extinction_nodes`, kept
# independent of `_WARNED_NO_GPU_MIE_NODES` (same rationale as that Ref's own
# comment: independent state so tests for one entry point don't clear the
# other's already-warned flag).
const _WARNED_NO_GPU_MIE_EXTINCTION_NODES = Ref(false)

@doc raw"""
    compute_aerosol_extinction_nodes(radii, weights, n_real, n_imag, wavelength;
        architecture = Architectures.CPU(),
        precision_policy = nothing) -> k

EXTINCTION-ONLY fast path for the caller-node Mie seam: returns just the
number-mean bulk extinction cross-section `k` **[μm²]** — the SAME quantity,
same convention, and the SAME arithmetic as `.k` from
[`compute_aerosol_optical_properties_nodes`](@ref) — without computing any
amplitude matrices, phase-matrix elements, or the six Greek coefficient
series. This equality is EXACT (not merely close) on the CPU path, on CUDA
(real hardware and the KA-CPU KernelAbstractions backend), and on the
KA-CPU backend under every precision policy — see "Numerics" below for
exactly what "exact" means on each path and how it is tested.

# Motivation
A downstream hybrid consumer (GCHPIO's `:parametric_exact_tau` mode) needs
only the number-mean extinction cross-section from ONE node set while taking
shapes (Greek coefficients) from another. Today it would have to call the
full [`compute_aerosol_optical_properties_nodes`](@ref) API and discard the
amplitude matrices and all six Greek series, which dominate the cost at large
angular orders. This function computes only what `k` actually depends on: the
per-node Mie coefficients (aₙ, bₙ) and a weighted sum over nodes — no
Gauss-Legendre angular quadrature, π/τ functions, amplitude functions, or
Greek-coefficient projection anywhere in this path.

# Arguments
Same node-ensemble arguments as
[`compute_aerosol_optical_properties_nodes`](@ref), EXCEPT there is no
`polarization`/`truncation`/`l_max` — none of these are relevant to a bare
extinction cross-section (there is no phase matrix or Greek series to
polarize, truncate, or length-cap).
- `radii::AbstractVector`: wet particle radii, one per node, **[μm]**.
- `weights::AbstractVector`: non-negative node weights (number-mixing-ratio or
  number-concentration weights). **Need not be normalized** — divided by
  their sum internally, so `k` is invariant to a uniform rescaling of
  `weights`.
- `n_real, n_imag`: the ONE complex refractive index shared by the whole node
  ensemble, convention `n = n_real - i·n_imag` (`n_imag ≥ 0`).
- `wavelength`: same units as `radii` (μm internally).

# Keyword arguments
- `architecture`: `Architectures.CPU()` (default) or a GPU architecture with a
  registered Mie pipeline (`Architectures.has_gpu_mie`) — identical routing to
  [`compute_aerosol_optical_properties_nodes`](@ref) (CPU / GPU / warn-once
  CPU fallback for an architecture without a GPU Mie pipeline).
- `precision_policy`: GPU-only, see [`MiePrecisionPolicy`](@ref)
  (`NativeFloat64`/`DSEmulated`/`NativeFloat32`); ignored on the CPU path.

# Returns
The number-mean bulk extinction cross-section `k` **[μm²]**, as a scalar of
the promoted float type of the inputs.

# Validation
Identical to [`compute_aerosol_optical_properties_nodes`](@ref):
`_validate_node_inputs` plus `sum(weights) > 0`.

# Numerics
**CPU path**: internal computation runs in `IC = Float64` (same IC/FT
narrowing convention as the rest of this module). Per node, this computes
exactly the SAME `C_ext` expression the shared NAI-2 core uses —
`2π/k² * dot(n_v, real.(an_v .+ bn_v))` (the internal `_nai2_bulk_optics`
core's `C_ext` formula, via the shared `_mie_node_ab_and_C_ext!` helper) — in
the SAME arithmetic order, weighted by the SAME normalized weights (same
normalization step), so `k` is BIT-IDENTICAL (`===`) to
`compute_aerosol_optical_properties_nodes(...).k` given identical
`(radii, weights, n_real, n_imag, wavelength)`.

**GPU path** ([`compute_aerosol_extinction_nodes_gpu`](@ref)): shares its
preamble (validation, `FT`/`RA`, per-node arrays, Kernel 1) with
[`compute_aerosol_optical_properties_nodes_gpu`](@ref) via the internal
`_prepare_node_mie_gpu` helper, so both entry points run `cross_sections_kernel!`
and the single-thread (`ndrange=1`) `weighted_sum_kernel!` reduction on
IDENTICAL per-node/weight device arrays. `k` is confirmed EXACT-EQUALITY
(`==`, not merely `isapprox`) tested for all three precision policies
(`NativeFloat64`, `DSEmulated`, `NativeFloat32`) on both the standard and
wide-range TOMAS-like node sets — CPU path (`===`), the KA-CPU
KernelAbstractions backend (`==`), and real CUDA hardware (`==`), all in
`test/test_mie_nodes.jl`. On `MetalGPU()` (`NativeFloat32` only, see
[`MiePrecisionPolicy`](@ref)) the same exact-equality test is wired in and
expected to hold by the same reasoning, but is unverified on real Metal
hardware (this package's CI has none); the test is written to activate
automatically the first time it runs on a Mac with a functional Metal device.
"""
function compute_aerosol_extinction_nodes(
    radii::AbstractVector, weights::AbstractVector,
    n_real::Real, n_imag::Real,
    wavelength::Real;
    architecture::Architectures.AbstractArchitecture = Architectures.CPU(),
    precision_policy::Union{Nothing,MiePrecisionPolicy} = nothing,
)
    _validate_node_inputs(radii, weights, n_imag)
    @assert sum(weights) > 0 "weights must not sum to zero"

    # Output float type: promoted user type of all numeric inputs (mirrors
    # compute_aerosol_optical_properties_nodes's FT convention).
    FT = float(promote_type(eltype(radii), eltype(weights),
                            typeof(n_real), typeof(n_imag), typeof(wavelength)))

    if architecture isa Architectures.CPU
        return _aerosol_extinction_nodes_cpu(radii, weights, n_real, n_imag, wavelength, FT)
    elseif Architectures.has_gpu_mie(architecture)
        backend = Architectures.ka_backend(architecture)
        policy = precision_policy === nothing ?
                 Architectures.default_mie_precision_policy(architecture, FT) :
                 precision_policy
        return compute_aerosol_extinction_nodes_gpu(radii, weights, n_real, n_imag,
                                                     wavelength, backend;
                                                     precision_policy = policy)
    else
        if !_WARNED_NO_GPU_MIE_EXTINCTION_NODES[]
            @warn "no GPU Mie pipeline for $(architecture); computing caller-node Mie extinction on CPU"
            _WARNED_NO_GPU_MIE_EXTINCTION_NODES[] = true
        end
        return _aerosol_extinction_nodes_cpu(radii, weights, n_real, n_imag, wavelength, FT)
    end
end

"""
    _aerosol_extinction_nodes_cpu(radii, weights, n_real, n_imag, λ, FT) -> k

CPU implementation backing [`compute_aerosol_extinction_nodes`](@ref).
Normalizes `weights` and derives `n_max` from the actual node radii (no cap),
exactly as `_aerosol_optical_properties_nodes_cpu` does, but computes ONLY the
per-node extinction cross-section (via the [`_nai2_bulk_optics`](@ref)-shared
helper `_mie_node_ab_and_C_ext!`, defined in `compute_NAI2.jl`) and the
weighted bulk sum — no scattering cross-section, amplitude, or phase-matrix
arrays are ever allocated. Same IC=Float64-internal convention; narrowed to
`FT` only at the end.
"""
function _aerosol_extinction_nodes_cpu(radii, weights, n_real, n_imag, λ, FT::Type)
    IC = FT <: AbstractFloat ? Float64 : FT

    r = IC.(radii)
    w = IC.(weights)          # broadcast always materializes a fresh copy, so
    w ./= sum(w)              # in-place normalization never aliases the input

    k = IC(2π) / IC(λ)
    x_size_param = k .* r
    n_max = get_n_max(maximum(x_size_param))

    nquad_radius = length(r)

    # m_ref, nmx_max, and the pre-allocated buffers below mirror
    # _nai2_bulk_optics exactly (same formulas, same order) so that C_ext[i]
    # comes out bit-identical.
    m_ref = IC(n_real) - IC(n_imag) * im
    y_max = maximum(x_size_param) * abs(m_ref)
    nmx_max = round(Int, max(n_max, y_max) + 51)
    an  = zeros(Complex{IC}, n_max)
    bn  = zeros(Complex{IC}, n_max)
    Dₙ  = zeros(Complex{IC}, nmx_max)
    n_  = IC.(2 .* collect(1:n_max) .+ 1)

    C_ext = zeros(IC, nquad_radius)

    for i = 1:nquad_radius

        n_max_i = get_n_max(x_size_param[i])
        @assert n_max_i ≤ n_max "supplied n_max=$n_max is too small for radius index $i (needs $n_max_i); n_max must come from the largest size parameter in the node set"

        # Same shared helper _nai2_bulk_optics uses for its own an_v/bn_v/C_ext_i
        # (see compute_NAI2.jl) -- literally the same FLOP sequence, not just a
        # textual copy, which is what makes k bit-identical (===) to
        # compute_aerosol_optical_properties_nodes(...).k.
        _, _, C_ext_i = _mie_node_ab_and_C_ext!(an, bn, Dₙ, n_, x_size_param[i], n_max_i, m_ref, k)
        @inbounds C_ext[i] = C_ext_i
    end

    bulk_C_ext = sum(w .* C_ext)

    return IC <: AbstractFloat && FT <: AbstractFloat ? FT(bulk_C_ext) : bulk_C_ext
end

@doc raw"""
    compute_aerosol_extinction_nodes_gpu(radii, weights, n_real, n_imag,
        wavelength, backend; precision_policy = NativeFloat64()) -> k

GPU/KernelAbstractions implementation backing
[`compute_aerosol_extinction_nodes`](@ref) — takes an explicit KA `backend`
(e.g. `KernelAbstractions.CPU()` or `CUDA.CUDABackend()`), mirroring
[`compute_aerosol_optical_properties_nodes_gpu`](@ref)'s calling convention so
it can be exercised directly in tests on a CPU backend without a real GPU
device.

Shares its preamble (validation, `FT`/`RA` derivation, per-node arrays, Kernel
1) with [`compute_aerosol_optical_properties_nodes_gpu`](@ref) via the
internal `_prepare_node_mie_gpu` helper. Launches ONLY the per-node
`cross_sections_kernel!` on top of that shared preamble (which computes
`C_ext` AND `C_sca` together in one fused pass — cheaper to keep fused than to
add a `C_ext`-only kernel variant; the unused `C_sca` output is simply not
reduced), then the device-resident `weighted_sum_kernel!` reduction across
nodes. UNLIKE [`compute_aerosol_optical_properties_nodes_gpu`](@ref), this path
launches no `amplitude_phase_kernel!`, `size_reduction_kernel!`, or
`greek_coefficients_kernel!`, and needs no angular grid at all (no
`gausslegendre`, `compute_mie_π_τ`, or `compute_legendre_poly` calls), since
none of them affect `k`. Only the single scalar `bulk_C_ext` is copied back to
the host.

**Precision**: identical convention to
[`compute_aerosol_optical_properties_nodes_gpu`](@ref) — the per-node kernels
run in the kernel's native float type `FT = eltype(radii)` (Float32 under
`DSEmulated`/`NativeFloat32`); the reduction accumulates in
`RA = _mie_reduction_type(precision_policy, FT)`, narrowed back to `FT` only
for the returned scalar. Under `NativeFloat32`, `RA === FT === Float32`, so
this function allocates ZERO Float64 device arrays end to end.
"""
function compute_aerosol_extinction_nodes_gpu(
    radii::AbstractVector, weights::AbstractVector,
    n_real::Real, n_imag::Real, wavelength::Real, backend;
    precision_policy::MiePrecisionPolicy = NativeFloat64(),
)
    (; FT, RA, w, k_wavenum, nquad_radius, nmax_dev, an_dev, bn_dev) =
        _prepare_node_mie_gpu(radii, weights, n_real, n_imag, wavelength, backend, precision_policy)

    AT = KernelAbstractions.allocate

    C_sca_dev = AT(backend, FT, nquad_radius)  # required output slot of the fused kernel; never reduced
    C_ext_dev = AT(backend, FT, nquad_radius)

    # --- Kernel: per-node cross-sections (unchanged, existing kernel) ---
    kernel4a = cross_sections_kernel!(backend)
    kernel4a(C_sca_dev, C_ext_dev, an_dev, bn_dev, k_wavenum, nmax_dev; ndrange=nquad_radius)
    KernelAbstractions.synchronize(backend)

    # --- Device-resident weighted reduction across nodes (RA precision) ---
    w_dev = AT(backend, RA, nquad_radius)
    KernelAbstractions.copyto!(backend, w_dev, RA.(w))
    bulk_Cext_dev = AT(backend, RA, 1)
    wsum = weighted_sum_kernel!(backend)
    wsum(bulk_Cext_dev, C_ext_dev, w_dev; ndrange=1)
    KernelAbstractions.synchronize(backend)

    # --- Copy back ONLY the one bulk scalar ---
    bulk_C_ext_host = Array(bulk_Cext_dev)[1]

    return RA !== FT ? FT(bulk_C_ext_host) : bulk_C_ext_host
end
