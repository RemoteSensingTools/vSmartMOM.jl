#=
 
This file wraps `compute_aerosol_optical_properties` so that autodiff users and non-autodiff 
users can call the same function with just a keyword argument change. 

=#

"""
    convert_jacobian_result_to_aerosol_optics(result) -> AerosolOptics

Convert a `DiffResults.JacobianResult` returned by ForwardDiff into an
[`AerosolOptics`](@ref) object.

The flattened value vector is interpreted as:
`[α; β; γ; δ; ϵ; ζ; ω̃; k]`, and `result.derivs[1]` is stored in
`AerosolOptics.derivs`.
"""
function convert_jacobian_result_to_aerosol_optics(result)
    
    value = result.value
    derivs = result.derivs[1]

    greek_length = Int64((length(value) - 2)/6)

    α = value[1:greek_length]
    β = value[greek_length + 1 : 2greek_length]
    γ = value[2greek_length + 1 : 3greek_length]
    δ = value[3greek_length + 1 : 4greek_length]
    ϵ = value[4greek_length + 1 : 5greek_length]
    ζ = value[5greek_length + 1 : 6greek_length]

    greek_coefs = GreekCoefs(α, β, γ, δ, ϵ, ζ)

    ω̃ = value[end-1]
    k = value[end]

    return AerosolOptics(greek_coefs=greek_coefs, ω̃=ω̃, k=k, derivs=derivs, fᵗ=eltype(ω̃)(1)) 
end

@doc raw"""
    compute_aerosol_optical_properties(model::MieModel; autodiff=false) -> AerosolOptics

Single-verb entry point for aerosol optics. Everything (Fourier-decomposition
method, compute architecture, GPU precision policy) is read off the model.

- `autodiff=false` (default): dispatches on `model.computation_type` **and**
  `model.architecture`:
    * `NAI2` + `CPU()` → the analytic CPU implementation
      ([`compute_aerosol_optical_properties`](@ref) 2-arg method).
    * `NAI2` + a GPU architecture with a Mie pipeline (`has_gpu_mie` true,
      i.e. CUDA `GPU()`) → the KernelAbstractions GPU pipeline
      ([`compute_aerosol_optical_properties_gpu`](@ref)).
    * `NAI2`/`PCW` + a non-CPU architecture **without** a GPU Mie pipeline
      (e.g. `MetalGPU()`) → one-time `@warn` then CPU fallback (Mie on CPU; RT
      arrays still run on that architecture).
    * `PCW` + GPU → one-time `@warn` then CPU fallback (no GPU PCW kernel).
    * `PCW` + `CPU()` → the analytic PCW implementation.
- `autodiff=true`: computes the Jacobian with respect to the 4 aerosol
  parameters ``\mathbf{x}=[r_m,\sigma,n_r,n_i]`` using ForwardDiff. AD always
  runs on the CPU analytic kernel (Dual numbers do not flow through the GPU
  kernels).

The AD Jacobian is stored in `AerosolOptics.derivs` with shape
`(6L + 2, 4)`, where `L` is the Greek coefficient length and rows are stacked as
`[α; β; γ; δ; ϵ; ζ; ω̃; k]`.
"""
function compute_aerosol_optical_properties(model::MieModel ; autodiff=false)

    # Closure that ForwardDiff differentiates: x = [r_m, σ_g, nᵣ, nᵢ] as Dual numbers.
    # The Dual numbers must propagate through the entire Mie computation.
    function compute_aerosol_optical_properties_autodiff(x)

        (; computation_type, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_A, wigner_B) = model

        # Construct aerosol with Dual-typed parameters so ForwardDiff can track derivatives
        aerosol_x = Aerosol(LogNormal(log(x[1]), log(x[2])), x[3], x[4])
        # The autodiff/Jacobian path always runs on the CPU analytic kernel:
        # ForwardDiff Dual numbers do not flow through the GPU KernelAbstractions
        # kernels, so we force CPU() here (and precision_policy is irrelevant on
        # the CPU path). The full 11-field positional constructor is used.
        model_x = MieModel(computation_type, aerosol_x, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_A, wigner_B, Architectures.CPU(), nothing)

        aerosol_optics = compute_aerosol_optical_properties(model_x)
    
        return [aerosol_optics.greek_coefs.α; 
                aerosol_optics.greek_coefs.β; 
                aerosol_optics.greek_coefs.γ; 
                aerosol_optics.greek_coefs.δ; 
                aerosol_optics.greek_coefs.ϵ; 
                aerosol_optics.greek_coefs.ζ;
                aerosol_optics.ω̃;
                aerosol_optics.k]
    end

    if (autodiff)

        x = [exp(model.aerosol.size_distribution.μ), 
            exp(model.aerosol.size_distribution.σ), 
            model.aerosol.nᵣ, 
            model.aerosol.nᵢ]

        # Get length of greek coefs
        r, wᵣ = gauleg(model.nquad_radius, 0.0, model.r_max ; norm=true)
        N_max = get_n_max(2 * π * model.r_max/ model.λ)
        greek_length = 2 * N_max - 1

        result = DiffResults.JacobianResult(zeros(6 * greek_length + 2), x)
        ForwardDiff.jacobian!(result, compute_aerosol_optical_properties_autodiff, x);
        return convert_jacobian_result_to_aerosol_optics(result);

    else
        return _dispatch_aerosol_optics(model.architecture, model)
    end

end

# ============================================================================
# Architecture / method router for the single-verb entry point
# ============================================================================
#
# `compute_aerosol_optical_properties(model)` (1-arg) lands in the method above
# and, with autodiff=false, calls this router. Dispatch is on the architecture
# object so the weak GPU dependency boundary is respected.
#
# The architecture is split into three capability tiers, NOT three concrete
# types, so new backends slot in without touching this file:
#   - CPU()                                  → analytic CPU path.
#   - non-CPU WITH    a GPU Mie pipeline      → GPU KernelAbstractions pipeline.
#     (`has_gpu_mie(arch)` true, e.g. CUDA `GPU()` via `vSmartMOMCUDAExt`.)
#   - non-CPU WITHOUT a GPU Mie pipeline      → warn-once + CPU fallback.
#     (`has_gpu_mie(arch)` false, e.g. `MetalGPU()`: Metal has no Mie kernels,
#      so Mie runs on the CPU exactly as it did before this router existed; the
#      RT arrays are unaffected and still run on Metal.)
#
# The output float type `FT2` rides on the model's `FT` type parameter
# (`_mie_output_type`), and the auto GPU precision policy is FT-aware
# (`default_mie_precision_policy(arch, FT)`), so a Float32 model auto-selects
# the Float32-native DSEmulated path instead of the Float64-only NativeFloat64.

"""
    _dispatch_aerosol_optics(arch, model) -> AerosolOptics

Internal router selecting the CPU vs GPU Mie path from `model.architecture`
and `model.computation_type`. Dispatching on both arguments keeps all four
(arch × method) combinations unambiguous. Not part of the public API.

The output float type `FT2` is taken from the model's `FT` type parameter
(`_mie_output_type`), so `compute_aerosol_optical_properties(model)` is fully
self-describing — the output precision rides on the `MieModel` exactly like the
architecture does.
"""
@inline _dispatch_aerosol_optics(::Architectures.CPU, model::MieModel{<:NAI2}) =
    compute_aerosol_optical_properties(model, _mie_output_type(model))

@inline _dispatch_aerosol_optics(::Architectures.CPU, model::MieModel{<:PCW}) =
    compute_aerosol_optical_properties(model, _mie_output_type(model))

# One-time warning state for the "non-CPU architecture without a GPU Mie
# pipeline → compute on CPU" fallback (e.g. Metal), keyed nothing-fancy: a single
# flag is enough because the message names the architecture and the situation is
# static for a given run.
const _WARNED_NO_GPU_MIE = Ref(false)

# Helper: emit the CPU-fallback warning at most once, then run the CPU path.
@inline function _cpu_mie_fallback(arch, model)
    if !_WARNED_NO_GPU_MIE[]
        @warn "no GPU Mie pipeline for $(arch); computing Mie on CPU, RT arrays unaffected"
        _WARNED_NO_GPU_MIE[] = true
    end
    return compute_aerosol_optical_properties(model, _mie_output_type(model))
end

# Non-CPU + NAI2. Branch on the capability trait, not the concrete type:
#   has_gpu_mie(arch) true  → GPU KernelAbstractions pipeline.
#   has_gpu_mie(arch) false → warn-once + CPU compute (restores old Metal path).
function _dispatch_aerosol_optics(arch::Architectures.AbstractArchitecture,
                                  model::MieModel{<:NAI2})
    if !Architectures.has_gpu_mie(arch)
        return _cpu_mie_fallback(arch, model)
    end
    backend = Architectures.ka_backend(arch)
    # Two distinct precision axes — must not be confused:
    #   FT_kernel : the GPU kernel's working float type, derived from the aerosol
    #               refractive-index element type (`eltype(nᵣ)`).  This is EXACTLY
    #               how `compute_aerosol_optical_properties_gpu` derives its `FT`
    #               (line `FT = eltype(nᵣ)` in compute_NAI2_gpu.jl).  The
    #               NativeFloat64 policy asserts FT_kernel === Float64; passing a
    #               Float32-kernel model there trips that assert.
    #   FT2       : the OUTPUT float type for the returned AerosolOptics.  It rides
    #               on the model's `FT` type parameter (set by `r_max`'s type) and
    #               is orthogonal to the kernel precision.  A Float32 aerosol with
    #               a Float64 r_max has FT2 = Float64 but FT_kernel = Float32.
    FT_kernel = _mie_kernel_type(model)   # governs policy selection
    FT2       = _mie_output_type(model)   # governs output conversion
    policy  = model.precision_policy === nothing ?
              Architectures.default_mie_precision_policy(arch, FT_kernel) :
              model.precision_policy
    return compute_aerosol_optical_properties_gpu(model, backend;
                                                  precision_policy = policy,
                                                  FT2 = FT2)
end

# One-time warning state for the PCW+GPU fallback (no GPU PCW kernel exists at
# all, independent of `has_gpu_mie`).
const _WARNED_PCW_GPU = Ref(false)

# Non-CPU + PCW → no GPU kernel; warn once and fall back to the CPU path. This
# covers both the "GPU has a NAI2 Mie pipeline but no PCW kernel" case and the
# "architecture has no GPU Mie pipeline at all" (Metal) case — either way PCW
# runs on the CPU.
function _dispatch_aerosol_optics(arch::Architectures.AbstractArchitecture,
                                  model::MieModel{<:PCW})
    if Architectures.has_gpu_mie(arch)
        if !_WARNED_PCW_GPU[]
            @warn "No GPU PCW Mie kernel exists; falling back to the CPU PCW path " *
                  "for this and all subsequent GPU PCW models." architecture = arch
            _WARNED_PCW_GPU[] = true
        end
        return compute_aerosol_optical_properties(model, _mie_output_type(model))
    else
        # No GPU Mie pipeline for this architecture at all (e.g. Metal): route
        # through the same warn-once CPU fallback used by the NAI2 branch.
        return _cpu_mie_fallback(arch, model)
    end
end

# Output float type carried by the model's `FT` type parameter. Falls back to
# Float64 when `FT` is not a concrete AbstractFloat (e.g. a ForwardDiff Dual or
# the autodiff closure's transient model) — those paths never reach this router.
@inline _mie_output_type(::MieModel{FDT,FT}) where {FDT,FT} =
    FT <: AbstractFloat ? FT : Float64

# Kernel input float type: the precision at which the GPU Mie kernels actually
# run. This mirrors EXACTLY the first thing `compute_aerosol_optical_properties_gpu`
# does: `FT = eltype(nᵣ)` (see compute_NAI2_gpu.jl). The GPU precision policy
# must be selected on this type, NOT on the output type, because the NativeFloat64
# kernel asserts `FT === Float64`. A model with Float32 aerosol parameters but a
# Float64 r_max has _mie_output_type → Float64 but _mie_kernel_type → Float32,
# and needs DSEmulated (the Float32-native policy), not NativeFloat64.
@inline _mie_kernel_type(model::MieModel) = eltype(model.aerosol.nᵣ)
