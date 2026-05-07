using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using vSmartMOM.Architectures
using Statistics
using LinearAlgebra

# Cap BLAS threads. Our RT matrices are 48×48 (Stokes_IQUV NSTREAMS=8) or
# smaller; 128-thread BLAS coordination overhead dominates the work at this
# size. 8 threads per process is empirically the sweet spot — leaves room for
# multiple processes when running specs in parallel via `runtests_parallel.sh`.
BLAS_THREADS_CAP = parse(Int, get(ENV, "VLIDORT_BASELINE_BLAS_THREADS", "8"))
LinearAlgebra.BLAS.set_num_threads(BLAS_THREADS_CAP)

const VLIDORT_BASELINE_DIR = @__DIR__
const REF_DIR     = joinpath(VLIDORT_BASELINE_DIR, "reference_data")
const CONFIG_DIR  = joinpath(VLIDORT_BASELINE_DIR, "configs")

# CUDA is a weak dep — load if available; failures are silent.
const CUDA_LOADED = try
    @eval using CUDA
    true
catch
    false
end
gpu_available() = CUDA_LOADED ? CUDA.functional() : false

"""
    axis_specs(; include_gpu = gpu_available())

Cross-product of (quadrature_type, float_type, architecture) specs we want
each baseline case to run under. Returns NamedTuples; the last field is the
architecture constructor (call it: `arch_ctor()`). Tolerance scaling for
F32 rows lives in `tol_scale(spec)`.
"""
function axis_specs(; include_gpu::Bool = gpu_available())
    # Single-spec filter for parallel-process invocation: when env var
    # `VLIDORT_BASELINE_ONLY_SPEC` is set (e.g.
    # "GaussLegQuad/Float64/CPU"), `axis_specs()` returns only that
    # spec. The parallel runner script `runtests_parallel.sh` uses this to
    # fan out across processes.
    only = get(ENV, "VLIDORT_BASELINE_ONLY_SPEC", "")
    # VLIDORT's `NSTREAMS` are Gauss-Legendre half-space streams. Keep this
    # baseline suite on the matching vSmartMOM quadrature; Radau has a
    # different weighted-node count and is a convergence experiment, not an
    # apples-to-apples VLIDORT comparison.
    quads  = (:GaussLegQuad,)
    floats = (Float64, Float32)
    archs  = include_gpu ? ((:CPU, () -> Architectures.CPU()),
                            (:GPU, () -> Architectures.GPU())) :
                           ((:CPU, () -> Architectures.CPU()),)
    specs = NamedTuple{(:quad, :float, :arch_name, :arch_ctor),
                       Tuple{Symbol,DataType,Symbol,Function}}[]
    for q in quads, f in floats, (an, ac) in archs
        s = (quad = q, float = f, arch_name = an, arch_ctor = ac)
        (only == "" || spec_tag(s) == only) && push!(specs, s)
    end
    return specs
end

"Tag string suitable for use as a `@testset` label."
spec_tag(s) = "$(s.quad)/$(s.float)/$(s.arch_name)"

"Multiplicative tolerance scale: F32 gets 1e3×, GPU gets 1× (BLAS-equivalent)."
tol_scale(s) = s.float === Float32 ? 1.0e3 : 1.0

"""
    apply_overrides!(params, spec)

Mutate a `vSmartMOM_Parameters` in place so a YAML-built config runs under
the (quadrature, float, architecture) tuple `spec`. Returns `params`.
"""
function apply_overrides!(params, spec)
    params.quadrature_type =
        spec.quad === :RadauQuad           ? CoreRT.RadauQuad() :
        spec.quad === :GaussLegQuad ? CoreRT.GaussLegQuad() :
        error("unknown quadrature $(spec.quad)")
    params.float_type   = spec.float
    params.architecture = spec.arch_ctor()
    return params
end

"""
    rel_error(modeled, truth; near_zero_atol=1e-9)

Relative error matrix; positions where `|truth| < near_zero_atol` are masked
(returned as `NaN`) so they don't dominate the maximum.
"""
function rel_error(modeled::AbstractArray, truth::AbstractArray;
                   near_zero_atol::Real = 1e-9)
    size(modeled) == size(truth) || error("size mismatch: $(size(modeled)) vs $(size(truth))")
    out = similar(modeled, Float64)
    for i in eachindex(modeled, truth)
        if abs(truth[i]) < near_zero_atol
            out[i] = NaN
        else
            out[i] = abs(modeled[i] - truth[i]) / abs(truth[i])
        end
    end
    return out
end

"""
    summarize(label, modeled, truth; rtol, name="value")

Print summary statistics and return `(passes::Bool, max_rel)` for use with
`@test`. Skips NaN-masked entries.
"""
function summarize(label::AbstractString, modeled, truth; rtol::Real,
                   name::AbstractString = "value")
    re = rel_error(modeled, truth)
    valid = filter(!isnan, vec(re))
    if isempty(valid)
        @warn "$label: no comparable (non-zero-truth) entries"
        return (true, 0.0)
    end
    mx, mn = maximum(valid), median(valid)
    pass = mx < rtol
    flag = pass ? "PASS" : "FAIL"
    @info "$label [$flag]" name n=length(valid) median_rel=round(mn, sigdigits=3) max_rel=round(mx, sigdigits=3) rtol
    return (pass, mx)
end
