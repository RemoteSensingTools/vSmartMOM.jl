using Printf
using Statistics
using Dates

using CUDA
using ForwardDiff
using vSmartMOM
using vSmartMOM.StandaloneSS

const RESULT_PATH = joinpath(@__DIR__, "standalone_ss_o2a_gpu_scalar_benchmark_results.txt")

function _check_cuda!()
    CUDA.functional() || error("CUDA is not functional in this Julia session")
    vSmartMOM.Architectures.has_cuda() ||
        error("vSmartMOM CUDA extension is not available; load CUDA before GPU()")
    return nothing
end

function _layer_weights(::Type{FT}, n_layers::Int; bottom_heavy::Bool) where {FT}
    x = collect(range(FT(0), FT(1); length=n_layers))
    raw = bottom_heavy ? exp.(FT(2.0) .* x) : exp.(-FT(1.5) .* x)
    return raw ./ sum(raw)
end

function _o2a_like_optics(::Type{FT}, n_layers::Int, n_spec::Int) where {FT}
    nu = collect(range(FT(12950), FT(13200); length=n_spec))
    nu0 = first(nu)
    mid = (nu0 + last(nu)) / FT(2)
    span = last(nu) - nu0

    line = zeros(FT, n_spec)
    centers = collect(range(first(nu) + FT(6), last(nu) - FT(6); length=56))
    @inbounds for (ic, c) in enumerate(centers)
        strength = FT(0.35) + FT(0.65) * abs(sin(FT(0.73) * FT(ic)))
        width = FT(0.035) + FT(0.065) * abs(cos(FT(0.41) * FT(ic)))
        inv_width2 = inv(width * width)
        for ispec in eachindex(nu)
            x = (nu[ispec] - c)
            line[ispec] += strength * exp(-FT(0.5) * x * x * inv_width2)
        end
    end
    line ./= maximum(line)

    broad = @. FT(0.18) + FT(0.82) * exp(-FT(0.5) * ((nu - mid) / (FT(0.18) * span))^2)
    ripple = @. FT(1) + FT(0.035) * sin(FT(2π) * (nu - nu0) / span * FT(7))

    abs_w = _layer_weights(FT, n_layers; bottom_heavy=true)
    scat_w = _layer_weights(FT, n_layers; bottom_heavy=false)

    tau_abs = zeros(FT, n_layers, n_spec)
    tau_ray = zeros(FT, n_layers, n_spec)
    tau_aer = zeros(FT, n_layers, n_spec)

    @inbounds for ispec in 1:n_spec
        ray_shape = (nu[ispec] / mid)^FT(4)
        aer_shape = (nu[ispec] / mid)^FT(1.1)
        abs_shape = FT(8.0) * line[ispec] + FT(1.3) * broad[ispec]
        for iz in 1:n_layers
            tau_abs[iz, ispec] = FT(0.0015) + abs_w[iz] * abs_shape
            tau_ray[iz, ispec] = scat_w[iz] * FT(0.055) * ray_shape * ripple[ispec]
            tau_aer[iz, ispec] = scat_w[iz] * FT(0.018) * aer_shape
        end
    end

    return (; nu, tau_abs, tau_ray, tau_aer)
end

function _make_config(::Type{FT};
                      n_spec::Int,
                      n_layers::Int,
                      n_geom::Int,
                      inner_nquad::Int,
                      azimuth_nquad::Int,
                      architecture=vSmartMOM.GPU(),
                      abs_scale=one(FT)) where {FT}
    optics = _o2a_like_optics(FT, n_layers, n_spec)
    mu_v = collect(range(FT(0.18), FT(0.92); length=n_geom))
    dphi = collect(range(FT(0.0), FT(pi); length=n_geom))
    geometry = SSGeometry(μ₀=FT(0.78), μv=mu_v, Δϕ=dphi)
    nu0 = first(optics.nu)
    nu_span = last(optics.nu) - nu0
    albedo = @. FT(0.21) + FT(0.035) * sin(FT(2π) *
        (optics.nu - nu0) / nu_span)
    surface = LambertianSSSurface(albedo=albedo)
    rayleigh = RayleighSSContributor(τ=optics.tau_ray)
    aerosol = HGAerosolSSContributor(g=FT(0.72), ϖ=FT(0.91), τ=optics.tau_aer)
    absorption = AbsorptionSSContributor(τ=optics.tau_abs .* abs_scale)
    return ExactSSConfig(
        geometry=geometry,
        surface=surface,
        contributors=(rayleigh, aerosol, absorption),
        I0=fill(one(FT), n_spec),
        n_stokes=1,
        inner_nquad=inner_nquad,
        azimuth_nquad=azimuth_nquad,
        architecture=architecture)
end

function _sync_elapsed_seconds(f)
    GC.gc(false)
    CUDA.reclaim()
    CUDA.synchronize()
    t0 = time_ns()
    out = f()
    CUDA.synchronize()
    return (time_ns() - t0) / 1e9, out
end

function _result_summary(result)
    total_h = Array(result.total)
    path1_h = Array(result.path1)
    path2_h = Array(result.path2)
    path3_h = Array(result.path3)
    path4_h = Array(result.path4)
    return (;
        total_sum = Float64(sum(total_h)),
        total_min = Float64(minimum(total_h)),
        total_max = Float64(maximum(total_h)),
        path1_sum = Float64(sum(path1_h)),
        path2_sum = Float64(sum(path2_h)),
        path3_sum = Float64(sum(path3_h)),
        path4_sum = Float64(sum(path4_h)),
        total_size = size(total_h),
        finite = all(isfinite, total_h))
end

function _benchmark_case(config; paths::Symbol, reps::Int)
    warm = run_exact_ss(config; paths=paths)
    CUDA.synchronize()
    warm = nothing
    GC.gc(false)
    CUDA.reclaim()

    times = Float64[]
    last_result = nothing
    for _ in 1:reps
        elapsed, result = _sync_elapsed_seconds() do
            run_exact_ss(config; paths=paths)
        end
        push!(times, elapsed)
        last_result = result
    end
    summary = _result_summary(last_result)
    last_result = nothing
    GC.gc(false)
    CUDA.reclaim()
    return (; paths, times, min=minimum(times), median=median(times),
            summary)
end

function _gpu_analytic_jacobian_probe(config)
    try
        run_exact_ss_with_jacobians(config; paths=:paths_1_2)
        CUDA.synchronize()
        return "unexpectedly succeeded"
    catch err
        return "$(typeof(err)): $(sprint(showerror, err))"
    end
end

function _gpu_forwarddiff_probe(; n_spec::Int, n_layers::Int, n_geom::Int,
                                inner_nquad::Int, azimuth_nquad::Int)
    warmup(x) = begin
        cfg = _make_config(Float32; n_spec=128, n_layers=3, n_geom=2,
                           inner_nquad=4, azimuth_nquad=4,
                           architecture=vSmartMOM.GPU(),
                           abs_scale=x[1])
        result = run_exact_ss(cfg; paths=:paths_1_2)
        CUDA.synchronize()
        vec(Array(result.total[:, :, 1:8]))
    end

    try
        ForwardDiff.jacobian(warmup, Float32[1])
        f(x) = begin
            cfg = _make_config(Float32; n_spec=n_spec, n_layers=n_layers,
                               n_geom=n_geom, inner_nquad=inner_nquad,
                               azimuth_nquad=azimuth_nquad,
                               architecture=vSmartMOM.GPU(),
                               abs_scale=x[1])
            result = run_exact_ss(cfg; paths=:paths_1_2)
            CUDA.synchronize()
            vec(Array(result.total))
        end
        elapsed = @elapsed jac = ForwardDiff.jacobian(f, Float32[1])
        return (status="ok", elapsed=elapsed, size=size(jac),
                finite=all(isfinite, jac),
                sample=Float64(jac[1, 1]))
    catch err
        return (status="failed", error="$(typeof(err)): $(sprint(showerror, err))")
    end
end

function _print_case(io, label, case)
    @printf(io, "%-14s paths=%-10s min=%.6f s median=%.6f s reps=%s\n",
            label, string(case.paths), case.min, case.median, repr(case.times))
    s = case.summary
    @printf(io, "  size=%s finite=%s total[min,max,sum]=[%.6e, %.6e, %.6e]\n",
            string(s.total_size), string(s.finite), s.total_min, s.total_max,
            s.total_sum)
    @printf(io, "  component sums: p1=%.6e p2=%.6e p3=%.6e p4=%.6e\n",
            s.path1_sum, s.path2_sum, s.path3_sum, s.path4_sum)
end

function main()
    _check_cuda!()
    FT = Float32
    n_spec = parse(Int, get(ENV, "SS_O2A_NSPEC", "20000"))
    n_layers = parse(Int, get(ENV, "SS_O2A_NLAYERS", "8"))
    n_geom = parse(Int, get(ENV, "SS_O2A_NGEOM", "16"))
    inner_nquad = parse(Int, get(ENV, "SS_O2A_INNER_NQUAD", "8"))
    azimuth_nquad = parse(Int, get(ENV, "SS_O2A_AZIMUTH_NQUAD", "12"))
    reps = parse(Int, get(ENV, "SS_O2A_REPS", "2"))

    config = _make_config(FT; n_spec, n_layers, n_geom, inner_nquad,
                          azimuth_nquad, architecture=vSmartMOM.GPU())

    cases = Any[]
    push!(cases, ("paths 1/2", _benchmark_case(config; paths=:paths_1_2,
                                                reps=reps)))
    push!(cases, ("all 1/2/3/4", _benchmark_case(config; paths=:all,
                                                  reps=reps)))

    jac_gpu = _gpu_analytic_jacobian_probe(config)
    ad_gpu = _gpu_forwarddiff_probe(; n_spec, n_layers, n_geom, inner_nquad,
                                    azimuth_nquad)

    open(RESULT_PATH, "w") do io
        println(io, "StandaloneSS scalar exact-pass GPU benchmark")
        println(io, "date=$(Dates.now())")
        println(io, "gpu=$(CUDA.name(CUDA.device()))")
        println(io, "julia=$(VERSION)")
        println(io, "n_spec=$n_spec n_layers=$n_layers n_geom=$n_geom")
        println(io, "inner_nquad=$inner_nquad azimuth_nquad=$azimuth_nquad reps=$reps")
        println(io, "note=individual :path3/:path4 timings are intentionally not run at 20k here; the combined :all path uses the GPU packed path3+path4 precompute.")
        println(io)
        for (label, case) in cases
            _print_case(io, label, case)
        end
        println(io)
        println(io, "GPU analytic Jacobian probe: $jac_gpu")
        println(io, "GPU ForwardDiff probe: $ad_gpu")
    end

    print(read(RESULT_PATH, String))
    return RESULT_PATH
end

main()
