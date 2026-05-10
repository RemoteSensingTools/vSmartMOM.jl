using Printf
using Statistics
using vSmartMOM
using vSmartMOM.StandaloneSS

const CUDA_LOADED = Ref(false)
try
    @eval using CUDA
    CUDA_LOADED[] = true
catch
end

function _cuda_ok()
    return CUDA_LOADED[] && CUDA.functional() &&
           vSmartMOM.Architectures.has_cuda()
end

function _spectral_config(::Type{FT}; n_geom::Int, n_spec::Int,
                          n_layers::Int, architecture=vSmartMOM.CPU()) where {FT}
    geometry = SSGeometry(
        μ₀ = FT(0.82),
        μv = collect(range(FT(0.25), FT(0.9); length=n_geom)),
        Δϕ = collect(range(FT(0.05), FT(2.9); length=n_geom)))
    surface = LambertianSSSurface(albedo=fill(FT(0.22), n_spec))
    rayleigh = RayleighSSContributor(τ=fill(FT(0.03), n_layers, n_spec))
    aerosol = HGAerosolSSContributor(g=FT(0.72), ϖ=FT(0.9),
                                     τ=fill(FT(0.01), n_layers, n_spec))
    absorption = AbsorptionSSContributor(τ=fill(FT(0.02), n_layers, n_spec))
    return ExactSSConfig(geometry=geometry, surface=surface,
                         contributors=(rayleigh, aerosol, absorption),
                         I0=fill(one(FT), n_spec), inner_nquad=8,
                         azimuth_nquad=16, architecture=architecture)
end

function _avg_elapsed(f, n::Int)
    samples = Float64[]
    for _ in 1:n
        GC.gc(false)
        push!(samples, @elapsed f())
    end
    return mean(samples), minimum(samples)
end

function _run_case(; n_geom::Int, n_spec::Int, n_layers::Int, nrep::Int)
    cpu_config = _spectral_config(Float64; n_geom, n_spec, n_layers)
    run_exact_ss(cpu_config; paths=:all)

    cpu_avg, cpu_min = _avg_elapsed(nrep) do
        run_exact_ss(cpu_config; paths=:all)
        nothing
    end

    if _cuda_ok()
        gpu_config = _spectral_config(Float64; n_geom, n_spec, n_layers,
                                      architecture=vSmartMOM.GPU())
        run_exact_ss(gpu_config; paths=:all)
        CUDA.synchronize()
        gpu_avg, gpu_min = _avg_elapsed(nrep) do
            run_exact_ss(gpu_config; paths=:all)
            CUDA.synchronize()
            nothing
        end
        cpu_result = run_exact_ss(cpu_config; paths=:all)
        gpu_result = run_exact_ss(gpu_config; paths=:all)
        CUDA.synchronize()
        err = maximum(abs.(Array(gpu_result.total) .- cpu_result.total))
        @printf("%8d %8d %8d  CPU avg/min %.4f %.4f  GPU avg/min %.4f %.4f  speedup %.2fx  maxerr %.3e\n",
                n_geom, n_spec, n_layers, cpu_avg, cpu_min, gpu_avg, gpu_min,
                cpu_min / gpu_min, err)
    else
        @printf("%8d %8d %8d  CPU avg/min %.4f %.4f  GPU unavailable\n",
                n_geom, n_spec, n_layers, cpu_avg, cpu_min)
    end
end

function main()
    println("StandaloneSS spectral scaling benchmark")
    println("Columns: nGeom nSpec nLayer timings")
    cases = ((16, 256, 6),
             (16, 1024, 6),
             (16, 4096, 6),
             (32, 1024, 8))
    for (n_geom, n_spec, n_layers) in cases
        _run_case(; n_geom, n_spec, n_layers, nrep=3)
    end
end

main()
