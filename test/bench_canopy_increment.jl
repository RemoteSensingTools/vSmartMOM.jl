# Run from test/: julia --project=.. bench_canopy_increment.jl
# Set VSMART_CANOPY_BENCH_GPU=true to select CUDA instead of CPU.
# Timings exclude phase-cache construction, elemental seeding and host transfers.
# The ordinary method is a timing baseline; its thin Float32 seeds lose accuracy.
using TOML, LinearAlgebra, Statistics, vSmartMOM, CanopyOptics
using vSmartMOM.CoreRT, vSmartMOM.InelasticScattering

const BENCH = TOML.parsefile(joinpath(@__DIR__,"canopy_fixtures/increment_benchmark.toml"))

function bench_context(::Type{FT},streams,npol,nspec,architecture) where {FT}
    g,o = BENCH["geometry"],BENCH["optics"]
    pol = npol==1 ? CoreRT.Stokes_I{FT}() : CoreRT.Stokes_IQUV{FT}()
    obs = CoreRT.ObsGeometry(FT(g["solar_zenith_degrees"]),FT[0],FT[0],FT[0])
    arr = vSmartMOM.Architectures.array_type(architecture)
    qp = CoreRT.rt_set_streams(CoreRT.GaussLegQuad(),2streams-1,obs,pol,arr)
    dim = length(qp.qp_μ)*npol
    added = CoreRT.make_added_layer(noRS{FT}(),FT,arr,(dim,dim),nspec)
    canopy = CoreRT.CanopySurface(;soil=CoreRT.LambertianSurfaceScalar(FT(o["soil_albedo"])),
        LAI=FT(g["lai"]),LAD=beta_leaves(FT(g["leaf_beta_alpha"]),FT(g["leaf_beta_beta"])),
        canopy_scattering=BiLambertianCanopyScattering(R=FT(o["leaf_reflectance"]),T=FT(o["leaf_transmittance"])),
        leaf_reflectance=FT(o["leaf_reflectance"]),leaf_transmittance=FT(o["leaf_transmittance"]),
        canopy_quadrature=CanopyQuadrature(n_leaf=g["leaf_nodes"],n_azimuth=g["azimuth_nodes"]))
    m = g["fourier_order"]
    CoreRT._init_canopy_cache!(canopy,added,pol,qp,architecture;m_max=m)
    cache = canopy._cache
    core = CoreRT.CoreDirectionalScatteringOpticalProperties(arr(fill(FT(g["lai"]),nspec)),
        arr(fill(FT(o["single_scattering_albedo"]),nspec)),
        arr(cache.Z_coarse_pp[:,:,1,m+1]),arr(cache.Z_coarse_mp[:,:,1,m+1]),cache.G)
    F0 = arr(repeat(FT.(o["solar_stokes"][1:npol]),1,nspec))
    (;pol,qp,dim,nspec,arr,architecture,cache,core,F0,m)
end

function bench_seed(ctx,::Type{FT}) where {FT}
    n = BENCH["numerics"]
    dt,nd,beam = CoreRT.init_layer(ctx.core,ctx.qp,ctx.pol,ctx.architecture;
        dτ_min_floor=FT(n["thin_floor"]),dτ_max_threshold=FT(n["thin_threshold"]))
    added = CoreRT.make_added_layer(noRS{FT}(),FT,ctx.arr,(ctx.dim,ctx.dim),ctx.nspec)
    CoreRT.elemental!(ctx.pol,true,ctx.arr(zeros(FT,ctx.nspec)),dt,ctx.F0,
        ctx.core,ctx.m,nd,true,ctx.qp,added,ctx.architecture)
    (;added,dt,nd,beam)
end

function bench_double!(ctx,s,increment)
    if increment
        CoreRT.canopy_increment_doubling!(ctx.pol,true,s.dt,ctx.core,ctx.m,s.nd,
            ctx.qp,s.added,ctx.cache.I_static,ctx.architecture)
    else
        CoreRT.doubling!(ctx.pol,true,s.beam,s.nd,s.added,ctx.cache.I_static,ctx.architecture)
    end
    nothing
end

function measure_doubling(ctx,::Type{FT},increment,synchronize) where {FT}
    n = BENCH["numerics"]
    times,bytes = Float64[],Int[]
    nd = 0
    for repetition in 1:n["warmup_calls"]+n["samples"]
        GC.gc()
        s = bench_seed(ctx,FT)
        nd = s.nd
        synchronize()
        measurement = @timed begin
            bench_double!(ctx,s,increment)
            synchronize()
        end
        all(isfinite,Array(s.added.t⁺⁺)) || error("Nonfinite benchmark transmission")
        all(isfinite,Array(s.added.r⁻⁺)) || error("Nonfinite benchmark reflection")
        if repetition>n["warmup_calls"]
            push!(times,measurement.time)
            push!(bytes,measurement.bytes)
        end
    end
    (;nd,minimum_seconds=minimum(times),median_seconds=median(times),minimum_host_bytes=minimum(bytes))
end

function main(architecture,synchronize)
    old_threads = BLAS.get_num_threads()
    try
        BLAS.set_num_threads(BENCH["numerics"]["blas_threads"])
        println("backend,precision,streams,stokes,spectral_points,operator_dimension,doublings,method,min_seconds,median_seconds,min_host_bytes")
        g = BENCH["geometry"]
        for FT in (Float32,Float64),streams in g["streams"],npol in g["stokes_components"],nspec in g["spectral_counts"]
            ctx = bench_context(FT,streams,npol,nspec,architecture)
            for increment in (false,true)
                r = measure_doubling(ctx,FT,increment,synchronize)
                println(join((nameof(typeof(architecture)),FT,streams,npol,nspec,ctx.dim,r.nd,
                    increment ? "increment" : "ordinary",r.minimum_seconds,r.median_seconds,r.minimum_host_bytes),','))
                flush(stdout)
            end
        end
    finally
        BLAS.set_num_threads(old_threads)
    end
end

if get(ENV,"VSMART_CANOPY_BENCH_GPU","false")=="true"
    using CUDA
    CUDA.functional() || error("CUDA requested but unavailable")
    CUDA.allowscalar(false)
    @info "Canopy doubling benchmark" device=CUDA.name(CUDA.device()) julia=VERSION julia_threads=Threads.nthreads()
    main(vSmartMOM.Architectures.GPU(),()->CUDA.synchronize())
else
    @info "Canopy doubling benchmark" cpu=Sys.CPU_NAME julia=VERSION julia_threads=Threads.nthreads()
    main(vSmartMOM.Architectures.CPU(),()->nothing)
end
