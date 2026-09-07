# Warmed complete-solve comparison of physical and local optical directions.
using TOML, SHA
using vSmartMOM, vSmartMOM.CoreRT, YAML, LinearAlgebra, Statistics, Logging, Profile
backend = get(ENV,"AUDIT_BACKEND","cpu")
n_spec = parse(Int,get(ENV,"AUDIT_NSPEC","64"))
const chunks = parse(Int,get(ENV,"AUDIT_CHUNKS","1"))
const external_solar = get(ENV,"AUDIT_EXTERNAL_SOLAR","false")=="true"
label = "$(backend)-n$(n_spec)-t$(Threads.nthreads())"
if backend == "cuda"
    @eval using CUDA
    CUDA.allowscalar(false)
    CUDA.functional() || error("CUDA required for this benchmark")
    @eval device_profile(f) = CUDA.@profile external=false f()
    @eval device_timed(f) = CUDA.@timed f()
end
sync_backend() = backend == "cuda" ? CUDA.synchronize() : nothing
BLAS.set_num_threads(1)
function quiet(f)
    redirect_stdout(devnull) do
        with_logger(NullLogger()) do
            f()
        end
    end
end
function measure(f)
    quiet(f); sync_backend() # warm this exact path
    results = map(1:3) do _
        GC.gc(); sync_backend()
        t = if backend == "cuda"
            device_timed(()->quiet(f))
        else
            @timed quiet(f)
        end
        sync_backend()
        (;time=t.time, bytes=backend == "cuda" ? t.cpu_bytes : t.bytes,
          gctime=backend == "cuda" ? t.cpu_gctime : t.gctime,
          gpu_bytes=backend == "cuda" ? t.gpu_bytes : 0)
    end
    return results
end

# Timing starts at rt_run: Mie, spectroscopy, and upstream optical derivatives
# are prepared before warming. Optical mixing and the complete surface-coupled
# atmospheric solve are included. Chunking is identical for both tangent paths.
function fixture(chunk)
    cfg=YAML.load_file(joinpath(pkgdir(vSmartMOM),"test/test_parameters/JacobianTestFast.yaml"))
    gases=filter(!isempty,split(get(ENV,"AUDIT_GASES",""),","))
    if isempty(gases)
        delete!(cfg,"absorption")
    else
        dry_gases=filter(!=("H2O"),gases)
        cfg["absorption"]["molecules"]=[dry_gases]
        cfg["absorption"]["variable_molecules"]=[dry_gases]
        if "H2O" in gases
            # Specific humidity, kg/kg: moist lower troposphere, dry aloft.
            pressure=cfg["atmospheric_profile"]["p"]
            cfg["atmospheric_profile"]["q"]=[0.007*((pressure[z]+pressure[z+1])/2000)^3
                for z in eachindex(cfg["atmospheric_profile"]["T"])]
        end
        cfg["absorption"]["vmr"]=Dict("O2"=>0.21,"CO2"=>420e-6,"H2O"=>0.005,"CH4"=>1.8e-6)
    end
    cfg["radiative_transfer"]["architecture"]=backend == "cuda" ? "GPU()" : "CPU()"
    cfg["radiative_transfer"]["polarization_type"]=get(ENV,"AUDIT_POL","Stokes_IQU()")
    cfg["radiative_transfer"]["nstreams"]=parse(Int,get(ENV,"AUDIT_STREAMS","3"))
    cfg["radiative_transfer"]["greek_beta_cutoff"]=nothing
    cfg["radiative_transfer"]["numerics"]=Dict("verbose"=>get(ENV,"AUDIT_TIMERS","false")=="true","blas_threads"=>1)
    cfg["atmospheric_profile"]["profile_reduction"]=parse(Int,get(ENV,"AUDIT_LAYERS","5"))
    cfg["geometry"]["vaz"]=[0.0,37.0]
    cfg["scattering"]["r_max"]=3.0
    cfg["scattering"]["aerosols"][1]["μ"]=0.15
    cfg["scattering"]["aerosols"][1]["σ"]=1.4
    p=read_parameters(cfg)
    lo=parse(Float64,get(ENV,"AUDIT_NU_MIN",string(first(p.spec_bands[1]))))
    hi=parse(Float64,get(ENV,"AUDIT_NU_MAX",string(last(p.spec_bands[1]))))
    grid=range(lo,hi;length=n_spec*chunks)
    p.spec_bands[1]=collect(grid[(chunk-1)*n_spec+1:chunk*n_spec])
    model,lin=quiet(()->model_from_parameters(LinMode(),p;external_solar))
    @assert model.quad_points.external_solar == external_solar
    return model,lin
end
function main()
    compare_medium=get(ENV,"AUDIT_COMPARE_MEDIUM","false")=="true"
    compare_vendor=get(ENV,"AUDIT_COMPARE_VENDOR","false")=="true"
    compare_tiles=get(ENV,"AUDIT_COMPARE_TILES","false")=="true"
    @assert count(identity,(compare_medium,compare_vendor,compare_tiles)) <= 1
    modes=compare_vendor ? (:blocked,:physical,:local,:forward) :
        (compare_medium || compare_tiles) ?
            (:reference,:physical,:local,:forward) : (:physical,:local,:forward)
    get(ENV,"AUDIT_PROFILE_ONLY","false")=="true" && (modes=(:local,:forward))
    totals=Dict(string(k)=>zeros(3) for k in modes if k != :blocked)
    records=Any[]
    for chunk in 1:chunks
        println("PREPARE chunk=$chunk/$chunks"); flush(stdout)
        model,lin=fixture(chunk)
        ng=size(lin.τ̇_abs[1],1)
        println("ABSORPTION minmax=$(extrema(model.τ_abs[1])) nonzero_gas_columns=$(count(p->any(!iszero,view(lin.τ̇_abs[1],p,:,:)),1:ng))/$ng")
        flush(stdout)
        baseline=nothing
        for mode in modes
            if compare_vendor
                CoreRT._MEDIUM_JACOBIANS_ENABLED[] = true
                CoreRT._VENDOR_JACOBIANS_ENABLED[] = mode != :blocked
            elseif compare_medium
                CoreRT._VENDOR_JACOBIANS_ENABLED[] = false
                CoreRT._MEDIUM_JACOBIANS_ENABLED[] = mode != :reference
            elseif compare_tiles
                CoreRT._TILED_JACOBIANS_ENABLED[] = mode != :reference
            end
            f=mode == :forward ? ()->(external_solar ? (rt_run_toa(model),nothing) : rt_run(model)) :
                ()->rt_run(model,lin,1,ng,1;jacobian_basis=mode in (:reference,:blocked) ? :physical : mode)
            println("MEASURE_BEGIN chunk=$chunk mode=$mode");flush(stdout)
            GC.gc(); backend == "cuda" && CUDA.reclaim()
            # The blocked parity reference is evaluated once. Its warmed
            # timings were measured separately; do not report a fake zero time.
            ts=mode == :blocked ? [] : measure(f)
            mode == :blocked || (totals[string(mode)] .+= [t.time for t in ts])
            out=quiet(f); sync_backend()
            err=0.0
            if baseline === nothing
                baseline=out
            else
                for i in 1:(mode == :forward ? 2 : 4)
                    a,b=out[i],baseline[i]
                    if a === nothing || b === nothing
                        @assert a === b
                        continue
                    end
                    @assert isapprox(a,b;rtol=1e-9,atol=1e-11)
                    err=max(err,maximum(abs,a-b))
                end
            end
            record=Dict("chunk"=>chunk,"mode"=>string(mode),"seconds"=>[t.time for t in ts],
                "host_bytes"=>[t.bytes for t in ts],"device_allocated_bytes"=>[t.gpu_bytes for t in ts],
                "max_parity_error"=>err)
            push!(records,record)
            println("MEASURE ",record);flush(stdout)
            if mode == :local && get(ENV,"AUDIT_PROFILE_ONLY","false")=="true"
                if backend == "cuda"
                    println("CUDA_PROFILE_BEGIN")
                    prof=device_profile(()->quiet(f));sync_backend()
                    show(stdout,MIME("text/plain"),prof);println()
                end
                Profile.clear()
                @profile begin quiet(f);sync_backend() end
                println("HOST_PROFILE_BEGIN")
                Profile.print(stdout;format=:flat,C=true,sortedby=:count,mincount=10)
            end
            if get(ENV,"AUDIT_TIMERS","false") == "true"
                println("STAGES chunk=$chunk mode=$mode")
                with_logger(f,NullLogger()); sync_backend()
            end
        end
        println("SCENE nSpec=$(n_spec*chunks) chunk_size=$n_spec layers=$(length(model.profile.p_full)) streams=$(model.quad_points.Nstreams) stokes=$(CoreRT.polarization_type(model).n) external_solar=$external_solar columns=$(size(baseline[3],4))")
    end
    source_files = sort([joinpath(root,file) for folder in ("src","ext")
        for (root,_,files) in walkdir(joinpath(pkgdir(vSmartMOM),folder))
        for file in files if endswith(file,".jl")])
    fingerprint = join([relpath(file,pkgdir(vSmartMOM))*":"*bytes2hex(sha256(read(file)))
        for file in source_files],"\n")
    result=Dict("julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),
        "source_sha256"=>bytes2hex(sha256(fingerprint)),
        "configuration"=>Dict(k=>v for (k,v) in ENV if startswith(k,"AUDIT_")),
        "backend"=>backend,"spectral_points"=>n_spec*chunks,"chunk_size"=>n_spec,
        "measurements"=>records,"totals_seconds"=>totals,
        "medians_seconds"=>Dict(k=>median(v) for (k,v) in totals))
    println("TOTAL ",result["medians_seconds"])
    file=get(ENV,"AUDIT_OUTPUT","/tmp/local-basis-benchmark.toml")
    open(io->TOML.print(io,result),file,"w")
end
main()
