# Warmed A/B benchmark of reference and batched Jacobian propagation.
using vSmartMOM, vSmartMOM.CoreRT, YAML, LinearAlgebra, Statistics, Logging, Profile
backend = get(ENV,"AUDIT_BACKEND","cpu")
n_spec = parse(Int,get(ENV,"AUDIT_NSPEC","64"))
label = "$(backend)-n$(n_spec)-t$(Threads.nthreads())"
if backend == "cuda"
    @eval using CUDA
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
function main()
    d=YAML.load_file(joinpath(pkgdir(vSmartMOM),"test","test_parameters","JacobianTestFast.yaml"))
    delete!(d,"absorption")
    d["radiative_transfer"]["architecture"] = backend == "cuda" ? "GPU()" : "CPU()"
    d["radiative_transfer"]["numerics"] = Dict("verbose"=>true,"blas_threads"=>1)
    d["radiative_transfer"]["greek_beta_cutoff"] = nothing
    d["radiative_transfer"]["polarization_type"] = get(ENV,"AUDIT_POL","Stokes_I()")
    d["scattering"]["r_max"]=3.0
    d["scattering"]["aerosols"][1]["μ"]=0.15
    d["scattering"]["aerosols"][1]["σ"]=1.4
    p=read_parameters(d)
    p.spec_bands[1]=collect(range(first(p.spec_bands[1]),last(p.spec_bands[1]);length=n_spec))
    println("BENCH_BEGIN ",label); flush(stdout)
    m,l=quiet(()->model_from_parameters(LinMode(),p))
    ng=size(l.τ̇_abs[1],1)
    native=ParameterLayout(n_aerosols=1,n_gases=ng,n_surface=1)
    keys=[ParameterKey(:atmosphere,:surface_pressure),
          ParameterKey(:aerosol,:tau_ref;component=1),
          ParameterKey(:aerosol,:profile_location;component=1),
          ParameterKey(:surface,:albedo;component=1,band=1)]
    names=["surface_pressure","aerosol_tau_ref","aerosol_profile_location","surface_albedo"]
    layout=ActiveParameterLayout(keys,names,keys,names,[1,2,7],3;surface_columns=4:4)
    plan=JacobianPlan(OCO_RRS_synth(),keys,names,[layout])
    selected=PlannedRTModelLin(l,plan)
    fwd=()->rt_run(m)
    full=()->rt_run(m,l,1,ng,1)
    sel=()->rt_run(m,selected)
    old = CoreRT._BATCHED_JACOBIANS_ENABLED[]
    results = []
    outputs = []
    try
        if get(ENV,"AUDIT_TIMER_ONLY","false") == "true"
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = true
            quiet(full); sync_backend()
            if get(ENV,"AUDIT_HOST_PROFILE","false") == "true"
                Profile.clear()
                @profile quiet(full)
                sync_backend()
                Profile.print(stdout; format=:flat, C=true, sortedby=:count, mincount=10)
                Profile.print(stdout; format=:tree, C=false, maxdepth=30, mincount=50)
            else
                with_logger(NullLogger()) do
                    full(); sync_backend()
                end
            end
            println("TIMER_DONE ",label)
            return
        end
        for enabled in (false,true)
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = enabled
            if !enabled && get(ENV,"AUDIT_FAST_REFERENCE","false") == "true"
                println("BENCH_PHASE reference parity only"); flush(stdout)
                push!(outputs, quiet(full))
                continue
            end
            GC.gc(); backend == "cuda" && CUDA.reclaim()
            println("BENCH_PHASE enabled=",enabled," kind=full"); flush(stdout)
            full_times=measure(full)
            GC.gc(); backend == "cuda" && CUDA.reclaim()
            println("BENCH_PHASE enabled=",enabled," kind=selected"); flush(stdout)
            selected_times=measure(sel)
            push!(outputs,quiet(full))
            push!(results,(enabled,full_times,selected_times))
        end
        CoreRT._BATCHED_JACOBIANS_ENABLED[]=true
        GC.gc(); backend == "cuda" && CUDA.reclaim()
        println("BENCH_PHASE kind=forward"); flush(stdout)
        forward_times=measure(fwd)
        rf=quiet(fwd)
        rs=quiet(sel)
        columns=[1,2,7,CoreRT.n_total(native)]
        for i in 1:4
            @assert isapprox(outputs[1][i],outputs[2][i];rtol=1e-9,atol=1e-11)
        end
        @assert isapprox(rs[3],outputs[2][3][:,:,:,columns];rtol=1e-9,atol=1e-11)
        println("SELECTED_PARITY dR=",maximum(abs,rs[3]-outputs[2][3][:,:,:,columns]))
        println("MODEL backend=",backend," nSpec=",n_spec," operator=",size(m.quad_points.qp_μ,1)*p.polarization_type.n,
            " full_columns=",size(outputs[1][3],4)," selected_columns=4")
        println("PARITY R=",maximum(abs,outputs[1][1]-outputs[2][1]),
            " T=",maximum(abs,outputs[1][2]-outputs[2][2]),
            " dR=",maximum(abs,outputs[1][3]-outputs[2][3]),
            " dT=",maximum(abs,outputs[1][4]-outputs[2][4]),
            " forward_R=",maximum(abs,rf[1]-outputs[2][1]))
        for (enabled,ft,st) in results
            for (kind,ts) in [("full",ft),("selected",st)]
                println("MEASURE enabled=",enabled," kind=",kind,
                    " median_seconds=",median(t.time for t in ts),
                    " median_bytes=",median(t.bytes for t in ts),
                    " median_gpu_bytes=",median(t.gpu_bytes for t in ts)," samples=",ts)
            end
        end
        println("FORWARD median_seconds=",median(t.time for t in forward_times))
        if get(ENV,"AUDIT_TIMERS","false")=="true"
            with_logger(NullLogger()) do
                full();sync_backend()
            end
        end
        if get(ENV,"AUDIT_CUPTI","false")=="true"
            result=device_profile(()->quiet(full));sync_backend()
            show(stdout,MIME("text/plain"),result);println()
        end
    finally
        CoreRT._BATCHED_JACOBIANS_ENABLED[]=old
    end
    println("BENCH_DONE ",label)
end
main()
