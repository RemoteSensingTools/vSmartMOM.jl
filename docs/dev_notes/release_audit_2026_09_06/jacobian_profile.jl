# Warmed forward/full/selected Jacobian benchmark; no solver modifications.
using vSmartMOM, vSmartMOM.CoreRT, YAML, LinearAlgebra, Statistics, Logging, Profile
backend = get(ENV,"AUDIT_BACKEND","cpu")
n_spec = parse(Int,get(ENV,"AUDIT_NSPEC","64"))
label = "$(backend)-n$(n_spec)-t$(Threads.nthreads())"
outdir = get(ENV,"AUDIT_OUTPUT","/tmp/vsmartmom-release-evidence")
if backend == "cuda"
    @eval using CUDA
    CUDA.functional() || error("CUDA required for this benchmark")
    @eval device_profile(f) = CUDA.@profile external=false f()
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
        t = @timed begin
            quiet(f)
            sync_backend()
        end
        (;time=t.time, bytes=t.bytes, gctime=t.gctime)
    end
    return results
end
function main()
    d=YAML.load_file(joinpath(pkgdir(vSmartMOM),"test","test_parameters","JacobianTestFast.yaml"))
    delete!(d,"absorption")
    d["radiative_transfer"]["architecture"] = backend == "cuda" ? "GPU()" : "CPU()"
    d["radiative_transfer"]["numerics"] = Dict("verbose"=>true,"blas_threads"=>1)
    d["radiative_transfer"]["greek_beta_cutoff"] = nothing
    d["scattering"]["r_max"]=3.0
    d["scattering"]["aerosols"][1]["μ"]=0.15
    d["scattering"]["aerosols"][1]["σ"]=1.4
    p=read_parameters(d)
    p.spec_bands[1]=collect(range(first(p.spec_bands[1]),last(p.spec_bands[1]);length=n_spec))
    println("BENCH_BEGIN ",label); flush(stdout)
    forward_diagnostic = get(ENV,"AUDIT_FORWARD_DIAGNOSTIC","false") == "true"
    device_only = get(ENV,"AUDIT_DEVICE_PROFILE","false") == "true" || forward_diagnostic
    build_fwd=device_only ? nothing : measure(()->model_from_parameters(p))
    build_lin=device_only ? nothing : measure(()->model_from_parameters(LinMode(),p))
    build_fixed=device_only ? nothing : measure(()->model_from_parameters(LinMode(),p;
        compute_aerosol_microphysics_jacobians=false,compute_h2o_jacobians=false))
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
    if forward_diagnostic
        quiet(fwd); sync_backend()
        println("FORWARD_TIMER_DIAGNOSTIC ",label)
        with_logger(NullLogger()) do
            fwd(); sync_backend()
        end
        Profile.clear()
        Profile.@profile for _ in 1:3
            quiet(fwd); sync_backend()
        end
        open(joinpath(outdir,"profile-forward-$(label).txt"),"w") do io
            Profile.print(io;format=:flat,sortedby=:count,mincount=5,C=true)
        end
        println("FORWARD_DIAGNOSTIC_DONE ",label)
        return
    end
    if device_only
        quiet(full); sync_backend()
        result = device_profile(()->quiet(full)); sync_backend()
        show(stdout, MIME("text/plain"), result); println()
        println("DEVICE_PROFILE_DONE ",label)
        return
    end
    forward=measure(fwd); linear=measure(full); compact=measure(sel)
    rf=quiet(fwd); rl=quiet(full); rs=quiet(sel)
    native_columns=[1,2,7,CoreRT.n_total(native)]
    println("MODEL nSpec=",n_spec," nLayers=",length(m.profile.p_full),
        " operator=",length(m.quad_points.qp_μ)*p.polarization_type.n,
        " full_columns=",size(rl[3],4)," selected_columns=",size(rs[3],4))
    println("PARITY forward=",maximum(abs,rf[1]-rl[1]),
        " selected_radiance=",maximum(abs,rs[1]-rl[1]),
        " selected_jacobian=",maximum(abs,rs[3]-rl[3][:,:,:,native_columns]))
    for (name,result) in [("build_forward",build_fwd),("build_full",build_lin),
                          ("build_fixed_microphysics",build_fixed),
                          ("solve_forward",forward),("solve_full",linear),("solve_selected",compact)]
        println("MEASURE ",label," ",name," median_seconds=",median(x.time for x in result),
            " median_host_bytes=",median(x.bytes for x in result),
            " median_gc_seconds=",median(x.gctime for x in result)," samples=",result)
    end
    println("RATIO full_forward=",median(x.time for x in linear)/median(x.time for x in forward),
        " full_selected=",median(x.time for x in linear)/median(x.time for x in compact))
    flush(stdout)
    println("TIMER_BREAKDOWN_FULL ",label)
    with_logger(NullLogger()) do
        full(); sync_backend()
    end
    Profile.clear()
    Profile.@profile for _ in 1:3
        quiet(full); sync_backend()
    end
    open(joinpath(outdir,"profile-$(label).txt"),"w") do io
        Profile.print(io;format=:flat,sortedby=:count,mincount=5,C=true)
    end
    if backend == "cpu" && Threads.nthreads()==1
        Profile.Allocs.clear()
        Profile.Allocs.@profile sample_rate=0.005 quiet(full)
        allocs=Profile.Allocs.fetch().allocs
        totals=Dict{String,Tuple{Int,Int}}()
        for a in allocs
            frame=findfirst(s->occursin("vsmartmom-release-review/src/",string(s.file)),a.stacktrace)
            key=frame===nothing ? "outside candidate source" : string(a.stacktrace[frame])
            count,bytes=get(totals,key,(0,0))
            totals[key]=(count+1,bytes+a.size)
        end
        open(joinpath(outdir,"allocs-$(label).txt"),"w") do io
            println(io,"Sampled allocations grouped by nearest candidate source frame; sizes are sampled bytes, not full allocation totals.")
            for (key,(count,bytes)) in sort(collect(totals);by=x->last(x)[2],rev=true)
                println(io,bytes," bytes ",count," allocations | ",key)
            end
        end
    end
    println("BENCH_DONE ",label)
end
main()
