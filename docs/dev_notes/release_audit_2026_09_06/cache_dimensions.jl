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
    build_fwd=nothing
    build_lin=nothing
    m,l=quiet(()->model_from_parameters(LinMode(),p))
    beta=m.aerosol_optics[1][1].greek_coefs.β
    println("GREEK_BETA shape=",size(beta)," length=",length(beta)," angular_rows=",size(beta,1))
    println("TABLE_BASE_BYTES actual=",6*length(m.quad_points.phase_qp_μ)*length(beta)^2*sizeof(Float64),
        " corrected=",6*length(m.quad_points.phase_qp_μ)*size(beta,1)^2*sizeof(Float64))
    on=quiet(()->rt_run(m)); sync_backend()
    old=CoreRT._Z_TABLES_ENABLED[]
    try
        CoreRT._Z_TABLES_ENABLED[]=false
        timings=measure(()->rt_run(m))
        off=quiet(()->rt_run(m)); sync_backend()
        println("CACHE_OFF median_seconds=",median(t.time for t in timings),
            " median_host_bytes=",median(t.bytes for t in timings)," samples=",timings)
        println("CACHE_PARITY R=",maximum(abs,on[1]-off[1])," T=",maximum(abs,on[2]-off[2]))
    finally
        CoreRT._Z_TABLES_ENABLED[]=old
    end
    println("CACHE_PROBE_DONE ",label)
end
main()
