# Reproduce the initial absorption-free benchmark's inputs and timing boundary.
include(ex -> ex == :(main()) ? nothing : ex,
    joinpath(@__DIR__, "benchmark.jl"))
using TOML
@assert backend == "cuda" "Set AUDIT_BACKEND=cuda for this historical comparison"
CUDA.allowscalar(false)
records = []
for pol in ("Stokes_I()","Stokes_IQU()")
    d=YAML.load_file(joinpath(pkgdir(vSmartMOM),"test/test_parameters/JacobianTestFast.yaml"))
    delete!(d,"absorption")
    d["radiative_transfer"]["architecture"]="GPU()"
    d["radiative_transfer"]["numerics"]=Dict("verbose"=>true,"blas_threads"=>1)
    d["radiative_transfer"]["greek_beta_cutoff"]=nothing
    d["radiative_transfer"]["polarization_type"]=pol
    d["scattering"]["r_max"]=3.0
    d["scattering"]["aerosols"][1]["μ"]=0.15
    d["scattering"]["aerosols"][1]["σ"]=1.4
    p=read_parameters(d)
    p.spec_bands[1]=collect(range(first(p.spec_bands[1]),last(p.spec_bands[1]);length=n_spec))
    println("PREPARE original scene polarization=$pol"); flush(stdout)
    m,l=quiet(()->model_from_parameters(LinMode(),p))
    ng=size(l.τ̇_abs[1],1)
    forward=()->rt_run(m)
    source=()->rt_run(m,l,1,ng,1;jacobian_basis=:local,jacobian_adding=:source)
    baseline=nothing
    for (name,f) in (("source",source),("forward",forward))
        GC.gc(); CUDA.reclaim()
        println("MEASURE_BEGIN polarization=$pol mode=$name"); flush(stdout)
        ts=measure(f)
        out=quiet(f); sync_backend()
        err=0.0
        if name=="source"
            baseline=out
            @assert size(out[3],4)==14
        else
            for i in 1:2
                @assert isapprox(out[i],baseline[i];rtol=1e-9,atol=1e-11)
                err=max(err,maximum(abs,out[i]-baseline[i]))
            end
        end
        record=Dict("polarization"=>pol,"mode"=>name,"seconds"=>[t.time for t in ts],
            "median_seconds"=>median(t.time for t in ts),"max_forward_parity_error"=>err,
            "host_bytes"=>[t.bytes for t in ts],"device_bytes"=>[t.gpu_bytes for t in ts])
        push!(records,record)
        println("RESULT ",record); flush(stdout)
    end
end
open(io->TOML.print(io,Dict("measurements"=>records,"commit"=>readchomp(`git -C $(pkgdir(vSmartMOM)) rev-parse HEAD`),
    "nSpec"=>n_spec,"layers"=>5,"columns"=>14,"absorption"=>false,
    "construction_timed"=>false,"chunking"=>"single spectral batch")),
    get(ENV,"AUDIT_OUTPUT","/tmp/original-fixture-current.toml"),"w")
