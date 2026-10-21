using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using YAML, LinearAlgebra, Statistics
BLAS.set_num_threads(1)
const ROOT = pkgdir(vSmartMOM)
function probe(f, name)
    println("\nAUDIT_PROBE_BEGIN ", name); flush(stdout)
    try
        f()
        println("AUDIT_PROBE_COMPLETE ", name)
    catch err
        println("AUDIT_PROBE_ERROR ", name, " ", sprint(showerror, err, catch_backtrace()))
    end
    flush(stdout)
end
probe("quickstart and documented external-solar default") do
    p = read_parameters(joinpath(ROOT,"config","quickstart.yaml"))
    m = model_from_parameters(p)
    R,T = rt_run(m)
    println("quickstart R=",R," T=",T," external_solar=",m.quad_points.external_solar)
    try
        rt_run_toa(m)
    catch err
        println("default model rt_run_toa: ",sprint(showerror,err))
    end
    me = model_from_parameters(p; external_solar=true)
    Re = rt_run_toa(me)
    println("explicit external TOA finite=",all(isfinite,Re)," R=",Re)
end
probe("Cox-Munk source linearity") do
    d = YAML.load_file(joinpath(ROOT,"config","quickstart.yaml"))
    d["radiative_transfer"]["surface"] = ["CoxMunkSurface(wind_speed=5.0)"]
    d["radiative_transfer"]["polarization_type"] = "Stokes_IQU()"
    d["geometry"]["sza"] = 30.0
    d["geometry"]["vza"] = [30.0]
    d["geometry"]["vaz"] = [0.0]
    p=read_parameters(d); m=model_from_parameters(p)
    F=zeros(3,1); F[1,:].=1
    R1=rt_run(m;sources=SolarBeam(F₀=F))[1]
    R2=rt_run(m;sources=SolarBeam(F₀=2F))[1]
    R0=rt_run(m;sources=NoSource())[1]
    println("R1=",R1," R2=",R2," dark R0=",R0)
    println("solar doubling relative defect=",maximum(abs,R2-2R1)/maximum(abs,2R1))
end
probe("reference-index forward/linearized optics parity") do
    d=YAML.load_file(joinpath(ROOT,"test","test_parameters","JacobianTestFast.yaml"))
    delete!(d,"absorption")
    d["scattering"]["n_ref"]="1.5 - 0.0im"
    d["scattering"]["aerosols"][1]["μ"]=0.15
    d["scattering"]["aerosols"][1]["σ"]=1.4
    d["scattering"]["r_max"]=3.0
    p=read_parameters(d)
    mf=model_from_parameters(p)
    ml,l=model_from_parameters(LinMode(),p)
    ms,ls=model_from_parameters(LinMode(),p;compute_aerosol_microphysics_jacobians=false)
    println("AOD forward=",vec(sum(mf.τ_aer[1];dims=3))," lin=",vec(sum(ml.τ_aer[1];dims=3))," fixed_microphysics=",vec(sum(ms.τ_aer[1];dims=3)))
    println("AOD relative defect=",maximum(abs,ml.τ_aer[1]-mf.τ_aer[1])/maximum(abs,mf.τ_aer[1]))
    Rf=rt_run(mf)[1]; Rl=rt_run(ml)[1]
    println("radiance relative defect=",maximum(abs,Rl-Rf)/maximum(abs,Rf))
end
probe("README ocean linearized snippet") do
    p=read_parameters(joinpath(ROOT,"config","ocean_coxmunk.yaml"))
    println("scattering_params=",p.scattering_params)
    println("README NAer=",length(p.scattering_params.rt_aerosols))
end
probe("retired delta-BGE angle") do
    println("requested angle=2.0; actual=",δBGE(30,2.0).Δ_angle)
end
