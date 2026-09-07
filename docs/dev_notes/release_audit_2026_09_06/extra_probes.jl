using vSmartMOM, vSmartMOM.CoreRT, YAML, LinearAlgebra, Statistics, Logging
BLAS.set_num_threads(1)
d=YAML.load_file(joinpath(pkgdir(vSmartMOM),"test","test_parameters","JacobianTestFast.yaml"))
delete!(d,"absorption")
d["scattering"]["r_max"]=3.0
d["scattering"]["aerosols"][1]["μ"]=0.15
d["scattering"]["aerosols"][1]["σ"]=1.4
push!(d["scattering"]["aerosols"],deepcopy(d["scattering"]["aerosols"][1]))
d["scattering"]["aerosols"][2]["nᵣ"]=1.5
p=read_parameters(d)
mf=model_from_parameters(p)
ml,l=model_from_parameters(LinMode(),p)
println("DEFAULT_MULTI_AER n_ref=",p.scattering_params.n_ref)
println("AOD_FWD=",sum(mf.τ_aer[1];dims=3))
println("AOD_LIN=",sum(ml.τ_aer[1];dims=3))
na=length(p.scattering_params.rt_aerosols); ng=size(l.τ̇_abs[1],1)
with_logger(NullLogger()) do
    Rf=rt_run(mf)[1]
    Rl=rt_run(ml)[1]
    rl=rt_run(ml,l,na,ng,1)
    println("MULTI_AER_RADIANCE_RELATIVE_DEFECT=",maximum(abs,Rl-Rf)/maximum(abs,Rf))
    println("SAME_MODEL_FWD_LIN_RELATIVE_DEFECT=",maximum(abs,rl[1]-Rl)/maximum(abs,Rl))
    println("JACOBIAN_SHAPE=",size(rl[3]))
    forward_times=[@elapsed(rt_run(ml)) for _ in 1:3]
    linear_times=[@elapsed(rt_run(ml,l,na,ng,1)) for _ in 1:3]
    println("FORWARD_SECONDS=",forward_times," LINEARIZED_SECONDS=",linear_times)
    println("WARM_SOLVE_RATIO=",median(linear_times)/median(forward_times))
end
