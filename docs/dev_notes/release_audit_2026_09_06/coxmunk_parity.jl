using vSmartMOM, vSmartMOM.CoreRT, YAML, LinearAlgebra, Logging
BLAS.set_num_threads(1)
d=YAML.load_file(joinpath(pkgdir(vSmartMOM),"config","quickstart.yaml"))
d["radiative_transfer"]["surface"]=["CoxMunkSurface(wind_speed=5.0)"]
d["radiative_transfer"]["polarization_type"]="Stokes_IQU()"
d["geometry"]["sza"]=30.0
d["geometry"]["vza"]=[30.0]
d["geometry"]["vaz"]=[0.0]
p=read_parameters(d)
m,l=model_from_parameters(LinMode(),p)
with_logger(NullLogger()) do
    R=rt_run(m)[1]
    result=rt_run(m,l,0,size(l.τ̇_abs[1],1),1)
    println("COXMUNK_SAME_MODEL_FORWARD=",R)
    println("COXMUNK_SAME_MODEL_LINEARIZED_FORWARD=",result[1])
    println("COXMUNK_SAME_MODEL_RELATIVE_DEFECT=",maximum(abs,R-result[1])/maximum(abs,R))
end
