using vSmartMOM, vSmartMOM.CoreRT, Test
p = read_parameters(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
m = model_from_parameters(p)
r = rt_run(m; sources=SolarBeam()).toa
r_sza = rt_run(m; sources=SolarBeam(sza=75)).toa
r_two = rt_run(m; sources=SolarBeam()+SolarBeam()).toa
r_double = rt_run(m; sources=SolarBeam(F₀=fill(2.0, 1, 1))).toa
println("sza ignored: ", r == r_sza)
println("two beams match one: ", r_two == r)
println("two beams / summed irradiance: ", r_two[1]/r_double[1])
p.brdf[1] = CoreRT.rpvSurfaceScalar(0.1, 1.0, 0.0, 0.1)
m = model_from_parameters(p)
r = rt_run(m; sources=SolarBeam()).toa
r_sif = rt_run(m; sources=SolarBeam()+SurfaceSIF(SIF₀=fill(0.1, 1, 1))).toa
println("RPV prescribed SIF ignored: ", r == r_sif)
