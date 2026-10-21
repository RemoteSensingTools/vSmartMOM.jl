include(ex -> ex == :(main()) ? nothing : ex,
        joinpath(@__DIR__, "local_basis_benchmark.jl"))
using vSmartMOM.Scattering
p = quiet(()->fixture_parameters(1))
scat = p.scattering_params
λ = 1e4/first(p.spec_bands[1])
mie = make_mie_model(scat.decomp_type,only(scat.rt_aerosols).aerosol,λ,
    p.polarization_type,CoreRT._resolved_truncation(p,Float64),scat.r_max,
    scat.nquad_radius; architecture=CPU())
for (name,f) in ((:forward,()->compute_aerosol_optical_properties(mie,Float64)),
                 (:linearized,()->compute_aerosol_optical_properties(LinMode(),mie,Float64)))
    f()
    results = [@timed(f()) for _ in 1:5]
    println("MIE_CPU mode=$name median=$(median(t.time for t in results)) bytes=$(median(t.bytes for t in results)) λ=$λ nquad=$(scat.nquad_radius)")
    Profile.clear()
    @profile for _ in 1:100; f(); end
    println("MIE_PROFILE mode=$name")
    Profile.print(stdout; format=:flat,C=false,sortedby=:count,mincount=10)
    flush(stdout)
end
