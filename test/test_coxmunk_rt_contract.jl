using Test, LinearAlgebra, YAML
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering

@testset "Cox–Munk radiance and Jacobian contract" begin
    d = YAML.load_file(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
    d["radiative_transfer"]["surface"] = ["CoxMunkSurface(wind_speed=5.0)"]
    d["radiative_transfer"]["polarization_type"] = "Stokes_IQU()"
    d["radiative_transfer"]["spec_bands"] = ["[12986.0, 12987.0, 12988.0]"]
    d["geometry"]["sza"] = 30.0
    d["geometry"]["vza"] = [30.0, 45.0]
    d["geometry"]["vaz"] = [0.0, 35.0]
    d["radiative_transfer"]["architecture"] = get(ENV, "VSMARTMOM_COXMUNK_ARCH", "CPU()")
    p = read_parameters(d)
    model, lin = model_from_parameters(LinMode(), p)
    ngas = size(lin.τ̇_abs[1], 1)
    F = [1.0 2.0 0.5; 0.0 0.0 0.0; 0.0 0.0 0.0]
    source = SolarBeam(F₀=F)
    R = rt_run(model; sources=source)[1]
    R2 = rt_run(model; sources=SolarBeam(F₀=2F))[1]
    @test R2 ≈ 2R rtol=2e-12
    @test iszero(rt_run(model; sources=NoSource())[1])
    linearized = rt_run(model, lin, 0, ngas, 1; sources=source)
    @test linearized[1] ≈ R rtol=2e-11 atol=2e-14
    @test iszero(rt_run(model, lin, 0, ngas, 1; sources=NoSource())[1])
    cache = rt_run_atmosphere(model; sources=source)
    @test rt_run_surface(cache, model.surfaces[1])[1] == R

    # Analytic wind derivative against independent complete forward solves.
    original = model.surfaces[1]
    h = 1e-4
    model.surfaces[1] = CoxMunkSurface(wind_speed=5.0+h)
    plus = rt_run(model; sources=source)[1]
    model.surfaces[1] = CoxMunkSurface(wind_speed=5.0-h)
    minus = rt_run(model; sources=source)[1]
    model.surfaces[1] = original
    wind = first(CoreRT.surface_range(linearized.layout))
    @test linearized[3][:,:,:,wind] ≈ (plus-minus)/(2h) rtol=2e-5 atol=2e-10

    # In a purely absorbing atmosphere the exact direct surface path is the
    # entire TOA field. This checks absolute BRDF normalization, Fourier m=0,
    # and both atmospheric transits, rather than merely solver-to-solver parity.
    fill!(model.τ_rayl[1], 0)
    for optical_depth in (0.0, 0.3)
        fill!(model.τ_abs[1], optical_depth / size(model.τ_abs[1], 2))
        actual = rt_run(model; sources=source)[1]
        expected = similar(actual)
        for iv in eachindex(p.vza), s in axes(F,2)
            μv, μ0 = cosd(p.vza[iv]), cosd(p.sza)
            M = CoreRT.coxmunk_brdf_mueller(original, 3, μv, μ0, deg2rad(p.vaz[iv]))
            expected[iv,:,s] = μ0 * F[1,s] * M[:,1] *
                              exp(-optical_depth * (inv(μ0) + inv(μv)))
        end
        @test actual ≈ expected rtol=2e-10 atol=2e-14
    end

    # Zero Fresnel contrast isolates the whitecap Lambertian term and locks
    # the π*BRDF convention required by the diffuse reflection operator.
    whitecap = CoxMunkSurface(wind_speed=12.0, n_water=1.0+0im,
                              include_whitecaps=true)
    albedo = CoreRT.whitecap_fraction(12.0) * whitecap.whitecap_albedo
    M0 = CoreRT.reflectance(whitecap, Stokes_I(), [0.4, 0.8], 0)
    @test M0 ≈ fill(albedo, 2, 2) rtol=1e-12

    # Optical-depth tangent of the same correction, without differentiating
    # through an implementation copy or loosening the forward tolerance.
    τ = [0.2, 0.3, 0.4]
    τdot = reshape([0.5, 0.8, 1.2], 3, 1)
    Rc = zeros(2,3,3); dRc = zeros(2,3,3,2)
    CoreRT.apply_ss_correction!(Rc, dRc, τdot, 2, original, Stokes_IQU(),
        p.vza, p.vaz, cosd(p.sza), τ, 5, 3; F₀=F)
    Rp = zero(Rc); Rm = zero(Rc)
    CoreRT.apply_ss_correction!(Rp, original, Stokes_IQU(), p.vza, p.vaz,
        cosd(p.sza), τ + h*vec(τdot), 5, 3; F₀=F)
    CoreRT.apply_ss_correction!(Rm, original, Stokes_IQU(), p.vza, p.vaz,
        cosd(p.sza), τ - h*vec(τdot), 5, 3; F₀=F)
    @test dRc[:,:,:,1] ≈ (Rp-Rm)/(2h) rtol=2e-7 atol=1e-12
end
