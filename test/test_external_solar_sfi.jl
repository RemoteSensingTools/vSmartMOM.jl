# Exact external-μ₀ SFI regression tests.
#
# The collimated solar direction is a source parameter, not a diffuse stream.
# These tests compare the historical zero-weight-μ₀ square operator against
# the opt-in external carrier.  The zero-weight block is mathematically
# decoupled, so the TOA Stokes field must agree while the generic IQU operator
# shrinks from 21×21 to 18×18 for five weighted streams plus one VZA. BOA,
# HDR, and BHR are intentionally outside this exoplanet output contract.

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Architectures
using Test

function _external_solar_model(model; nstreams=5)
    q = CoreRT.rt_set_streams(
        model.solver.quadrature_type, model.geometry,
        model.solver.polarization_type, array_type(model.architecture);
        nstreams, external_solar=true)
    return CoreRT.RTModel(
        model.architecture, model.solver, model.numerics,
        model.geometry, q, model.atmosphere, model.optics,
        model.surfaces, model.sources)
end

function _parity_model(; pressure_bar, albedo, sza, vza, vaz)
    p = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
    p.architecture = CPU()
    p.polarization_type = CoreRT.Stokes_IQU()
    p.nstreams = 5
    p.l_trunc = 9
    p.sza = Float64(sza)
    p.vza = Float64[vza]
    p.vaz = Float64[vaz]
    p.p = collect(range(0.01, Float64(pressure_bar) * 1000; length=6))
    p.T = fill(250.0, 5)
    p.q = zeros(5)
    p.profile_reduction_n = -1
    p.brdf = CoreRT.AbstractSurfaceType[
        CoreRT.LambertianSurfaceScalar(Float64(albedo))]
    return model_from_parameters(p)
end

@testset "External solar SFI (mu0 is not a stream)" begin
    cases = (
        (pressure_bar=0.1, albedo=0.0, sza=41.0, vza=23.0, vaz=117.0),
        (pressure_bar=10.0, albedo=1.0, sza=67.0, vza=38.0, vaz=53.0),
        (pressure_bar=1.0, albedo=0.3, sza=31.0, vza=31.0, vaz=180.0),
    )

    for case in cases
        legacy = _parity_model(; case...)
        external = _external_solar_model(legacy)
        q = external.quad_points

        @test q.external_solar
        @test q.iμ₀ == 0
        @test q.iμ₀Nstart == 0
        @test q.phase_qp_μ[q.iμ₀_phase] == q.μ₀
        @test !(q.μ₀ in q.qp_μ) || case.vza == case.sza

        if case.vza != case.sza
            @test legacy.quad_points.Nquad == 7
            @test external.quad_points.Nquad == 6
            @test length(legacy.quad_points.qp_μN) == 21
            @test length(external.quad_points.qp_μN) == 18
        end

        old = rt_run(legacy)
        new_toa = rt_run_toa(external)
        @test new_toa ≈ old.toa rtol=2e-13 atol=2e-15

        # The external-μ₀ representation is deliberately TOA-only; it must
        # never fall back to allocating BOA/HDR/BHR through the generic driver.
        @test_throws ArgumentError rt_run(external)
        @test_throws ArgumentError rt_run_ss(external)
    end
end
