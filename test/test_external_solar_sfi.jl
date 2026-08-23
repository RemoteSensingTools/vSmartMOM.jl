# Exact external-μ₀ SFI regression tests.
#
# The collimated solar direction is a source parameter, not a diffuse stream.
# These tests compare the explicitly requested historical zero-weight-μ₀
# square operator against the default external carrier. The zero-weight block is mathematically
# decoupled, so the TOA Stokes field must agree while the generic IQU operator
# shrinks from 21×21 to 18×18 for five weighted streams plus one VZA. BOA,
# HDR, and BHR are intentionally outside this exoplanet output contract.

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Architectures
using vSmartMOM.Scattering
using Test

# Qualified: the shared runtests session also loads AtmosphericAbsorption,
# whose Architectures module exports CPU/GPU too — the bare names are
# ambiguous in Main.
const EXTERNAL_SOLAR_TEST_ARCH =
    get(ENV, "EXTERNAL_SOLAR_TEST_GPU", "0") == "1" ?
    vSmartMOM.Architectures.GPU() : vSmartMOM.Architectures.CPU()

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

function _parity_parameters(; pressure_bar, albedo, sza, vza, vaz)
    p = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
    p.architecture = EXTERNAL_SOLAR_TEST_ARCH
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
    return p
end

@testset "External solar SFI (mu0 is not a stream)" begin
    cases = (
        (pressure_bar=0.1, albedo=0.0, sza=41.0, vza=23.0, vaz=117.0),
        (pressure_bar=10.0, albedo=1.0, sza=67.0, vza=38.0, vaz=53.0),
        (pressure_bar=1.0, albedo=0.3, sza=31.0, vza=31.0, vaz=180.0),
    )

    for case in cases
        p = _parity_parameters(; case...)
        legacy = model_from_parameters(p; external_solar=false)
        external = model_from_parameters(p; external_solar=true)
        q = external.quad_points

        @test q.external_solar
        @test q.iμ₀ == 0
        @test q.iμ₀Nstart == 0
        phase_qp_μ_host = Array(q.phase_qp_μ)
        qp_μ_host = Array(q.qp_μ)
        @test phase_qp_μ_host[q.iμ₀_phase] == q.μ₀
        @test !(q.μ₀ in qp_μ_host) || case.vza == case.sza

        # A directly evaluated rectangular solar block must equal the exact
        # μ₀ columns sliced from the historical augmented square matrix.
        greek = external.optics.rayleigh.greek_rayleigh[1]
        nstokes = p.polarization_type.n
        solar_cols = q.iμ₀Nstart_phase:(q.iμ₀Nstart_phase + nstokes - 1)
        diffuse_rows = 1:length(q.qp_μN)
        for m in 0:2
            Zpp_aug, Zmp_aug = Scattering.compute_Z_moments(
                p.polarization_type, phase_qp_μ_host, greek, m)
            Z0pp, Z0mp = Scattering.compute_Z_source_moments(
                p.polarization_type, qp_μ_host, q.μ₀, greek, m)
            @test Z0pp ≈ Zpp_aug[diffuse_rows, solar_cols] rtol=5e-15 atol=5e-16
            @test Z0mp ≈ Zmp_aug[diffuse_rows, solar_cols] rtol=5e-15 atol=5e-16

            # The phase mapping is linear in every Greek coefficient.  Use
            # nonzero, parameter-distinct tangents to verify that the direct
            # rectangular tangent has exactly the same layout and values as
            # the corresponding columns of the legacy augmented derivative.
            ncoef = length(greek.β)
            tangent = reshape(collect(1.0:(4ncoef)), 4, ncoef) ./ (10ncoef)
            lin_greek = Scattering.linGreekCoefs(
                tangent, 2tangent, 3tangent, 4tangent, 5tangent, 6tangent)
            _, _, dZpp_aug, dZmp_aug = Scattering.compute_Z_moments(
                p.polarization_type, phase_qp_μ_host, greek, lin_greek, m)
            _, _, dZ0pp, dZ0mp = Scattering.compute_Z_source_moments(
                p.polarization_type, qp_μ_host, q.μ₀, greek, lin_greek, m)
            @test dZ0pp ≈ dZpp_aug[:, diffuse_rows, solar_cols] rtol=7e-15 atol=7e-16
            @test dZ0mp ≈ dZmp_aug[:, diffuse_rows, solar_cols] rtol=7e-15 atol=7e-16
        end

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
        if case == first(cases)
            legacy_fwd, legacy_lin = model_from_parameters(
                LinMode(), p; external_solar=false)
            external_fwd, external_lin = model_from_parameters(LinMode(), p;
                                                               external_solar=true)
            ngas = size(legacy_lin.τ̇_abs[1], 1)
            old_lin = rt_run(legacy_fwd, legacy_lin, 0, ngas, 1; i_band=1)
            new_lin = rt_run(external_fwd, external_lin, 0, ngas, 1; i_band=1)
            @test new_lin.toa ≈ old_lin.toa rtol=3e-13 atol=3e-15
            @test new_lin.toa_jacobian ≈ old_lin.toa_jacobian rtol=2e-11 atol=2e-13
        end
    end
end
