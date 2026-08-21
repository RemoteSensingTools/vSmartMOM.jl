using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.Architectures
using Statistics
using Test

const EXTERNAL_RAMAN_TEST_ARCH =
    get(ENV, "EXTERNAL_SOLAR_TEST_GPU", "0") == "1" ? GPU() : CPU()

function _small_rrs(model)
    FT = CoreRT.float_type(model)
    ν = model.atmosphere.spec_bands[1]
    ν̃ = mean(ν)
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)
    nPol = CoreRT.polarization_type(model).n
    F₀ = zeros(FT, nPol, length(ν)); F₀[1,:] .= one(FT)
    rs = InelasticScattering.RRS(
        n2=n2, o2=o2,
        greek_raman=InelasticScattering.GreekCoefs(
            [one(FT)], [one(FT)], [one(FT)], [one(FT)], [one(FT)], [one(FT)]),
        fscattRayl=[one(FT)], ϖ_Cabannes=[one(FT)],
        ϖ_λ₁λ₀=zeros(FT,1), i_λ₁λ₀=zeros(Int,1),
        Z⁻⁺_λ₁λ₀=zeros(FT,1,1), Z⁺⁺_λ₁λ₀=zeros(FT,1,1),
        i_ref=argmin(abs.(ν .- ν̃)), n_Raman=0,
        F₀=F₀, SIF₀=zeros(FT,nPol,length(ν)))
    CoreRT.getRamanSSProp!(rs, 1e7/ν̃, ν)
    rs
end

@testset "External-solar rotational Raman columns" begin
    p = parameters_from_yaml("test_parameters/Phase1b_RRS_761-764nm.yaml")
    p.architecture = EXTERNAL_RAMAN_TEST_ARCH
    legacy = model_from_parameters(p; external_solar=false)
    external = model_from_parameters(p)
    rs_legacy = _small_rrs(legacy)
    rs_external = _small_rrs(external)

    q = external.quad_points
    phase_qp_μ_host = Array(q.phase_qp_μ)
    qp_μ_host = Array(q.qp_μ)
    n = CoreRT.polarization_type(external).n
    solar = q.iμ₀Nstart_phase:(q.iμ₀Nstart_phase+n-1)
    rows = 1:length(q.qp_μN)
    for m in 0:2
        Zpp, Zmp = vSmartMOM.Scattering.compute_Z_moments(
            CoreRT.polarization_type(external), phase_qp_μ_host,
            rs_external.greek_raman, m)
        Z0pp, Z0mp = vSmartMOM.Scattering.compute_Z_source_moments(
            CoreRT.polarization_type(external), qp_μ_host, q.μ₀,
            rs_external.greek_raman, m)
        @test Z0pp ≈ Zpp[rows,solar] rtol=8e-6 atol=2e-7
        @test Z0mp ≈ Zmp[rows,solar] rtol=8e-6 atol=2e-7
    end

    old_R, _, old_ieR, _, _, _, _ = CoreRT.rt_run_test(rs_legacy, legacy, 1)
    new = CoreRT.rt_run_toa(rs_external, external; i_band=1)
    # The two paths contract in a different order (square-column slice versus
    # rectangular column), so Float32 roundoff is expected at a few 1e-7
    # relative even though the underlying operators are algebraically equal.
    @test new.elastic ≈ old_R rtol=2e-6 atol=2e-11
    @test new.inelastic ≈ old_ieR rtol=2e-6 atol=2e-11
    @test all(isfinite, new.inelastic)
end
