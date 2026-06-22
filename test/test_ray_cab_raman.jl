using Test
using vSmartMOM
using vSmartMOM.InelasticScattering

const _CAB_NM_PER_CM = 1e7
const _CAB_TEST_T = 300.0

_gamma_cabannes_eq4(γ_ray, ϖ) =
    0.5 * (2γ_ray * (2 + 3ϖ) - 3 * (1 - ϖ)) /
    ((3 + 2ϖ) - 4γ_ray * (1 - ϖ))

_phase_p(γ) = (1 - γ) / (1 + 2γ)
_phase_q(γ) = (1 - 3γ) / (1 + 2γ)

_sigma_rayleigh(mol, λ_nm) =
    mol.effCoeff.σ_Rayl_coeff * (_CAB_NM_PER_CM / λ_nm)^4

function _air_gamma_cabannes_from_components(λ_nm, n2, o2)
    numerator_base = 0.0
    gamma_term = 0.0
    for mol in (n2, o2)
        ϖ, γ = InelasticScattering.compute_γ_mol_Cabannes!(λ_nm, mol)[1:2]
        κ = ϖ * mol.effCoeff.σ_Rayl_coeff * (3 - 4γ) / (1 + 2γ)
        numerator_base += mol.vmr * κ
        gamma_term += mol.vmr * κ * γ / (3 - 4γ)
    end
    return 3 / (4 + numerator_base / gamma_term)
end

const _CAB_EXPECTED = Dict(
    395.0 => Dict(
        :N2 => (γ_cab = 0.002651990787, γ_ray = 0.010209333405,
                σ_ray_e26 = 1.801590745153, ϖ_cab = 0.975225331974),
        :O2 => (γ_cab = 0.007459756259, γ_ray = 0.027688455218,
                σ_ray_e26 = 1.589698665998, ϖ_cab = 0.935467219870),
        :air => (γ_cab = 0.003485077346, γ_ray = 0.013282065179,
                 σ_ray_e26 = 1.759212329322, ϖ_cab = 0.968039909644),
    ),
    760.0 => Dict(
        :N2 => (γ_cab = 0.002804804170, γ_ray = 0.010725167743,
                σ_ray_e26 = 0.125262856093, ϖ_cab = 0.974056187314),
        :O2 => (γ_cab = 0.008045467574, γ_ray = 0.029674930609,
                σ_ray_e26 = 0.108639976534, ϖ_cab = 0.931203072893),
        :air => (γ_cab = 0.003696898814, γ_ray = 0.014001882706,
                 σ_ray_e26 = 0.121938280182, ϖ_cab = 0.966420256516),
    ),
)

@testset "Rayleigh/Cabannes/RRS corrected gamma" begin
    γ_rrs = 3 / 4

    for λ_nm in (395.0, 760.0)
        n2, o2 = InelasticScattering.getRamanAtmoConstants(
            _CAB_NM_PER_CM / λ_nm, _CAB_TEST_T)

        for (label, mol) in ((:N2, n2), (:O2, o2))
            expected = _CAB_EXPECTED[λ_nm][label]
            ϖ, γ_cab, γ_ray = InelasticScattering.compute_γ_mol_Cabannes!(
                λ_nm, mol)

            @test γ_cab ≈ _gamma_cabannes_eq4(γ_ray, ϖ) rtol = 1e-13
            @test _phase_p(γ_ray) ≈
                  ϖ * _phase_p(γ_cab) + (1 - ϖ) * _phase_p(γ_rrs) rtol = 1e-13
            @test _phase_q(γ_ray) ≈
                  ϖ * _phase_q(γ_cab) + (1 - ϖ) * _phase_q(γ_rrs) rtol = 1e-13

            @test γ_cab ≈ expected.γ_cab atol = 5e-9
            @test γ_ray ≈ expected.γ_ray atol = 5e-9
            @test ϖ ≈ expected.ϖ_cab atol = 5e-9
            @test _sigma_rayleigh(mol, λ_nm) / 1e-26 ≈
                  expected.σ_ray_e26 rtol = 5e-6
        end

        expected_air = _CAB_EXPECTED[λ_nm][:air]
        γ_air_cab, ϖ_air = InelasticScattering.compute_γ_air_Cabannes!(
            λ_nm, n2, o2)
        γ_air_ray, σ_air = InelasticScattering.compute_γ_air_Rayleigh!(
            λ_nm, n2, o2)

        γ_air_cab_from_components = _air_gamma_cabannes_from_components(
            λ_nm, n2, o2)
        @test γ_air_cab ≈ γ_air_cab_from_components rtol = 1e-13
        @test γ_air_cab ≈ expected_air.γ_cab atol = 5e-9
        @test γ_air_ray ≈ expected_air.γ_ray atol = 5e-9
        @test ϖ_air ≈ expected_air.ϖ_cab atol = 5e-9
        @test σ_air / 1e-26 ≈ expected_air.σ_ray_e26 rtol = 5e-6
    end
end
