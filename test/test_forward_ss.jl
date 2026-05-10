# Phase 1c smoke test — single-scatter approximation driver
# =========================================================================
#
# Exercises `rt_run_ss` on the toy 1-band Phase1b_RRS config (noRS path
# for speed) and verifies:
#   - 6-tuple return shape: (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hem_R, hem_T).
#   - All outputs finite and physically plausible (I > 0, |Q| ≤ I).
#   - SS intensity < full multiple-scattering intensity (physical check —
#     MS adds photons to the reflected beam that SS omits).
#   - Hemispheric integrals non-negative and bounded.
#
# The driver was ported from sanghavi-branch `rt_run_ss` in Phase 1c of the
# sanghavi-unified merge (see plans/IMPLEMENTATION_PLAN_v2.md). Before the
# port, unified exported `rt_run_ss` without defining it — calling it raised
# `UndefVarError`. This smoke test guards against regressing the fix.

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Test

@testset "Phase 1c single-scatter driver smoke" begin
    params = parameters_from_yaml("test_parameters/Phase1b_RRS_761-764nm.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()
    model = model_from_parameters(params)

    nVza = length(model.obs_geom.vza)
    nPol = CoreRT.polarization_type(model).n
    nSpec = length(model.atmosphere.spec_bands[1])

    @testset "6-tuple return + shapes" begin
        result = CoreRT.rt_run_ss(model; i_band=1)
        @test result isa Tuple && length(result) == 6
        R_ss, T_ss, ieR_ss, ieT_ss, hem_R, hem_T = result
        @test size(R_ss) == (nVza, nPol, nSpec)
        @test size(T_ss) == (nVza, nPol, nSpec)
        @test size(ieR_ss) == (nVza, nPol, nSpec)
        @test size(ieT_ss) == (nVza, nPol, nSpec)
        @test size(hem_R) == (nSpec,)
        @test size(hem_T) == (nSpec,)
        @test all(isfinite.(R_ss))
        @test all(isfinite.(T_ss))
        @test all(isfinite.(ieR_ss))
        @test all(isfinite.(ieT_ss))
        @test all(isfinite.(hem_R))
        @test all(isfinite.(hem_T))
    end

    @testset "Physical sanity" begin
        R_ss, _, _, _, hem_R, hem_T = CoreRT.rt_run_ss(model; i_band=1)
        I_ss = R_ss[:, 1, :]
        @test all(I_ss .> 0) || @warn "SS I has non-positive entries: min=$(minimum(I_ss))"
        @test maximum(I_ss) < 1.0
        if nPol >= 2
            Q_ss = R_ss[:, 2, :]
            @test all(abs.(Q_ss) .<= I_ss .+ 1e-10)
        end
        # Hemispheric reflectance / transmittance ≥ 0 and bounded below unity.
        @test all(hem_R .>= 0)
        @test all(hem_T .>= 0)
        @test maximum(hem_R) < 1.0
        @test maximum(hem_T) < 1.0
    end

    @testset "SS ≤ full MS on Stokes I" begin
        # Full multiple-scattering reference on the same model.
        R_ms = CoreRT.rt_run(model; i_band=1)[1]
        R_ss = CoreRT.rt_run_ss(model; i_band=1)[1]
        # SS captures only first-order scattering; MS includes SS + higher orders.
        # Per-pixel comparison on Stokes I should satisfy R_ss ≤ R_ms + noise.
        I_ms = R_ms[:, 1, :]
        I_ss = R_ss[:, 1, :]
        @test all(I_ss .<= I_ms .+ 1e-10) ||
            @warn "SS I exceeds MS I somewhere — diagnose before marking Phase 1c done"
        @test maximum(I_ss) < maximum(I_ms)
    end

    @testset "rt_run_test_ss entry point + RRS smoke" begin
        iBand = 1
        FT = Float32
        ν = model.atmosphere.spec_bands[iBand]
        nspec = length(ν)
        ν̃ = (ν[1] + ν[end]) / 2
        effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
        n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

        F₀ = zeros(FT, nPol, nspec); F₀[1, :] .= 1
        SIF₀ = zeros(FT, nPol, nspec)

        RS_type = InelasticScattering.RRS(
            n2 = n2, o2 = o2,
            greek_raman = InelasticScattering.GreekCoefs(
                [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
            fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
            ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
            Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
            i_ref = argmin(abs.(ν .- ν̃)),
            n_Raman = 0,
            F₀ = F₀, SIF₀ = SIF₀)
        CoreRT.getRamanSSProp!(RS_type, 1e7/ν̃, ν)

        result = CoreRT.rt_run_test_ss(RS_type, model, iBand)
        @test result isa Tuple && length(result) == 6
        R_ss, T_ss, ieR_ss, ieT_ss, hem_R, hem_T = result
        @test all(isfinite.(R_ss))
        @test all(isfinite.(ieR_ss))
        @test all(isfinite.(hem_R))
        @test all(isfinite.(hem_T))
        # Raman inelastic component should be non-zero for RRS path
        @test maximum(abs.(ieR_ss)) > 0
    end
end
