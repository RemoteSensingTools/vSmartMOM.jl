#!/usr/bin/env julia
# ==========================================================================
# GPU-only Raman Forward Model Test
# ==========================================================================
# Runs the RRS forward model specifically on GPU with the small-grid
# O2Parameters_GPU.yaml config (~60 spectral points).
# Mirrors test_forward_raman.jl but forces GPU execution and skips
# if CUDA is not available.
# ==========================================================================

using Test
using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Statistics

CUDA_OK = false
try
    using CUDA
    global CUDA_OK = CUDA.functional()
catch
end

@testset "Raman GPU Forward" begin

if !CUDA_OK
    @test_skip "CUDA not available - skipping Raman GPU tests"
else

    params = parameters_from_yaml("test_parameters/O2Parameters_GPU.yaml")
    params.architecture = vSmartMOM.Architectures.GPU()
    model = model_from_parameters(params)

    iBand = 1
    FT = Float64
    ν = model.params.spec_bands[iBand]
    nSpec = length(ν)
    ν̃ = mean(ν)

    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

    nPol = model.params.polarization_type.n
    F₀ = zeros(FT, nPol, nSpec); F₀[1, :] .= 1.0
    SIF₀ = zeros(FT, nPol, nSpec)

    RS_type = InelasticScattering.RRS(
        n2          = n2,
        o2          = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)],
        ϖ_λ₁λ₀     = zeros(FT, 1),
        i_λ₁λ₀     = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀   = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀   = zeros(FT, 1, 1),
        i_ref       = argmin(abs.(ν .- ν̃)),
        n_Raman     = 0,
        F₀          = F₀,
        SIF₀        = SIF₀
    )
    CoreRT.getRamanSSProp!(RS_type, 1e7/ν̃, ν)

    println("  Running Raman forward on GPU (nSpec=$nSpec)...")

    @testset "RRS forward GPU" begin
        R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(RS_type, model, iBand)

        @test ndims(R_rrs) == 3
        @test all(isfinite.(R_rrs))
        @test all(isfinite.(ieR))

        nVza = length(model.obs_geom.vza)
        @test size(R_rrs) == (nVza, nPol, nSpec)

        I_elastic = R_rrs[:, 1, :]
        @test all(I_elastic .> 0) || @warn "Some elastic I < 0: min=$(minimum(I_elastic))"
        println("    R shape: $(size(R_rrs)), all finite and positive.")
    end

    @testset "noRS GPU baseline" begin
        noRS_type = InelasticScattering.noRS(
            fscattRayl  = [FT(1)],
            ϖ_Cabannes  = [FT(1)],
            bandSpecLim = UnitRange{Int}[],
            iBand       = [1],
            F₀          = zeros(FT, nPol, nSpec)
        )
        noRS_type.F₀[1, :] .= 1.0

        R_noRS, T_noRS, _, _ = CoreRT.rt_run_test(noRS_type, model, iBand)

        @test all(isfinite.(R_noRS))
        I_noRS = R_noRS[:, 1, :]
        @test all(I_noRS .> 0)
        println("    noRS baseline on GPU: all finite and positive.")
    end

    @testset "Ring effect GPU" begin
        noRS_type = InelasticScattering.noRS(
            fscattRayl  = [FT(1)],
            ϖ_Cabannes  = [FT(1)],
            bandSpecLim = UnitRange{Int}[],
            iBand       = [1],
            F₀          = zeros(FT, nPol, nSpec)
        )
        noRS_type.F₀[1, :] .= 1.0

        R_noRS, _, _, _ = CoreRT.rt_run_test(noRS_type, model, iBand)
        R_rrs, _, ieR, _ = CoreRT.rt_run_test(RS_type, model, iBand)

        I_total = R_rrs[:, 1, :] .+ ieR[:, 1, :]
        I_noRS  = R_noRS[:, 1, :]

        ring_pct = 100.0 .* (I_total .- I_noRS) ./ I_noRS
        mean_ring = mean(abs.(ring_pct))
        max_ring  = maximum(abs.(ring_pct))

        @test mean_ring < 50.0
        @test max_ring < 200.0
        @test all(isfinite.(I_total))
        @test all(I_total .> 0)
        println("    Ring effect: mean=$(round(mean_ring, digits=2))%, max=$(round(max_ring, digits=2))%")
    end

end # CUDA_OK

end # testset
println("\nRaman GPU tests complete.")
