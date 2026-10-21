#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "Round4SIFTruthConvention.jl"))
using .Round4SIFTruthConvention

@testset "round-4 corrected-v2 SIF truth convention" begin
    convention = round4_sif_truth_convention()
    @test validate_round4_sif_truth_convention(convention) === convention

    at_759 = truth_sif_at_wavelength_nm(759.0)
    @test at_759.wavenumber_cm1 == 1.0e7 / 759.0
    @test isapprox(at_759.Lnu, 0.004818031987713776;
                   atol=3e-18, rtol=0)
    @test isapprox(at_759.Llambda, 0.0836346275560863;
                   atol=3e-17, rtol=0)
    @test isapprox(at_759.linear_Lnu, 0.004809473438302729;
                   atol=3e-18, rtol=0)
    @test isapprox(at_759.linear_Llambda, 0.08348606252076927;
                   atol=3e-17, rtol=0)
    @test at_759.Lnu != at_759.linear_Lnu

    at_760 = truth_sif_at_wavelength_nm(760.0)
    @test isapprox(at_760.Lnu, convention.SIF760;
                   atol=2e-18, rtol=0)
    @test isapprox(at_760.Llambda, 0.5 / (2pi);
                   atol=2e-16, rtol=0)
    @test isapprox(llambda_to_lnu(at_759.Llambda, 759.0), at_759.Lnu;
                   atol=3e-18, rtol=0)
    @test isapprox(lnu_to_llambda(at_759.Lnu, 759.0), at_759.Llambda;
                   atol=3e-17, rtol=0)
    @test_throws ArgumentError truth_sif_at_wavelength_nm(0.0)
end

@testset "round-4 corrected-v2 SIF provenance" begin
    on = expected_round4_sif_provenance(; enabled=true)
    checked = validate_round4_sif_provenance(
        on; enabled=true, source="synthetic SIF-on truth")
    @test Set(keys(checked)) == Set(keys(on))
    @test checked["sif_definition_version"] == Int32(2)
    @test checked["sif_reference_wavelength_nm"] == 760.0
    @test isapprox(
        checked["sif_radiance_760_mW_m-2_sr-1_nm-1"],
        0.5 / (2pi); atol=2e-16, rtol=0)
    @test isapprox(
        checked["sif_SIF760_mW_m-2_sr-1_per_cm-1"],
        0.004596394756493938; atol=2e-18, rtol=0)
    @test isapprox(
        checked["sif_mSIF_mW_m-2_sr-1_per_cm-2"],
        1.2291230681458325e-5; atol=2e-19, rtol=0)

    off = expected_round4_sif_provenance(; enabled=false)
    @test validate_round4_sif_provenance(
        off; enabled=false, source="synthetic SIF-off truth") == off
    @test iszero(off["sif_angular_integral_760_mW_m-2_nm-1"])
    @test iszero(off["sif_SIF760_mW_m-2_sr-1_per_cm-1"])
    @test iszero(off["sif_mSIF_mW_m-2_sr-1_per_cm-2"])

    stale = copy(on)
    stale["sif_case_on_label"] = "total_0p5"
    @test_throws ErrorException validate_round4_sif_provenance(
        stale; enabled=true, source="legacy truth")

    wrong_units = copy(on)
    wrong_units["sif_SIF760_mW_m-2_sr-1_per_cm-1"] =
        wrong_units["sif_radiance_760_mW_m-2_sr-1_nm-1"]
    @test_throws ErrorException validate_round4_sif_provenance(
        wrong_units; enabled=true, source="unit-corrupted truth")

    missing = copy(on)
    delete!(missing, "sif_definition_version")
    @test_throws ErrorException validate_round4_sif_provenance(
        missing; enabled=true, source="incomplete truth")
end
