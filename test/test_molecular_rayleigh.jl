using Test
using vSmartMOM

# Regression anchors for the existing pure-gas RI helper. The HWO provenance
# review is scoped to Ar and CO2; He/N2/O2 remain compatibility behavior.
@testset "molecular Rayleigh species" begin
    ri_method = RefractiveIndexRayleigh()
    species = (:He, :Ar, :N2, :O2, :CO2)
    expected_sigma_400 = (
        He = 2.4579195071330603e-28,
        Ar = 1.4732404564376880e-26,
        N2 = 1.7138759042453256e-26,
        O2 = 1.5190047653546414e-26,
        CO2 = 4.3523045188182340e-26,
    )
    expected_ratio_360_400 = (
        He = 1.5349157992842803,
        Ar = 1.5492964716715354,
        N2 = 1.5512982516210456,
        O2 = 1.5661864589978960,
        CO2 = 1.5593995107355967,
    )

    for gas in species
        p360 = molecular_rayleigh_properties(ri_method, gas, 360.0)
        p400 = molecular_rayleigh_properties(ri_method, gas, 400.0)
        @test p400.species == gas
        # The compatibility API is still RI based, but new science code must
        # make the method visible at the call site and in provenance.
        @test p360 == molecular_rayleigh_properties(gas, 360.0)
        @test p400 == molecular_rayleigh_properties(gas, 400.0)
        # The implementation uses cancellation- and overflow-safe algebra for
        # both Float64 and Float32. Algebraically equivalent evaluation orders
        # differ by a few parts in 10⁻¹² from the original pinned values.
        @test p400.cross_section_cm2 ≈ getproperty(expected_sigma_400, gas) rtol = 1e-11
        explicit_ratio = molecular_rayleigh_cross_section_ratio(
            ri_method, gas, 360.0, 400.0)
        @test explicit_ratio ==
              molecular_rayleigh_cross_section_ratio(gas, 360.0, 400.0)
        @test explicit_ratio ≈
              getproperty(expected_ratio_360_400, gas) rtol = 1e-11
        @test p360.cross_section_cm2 > p400.cross_section_cm2 > 0
        @test 0 ≤ p400.depolarization < 1

        p400_f32 = molecular_rayleigh_properties(ri_method, gas, 400.0f0)
        ratio_f32 = molecular_rayleigh_cross_section_ratio(
            ri_method, gas, 360.0f0, 400.0f0)
        @test p400_f32.cross_section_cm2 isa Float32
        @test isfinite(p400_f32.cross_section_cm2)
        @test p400_f32.cross_section_cm2 > 0.0f0
        @test ratio_f32 isa Float32
        @test isfinite(ratio_f32)
        @test ratio_f32 > 1.0f0
        @test p400_f32.cross_section_cm2 ≈
              Float32(getproperty(expected_sigma_400, gas)) rtol = 5f-6
        @test ratio_f32 ≈
              Float32(getproperty(expected_ratio_360_400, gas)) rtol = 5f-6
    end

    @test molecular_rayleigh_properties(:He, 400.0).depolarization == 0
    @test molecular_rayleigh_properties(:Ar, 400.0).depolarization == 0
    @test molecular_rayleigh_properties(:N₂, 400.0).species === :N2
    @test molecular_rayleigh_properties("CO2", 400.0).species === :CO2
    @test molecular_rayleigh_properties(:N2, 400.0f0).wavelength_nm isa Float32
    @test_throws ArgumentError molecular_rayleigh_properties(:H2O, 400.0)
    @test_throws ArgumentError molecular_rayleigh_properties(:N2, 0.0)

    for gas in (:Ar, :CO2)
        provenance = molecular_rayleigh_provenance(ri_method, gas)
        @test provenance.method === :refractive_index_king_factor
        @test provenance.species === gas
        @test provenance.scope === :no_frequency_redistribution
        @test provenance.cross_section_convention ===
              :total_molecular_rayleigh_with_effective_king_factor
        @test provenance.raman_scattering === :disabled
        @test provenance.wavelength_medium === :vacuum
        @test provenance.adopted_formula_wavenumber_range_cm1[1] <
              provenance.adopted_formula_wavenumber_range_cm1[2]
        @test provenance.primary_dispersion_measurement_wavelength_range_nm[1] <
              provenance.primary_dispersion_measurement_wavelength_range_nm[2]
        @test !isempty(provenance.coefficient_compilation_doi)
        @test !isempty(provenance.dispersion_reference_doi)
        @test first(provenance.direct_validation).wavelength_range_nm ==
              (333.0, 363.0)
        @test all(!isempty(v.doi) for v in provenance.direct_validation)
        @test all(!isempty(v.doi)
                  for v in provenance.precision_refractivity_validation)
        @test provenance.reference_temperature_k == 288.15
        @test provenance.reference_pressure_pa == 101325.0
    end

    ar = molecular_rayleigh_provenance(ri_method, :Ar)
    @test ar.primary_dispersion_measurement_wavelength_range_nm ==
          (467.9465, 2058.6909)
    @test ar.adopted_formula_wavenumber_range_cm1 == (5000.0, 33000.0)
    @test all(isapprox.(ar.adopted_formula_wavelength_range_nm,
                        (1.0e7 / 33000.0, 2000.0)))
    @test ar.king_factor_kind === :monatomic_identity
    @test isempty(ar.king_factor_measurement_wavelengths_nm)
    @test ar.king_factor_interpolation_wavelength_range_nm === nothing
    @test ar.king_factor_fit_source == ar.king_factor_source
    @test ar.king_factor_fit_reference_doi == ar.king_factor_reference_doi
    @test length(ar.direct_validation) == 2
    ar_thalman = ar.direct_validation[2]
    @test ar_thalman.coverage_kind === :selected_regions
    @test ar_thalman.wavelength_range_nm === nothing
    @test ar_thalman.reported_outer_bounds_nm == (350.0, 660.0)
    @test ar_thalman.tabulated_wavelengths_nm ==
          (370.0, 405.8, 532.2, 660.0)
    @test ar_thalman.reported_absolute_accuracy_bound_percent == 1.3
    @test ar_thalman.correction_doi ==
          "https://doi.org/10.1016/j.jqsrt.2016.12.014"

    co2 = molecular_rayleigh_provenance(ri_method, :CO2)
    @test co2.primary_dispersion_measurement_wavelength_range_nm ==
          (180.7, 253.7)
    @test co2.adopted_formula_wavelength_range_nm == (180.7, 1694.5)
    @test all(isapprox.(co2.adopted_formula_wavenumber_range_cm1,
                        (1.0e7 / 1694.5, 1.0e7 / 180.7)))
    @test co2.adopted_formula_range_source_commit ==
          "37773dce8128ccabcd9d5108d9958d6041244661"
    @test co2.adopted_formula_range_source_path ==
          "database/data/main/CO2/nk/Bideau-Mehu.yml"
    @test co2.adopted_formula_range_source_sha256 ==
          "d3b52d2a0b667a2e95aec110f35a7fced0c2ab3f02583d2e604aeb58c204f55c"
    @test co2.adopted_formula_range_reference_doi ==
          "https://doi.org/10.1038/s41597-023-02898-2"
    @test occursin(co2.adopted_formula_range_source_commit,
                   co2.adopted_formula_range_source_url)
    @test co2.king_factor_kind === :empirical_fit
    @test co2.king_factor_measurement_wavelengths_nm ==
          (457.9, 488.0, 514.5, 568.2, 647.1)
    @test co2.king_factor_interpolation_wavelength_range_nm ==
          (457.9, 647.1)
    @test co2.king_factor_fit_residual_bound_percent == 0.5
    @test co2.king_factor_measurement_reference_doi ==
          "https://doi.org/10.1063/1.431821"
    @test co2.king_factor_fit_reference_doi ==
          "https://doi.org/10.1016/j.jqsrt.2004.07.025"
    @test co2.direct_validation[2].coverage_kind === :continuous
    @test co2.direct_validation[2].wavelength_range_nm == (307.0, 725.0)

    @test_throws ArgumentError molecular_rayleigh_provenance(ri_method, :O2)
    @test_throws ArgumentError molecular_rayleigh_provenance(ri_method, :N2)
end
