using Test
using vSmartMOM

# Published pure-gas Rayleigh dispersion and King-factor fits used by the
# ExoOptics Horak composition comparison. These tests pin the 360--400 nm
# interval without conflating absolute cross section with the matched-tau
# experiment (which consumes only sigma(lambda)/sigma(400)).
@testset "molecular Rayleigh species" begin
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
        p360 = molecular_rayleigh_properties(gas, 360.0)
        p400 = molecular_rayleigh_properties(gas, 400.0)
        @test p400.species == gas
        # The implementation uses cancellation- and overflow-safe algebra for
        # both Float64 and Float32. Algebraically equivalent evaluation orders
        # differ by a few parts in 10⁻¹² from the original pinned values.
        @test p400.cross_section_cm2 ≈ getproperty(expected_sigma_400, gas) rtol = 1e-11
        @test molecular_rayleigh_cross_section_ratio(gas, 360.0, 400.0) ≈
              getproperty(expected_ratio_360_400, gas) rtol = 1e-11
        @test p360.cross_section_cm2 > p400.cross_section_cm2 > 0
        @test 0 ≤ p400.depolarization < 1

        p400_f32 = molecular_rayleigh_properties(gas, 400.0f0)
        ratio_f32 = molecular_rayleigh_cross_section_ratio(gas, 360.0f0, 400.0f0)
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
end
