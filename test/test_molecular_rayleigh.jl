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

@testset "molecular Rayleigh method types" begin
    @test RefractiveIndexRayleigh === UsingRefInd
    @test RefractiveIndexRayleigh() isa UsingRefInd
    @test RefractiveIndexRayleigh() === UsingRefInd()
    @test UsingRefInd() isa AbstractMolecularRayleighMethod
    @test UsingPolTen() isa AbstractMolecularRayleighMethod
    @test !(UsingPolTen() isa UsingRefInd)
    # The same names are reachable through the Scattering submodule.
    @test vSmartMOM.Scattering.UsingRefInd === UsingRefInd
    @test vSmartMOM.Scattering.UsingPolTen === UsingPolTen
    @test vSmartMOM.Scattering.RefractiveIndexRayleigh === UsingRefInd

    # Two-argument forms remain the refractive-index method for every
    # species ExoOptics calls them with.
    for gas in (:He, :Ar, :N2, :O2, :CO2), λ in (350.0, 532.0, 800.0)
        @test molecular_rayleigh_properties(gas, λ) ==
              molecular_rayleigh_properties(UsingRefInd(), gas, λ)
        @test molecular_rayleigh_cross_section_ratio(gas, λ, 400.0) ==
              molecular_rayleigh_cross_section_ratio(UsingRefInd(), gas, λ, 400.0)
    end
    @test molecular_rayleigh_provenance(UsingRefInd(), :CO2) ==
          molecular_rayleigh_provenance(RefractiveIndexRayleigh(), :CO2)
end

@testset "molecular Rayleigh polarizability tensor" begin
    pt = UsingPolTen()
    fields = (:species, :wavelength_nm, :refractive_index, :king_factor,
              :depolarization, :cross_section_cm2)
    @test keys(molecular_rayleigh_properties(pt, :N2, 400.0)) == fields
    @test keys(molecular_rayleigh_properties(UsingRefInd(), :N2, 400.0)) == fields

    for gas in (:N2, :O2, :H2)
        p = molecular_rayleigh_properties(pt, gas, 400.0)
        @test p.species === gas
        @test p.wavelength_nm == 400.0
        @test p.cross_section_cm2 > 0
        @test 0 < p.depolarization < 1
        # King factor and depolarization obey the same identity as the RI
        # method and get_greek_rayleigh.
        @test p.king_factor ≈ (6 + 3p.depolarization) / (6 - 7p.depolarization) rtol = 1e-12
        @test p.refractive_index > 1
        @test molecular_rayleigh_properties(pt, gas, 360.0).cross_section_cm2 >
              p.cross_section_cm2

        # Per-molecule quantities do not depend on the number density used
        # for the Lorentz–Lorenz refractive index.
        q = molecular_rayleigh_properties(pt, gas, 400.0;
                                          number_density_cm3 = nothing)
        @test q.refractive_index === nothing
        @test q.cross_section_cm2 == p.cross_section_cm2
        @test q.depolarization == p.depolarization
        dense = molecular_rayleigh_properties(pt, gas, 400.0;
                                              number_density_cm3 = 2 * 2.546899e19)
        @test dense.cross_section_cm2 == p.cross_section_cm2
        @test dense.refractive_index - 1 ≈ 2 * (p.refractive_index - 1) rtol = 1e-3

        # Lorentz–Lorenz round trip: the RI formula applied to the PolTen
        # refractive index at the same N reproduces the tensor cross section.
        N = 2.546899e19
        x = p.refractive_index - 1
        L = x * (2 + x) / (3 + x * (2 + x))
        ν = 1.0e7 / 400.0
        σ_LL = 24π^3 * (ν^2 / N)^2 * L^2 * p.king_factor
        @test σ_LL ≈ p.cross_section_cm2 rtol = 1e-10

        # Ratio API dispatches on the method and holds T fixed; N cancels.
        r = molecular_rayleigh_cross_section_ratio(pt, gas, 360.0, 400.0)
        @test r ≈ molecular_rayleigh_properties(pt, gas, 360.0).cross_section_cm2 /
                  p.cross_section_cm2
        @test molecular_rayleigh_cross_section_ratio(pt, gas, 360.0, 400.0;
                                                     temperature_K = 250.0) ≈
              molecular_rayleigh_properties(pt, gas, 360.0; temperature_K = 250.0).cross_section_cm2 /
              molecular_rayleigh_properties(pt, gas, 400.0; temperature_K = 250.0).cross_section_cm2

        # Float32: the tensor products are evaluated in Float64 internally
        # (α̅² ~ 1e-48 underflows Float32), and returned in Float32.
        p32 = molecular_rayleigh_properties(pt, gas, 400.0f0)
        @test p32.wavelength_nm isa Float32
        @test p32.cross_section_cm2 isa Float32
        @test p32.depolarization isa Float32
        @test p32.king_factor isa Float32
        @test p32.refractive_index isa Float32
        @test isfinite(p32.cross_section_cm2) && p32.cross_section_cm2 > 0
        @test p32.cross_section_cm2 ≈ Float32(p.cross_section_cm2) rtol = 1f-6
        @test p32.depolarization ≈ Float32(p.depolarization) rtol = 1f-6
        @test molecular_rayleigh_cross_section_ratio(pt, gas, 360.0f0, 400.0f0) isa Float32
    end

    # Aliases share the canonical species.
    @test molecular_rayleigh_properties(pt, :N₂, 500.0) ==
          molecular_rayleigh_properties(pt, "N2", 500.0)
    @test molecular_rayleigh_properties(pt, :H₂, 500.0).species === :H2
    @test molecular_rayleigh_properties(pt, :oxygen, 500.0).species === :O2

    # Exactly the Raman-module numbers that feed the RT model's default
    # (`depol: -1`) molecular depolarization path.
    IS = vSmartMOM.InelasticScattering
    for λ in (350.0, 532.0, 760.0), T in (220.0, 288.15)
        ν = 1.0e7 / λ
        n2, o2 = IS.getRamanAtmoConstants(ν, T)
        for (gas, mol) in ((:N2, n2), (:O2, o2))
            p = molecular_rayleigh_properties(pt, gas, λ; temperature_K = T)
            @test p.depolarization == mol.effCoeff.rho_depol_Rayl
            @test p.cross_section_cm2 ≈ mol.effCoeff.σ_Rayl_coeff * ν^4 rtol = 1e-14
        end
    end

    # Independent anchor: closed form evaluated here from the raw tensor
    # constants as written in src/Inelastic/src/molecular_constructors.jl
    # (α̅₀₀, ω₀, α_b, α_c, γ̅₀₀; cgs units). The first check catches an edit
    # to those constants; the second checks the Raman-module algebra and the
    # adapter against the textbook expressions
    #   α̅ = α̅₀₀(1 + α_b T + α_c T²)/(1 - (2πcν/ω₀)²),
    #   ρ = 6γ²/(45α̅² + 7γ²),  F_K = 1 + 2γ²/(9α̅²),
    #   σ = (128/3) π⁵ ν⁴ α̅² F_K.
    # The pinned 400 nm, 288.15 K values (5 significant figures) guard
    # against silent drift of both.
    raw = (N2 = (α̅₀₀ = 1.7406e-24, ω₀ = 2.6049e16, α_b = 1.8e-6,
                 α_c = 0.0, γ̅₀₀ = 0.71e-24, mol = IS.N₂()),
           O2 = (α̅₀₀ = 1.5658e-24, ω₀ = 2.1801e16, α_b = -2.369e-6,
                 α_c = 8.687e-9, γ̅₀₀ = 1.080e-24, mol = IS.O₂()),
           H2 = (α̅₀₀ = 0.8032e-24, ω₀ = 2.1399e16, α_b = 5.870e-6,
                 α_c = 7.544e-9, γ̅₀₀ = 0.288e-24, mol = IS.H₂()))
    pinned_400 = (N2 = (σ = 1.7103e-26, ρ = 0.020247),
                  O2 = (σ = 1.5082e-26, ρ = 0.054014),
                  H2 = (σ = 3.7447e-27, ρ = 0.015177))
    c_cgs = 2.99792458e10   # cm s⁻¹
    for gas in (:N2, :O2, :H2)
        r = getproperty(raw, gas)
        tensor = IS.getMolecularConstants(r.mol, 1.0).PolTensor
        @test (tensor.α̅₀₀, tensor.ω₀, tensor.α_b, tensor.α_c, tensor.γ̅₀₀) ==
              (r.α̅₀₀, r.ω₀, r.α_b, r.α_c, r.γ̅₀₀)
        for λ in (400.0, 650.0), T in (250.0, 288.15)
            ν = 1.0e7 / λ
            α̅ = r.α̅₀₀ * (1 + r.α_b * T + r.α_c * T^2) /
                 (1 - (2π * c_cgs * ν / r.ω₀)^2)
            γ = r.γ̅₀₀
            ρ = 6γ^2 / (45α̅^2 + 7γ^2)
            F = 1 + 2γ^2 / (9α̅^2)
            σ = 128 / 3 * π^5 * ν^4 * α̅^2 * F
            p = molecular_rayleigh_properties(pt, gas, λ; temperature_K = T)
            @test p.cross_section_cm2 ≈ σ rtol = 1e-12
            @test p.depolarization ≈ ρ rtol = 1e-12
            @test p.king_factor ≈ F rtol = 1e-12
        end
        p400 = molecular_rayleigh_properties(pt, gas, 400.0)
        @test p400.cross_section_cm2 ≈ getproperty(pinned_400, gas).σ rtol = 1e-4
        @test p400.depolarization ≈ getproperty(pinned_400, gas).ρ rtol = 1e-4
    end

    # The tensor anisotropy is λ- and T-independent; only α̅ carries the
    # Buldakov dispersion and temperature dependence.
    hook = vSmartMOM.Scattering._poltensor_rayleigh_coefficients
    for gas in (:N2, :O2, :H2)
        samples = [hook(gas, ν, T) for (ν, T) in ((1.0e7 / 350.0, 200.0),
                                                   (1.0e7 / 532.0, 288.15),
                                                   (1.0e7 / 800.0, 300.0),
                                                   (25_000.0, 1000.0))]
        @test all(c.γ̅ == first(samples).γ̅ for c in samples)
        @test first(samples).γ̅ == getproperty(raw, gas).γ̅₀₀
        @test length(unique(c.α̅ for c in samples)) == length(samples)
    end
    @test molecular_rayleigh_properties(pt, :N2, 500.0; temperature_K = 200.0).cross_section_cm2 !=
          molecular_rayleigh_properties(pt, :N2, 500.0; temperature_K = 300.0).cross_section_cm2

    # Unsupported species: clear ArgumentError naming the species and the
    # RI alternative.
    for gas in (:He, :Ar, :CO2, :H2O)
        err = try
            molecular_rayleigh_properties(pt, gas, 400.0)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(string(gas), err.msg)
        @test occursin("UsingRefInd()", err.msg)
    end
    @test_throws ArgumentError molecular_rayleigh_cross_section_ratio(pt, :Ar, 360.0, 400.0)
    @test_throws ArgumentError molecular_rayleigh_properties(UsingRefInd(), :H2, 400.0)
    @test_throws ArgumentError molecular_rayleigh_properties(pt, :N2, 0.0)
    @test_throws ArgumentError molecular_rayleigh_properties(pt, :N2, NaN)
    @test_throws ArgumentError molecular_rayleigh_properties(pt, :N2, 400.0; temperature_K = 0.0)
    @test_throws ArgumentError molecular_rayleigh_properties(pt, :N2, 400.0; number_density_cm3 = -1.0)
    # Provenance records are currently defined for the RI method only.
    @test_throws MethodError molecular_rayleigh_provenance(pt, :N2)
end

@testset "molecular Rayleigh PolTen vs RefInd intercomparison" begin
    # The two methods use independent data: RI fits (N2: Bates 1984 /
    # Peck & Khanna UV dispersion with the Bates King factor; O2: He et al.
    # 2021 Table 1, dispersion fit at 273.15 K and Bates-type King factor)
    # versus the Raman-module polarizability tensor (Buldakov-type α̅(ν,T),
    # wavelength-independent γ). At 288.15 K and 350-800 nm the observed
    # differences are (PolTen/RI - 1):
    #   N2: σ  -0.51 % (350 nm) ... +0.46 % (800 nm); n-1 within 0.14 %;
    #       F_K within 0.27 %; ρ -7.3 % (350 nm) ... +5.2 % (800 nm).
    #   O2: σ  -1.54 % (350 nm) ... +0.05 % (800 nm); n-1 within 0.46 %;
    #       F_K -1.60 % ... +0.98 %; ρ -15 % (350 nm) ... +11 % (800 nm).
    # The thresholds below are sanity/regression bounds informed by the
    # literature; they are NOT derived error bounds of either method. He et
    # al. (2021) state that the N2 Rayleigh cross section "is validated up
    # to 468 nm with an uncertainty of 1 %", and report average agreement of
    # measured O2 cross sections with n-based values of 1.2 % (pure O2) and
    # (2.8 ± 1.2) % (extinction), over a range ending at 725 nm. The test
    # extends to 800 nm, beyond both RI fit domains (N2 254-468 nm, O2
    # 288-546 nm per He et al. Table 1), so it compares two extrapolations
    # there. We therefore require σ to agree within 1 % (N2) and 2 % (O2),
    # i.e. within the size of the literature RI-vs-measurement differences.
    # The same thresholds are applied to F_K. The refractivity threshold is
    # half of that (σ ∝ (n-1)² F_K to first order), but the n-1 and F_K
    # contributions to σ can partly compensate, so these are independent
    # checks rather than a budget. The depolarization ratio is only
    # loosely constrained (20 %): the Raman-module anisotropy has no
    # dispersion, so ρ_PolTen rises with λ while the RI King factors fall,
    # and ρ ∝ F_K - 1 magnifies a small F_K difference into a
    # ~10 % relative ρ difference. That bound is a sanity check, not a
    # validation of either depolarization.
    tolerances = (N2 = (σ = 0.01, n = 0.005, F = 0.01),
                  O2 = (σ = 0.02, n = 0.010, F = 0.02))
    # Each RI fit carries its own reference density; compare n at that N.
    reference_density = (N2 = 2.546899e19, O2 = 2.68678e19)
    for gas in (:N2, :O2), λ in 350.0:50.0:800.0
        tol = getproperty(tolerances, gas)
        ri = molecular_rayleigh_properties(UsingRefInd(), gas, λ)
        pt = molecular_rayleigh_properties(UsingPolTen(), gas, λ;
                                           temperature_K = 288.15,
                                           number_density_cm3 =
                                               getproperty(reference_density, gas))
        @test abs(pt.cross_section_cm2 / ri.cross_section_cm2 - 1) < tol.σ
        @test abs(pt.king_factor / ri.king_factor - 1) < tol.F
        @test abs((pt.refractive_index - 1) / (ri.refractive_index - 1) - 1) < tol.n
        @test abs(pt.depolarization / ri.depolarization - 1) < 0.2
    end
end
