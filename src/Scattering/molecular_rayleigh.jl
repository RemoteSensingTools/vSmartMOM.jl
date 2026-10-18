"""
    AbstractMolecularRayleighMethod

Dispatch category for independently implemented molecular-Rayleigh
cross-section methods with no explicit frequency redistribution.  The method
is an explicit argument so a
refractive-index calculation can never be confused with a polarizability-
tensor calculation in provenance records or calling code.

This abstraction is deliberately part of `Scattering`, not `Inelastic`: the
HWO calculation must not alter or depend on the RRS molecular state.
"""
abstract type AbstractMolecularRayleighMethod end

"""
    RefractiveIndexRayleigh()

Select the conventional refractive-index plus effective-King-factor Rayleigh
formulation without constructing frequency-shifted Raman branches.
The HWO A--D experiment uses this explicit method for the newly required Ar
and CO₂ cross sections.  A future polarizability-tensor implementation must
use a different dispatch type; it must not silently replace this method.
"""
struct RefractiveIndexRayleigh <: AbstractMolecularRayleighMethod end

"""
    molecular_rayleigh_properties(::RefractiveIndexRayleigh,
                                   species, wavelength_nm)
    molecular_rayleigh_properties(species, wavelength_nm)

Return no-frequency-redistribution Rayleigh properties for a pure gas at
vacuum wavelength `wavelength_nm` (nm). Supported species are `:He`, `:Ar`,
`:N2`, `:O2`, and `:CO2` (Unicode aliases such as `:N₂` are accepted).

The HWO workflow uses the method-explicit form for Ar and CO₂. The two-
argument form and the He/N₂/O₂ coefficient sets predate that workflow and
remain compatibility behavior; they are not part of the new provenance review.

The returned named tuple contains

- `refractive_index`: refractive index at the source reference number density;
- `king_factor`: molecular-anisotropy correction ``F_K``;
- `depolarization`: depolarization ratio ``rho`` in the convention consumed by
  [`get_greek_rayleigh`](@ref), related by
  ``F_K = (6 + 3rho)/(6 - 7rho)``;
- `cross_section_cm2`: Rayleigh scattering cross section in cm² molecule⁻¹.

The individual dispersion and King-factor sources are returned by
[`molecular_rayleigh_provenance`](@ref) for the HWO additions Ar and CO₂.
Wilmouth and Sayres (2019) provide a direct 333--363 nm check of both gases.
The cross-section expression is

```math
sigma_R = \\frac{24 pi^3 nu^4}{N^2}
          \\left(\\frac{n^2-1}{n^2+2}\\right)^2 F_K,
```

where `nu` is in cm⁻¹ and `N` in molecule cm⁻³. The implementation can
evaluate any positive wavelength, but callers must use the distinct
dispersion, King-factor, and validation evidence returned by
`molecular_rayleigh_provenance` and mark extrapolation. The returned
conventional total molecular-Rayleigh cross section uses an effective King
factor; it is not a separately computed Cabannes-only cross section. No
species-specific rotational or vibrational Raman redistribution is built.

# References

- Wilmouth and Sayres (2019), *JQSRT* 234, 32--39,
  https://doi.org/10.1016/j.jqsrt.2019.05.031.
- Sneep and Ubachs (2005), *JQSRT* 92, 293--310,
  https://doi.org/10.1016/j.jqsrt.2004.07.025.
- He et al. (2021), *Atmospheric Chemistry and Physics* 21, 14927--14940,
  https://doi.org/10.5194/acp-21-14927-2021.
"""
function molecular_rayleigh_properties(::RefractiveIndexRayleigh,
                                        species, wavelength_nm::Real)
    gas = _canonical_rayleigh_species(species)
    λ = float(wavelength_nm)
    isfinite(λ) && λ > zero(λ) ||
        throw(ArgumentError("wavelength_nm must be finite and positive; got $wavelength_nm"))

    ν = oftype(λ, 1.0e7) / λ
    refractivity, number_density = _rayleigh_refractivity(gas, ν)
    refractivity_fraction = refractivity * oftype(λ, 1.0e-8)
    refractive_index = one(λ) + refractivity_fraction
    king_factor = _rayleigh_king_factor(gas, ν)
    king_factor >= one(king_factor) ||
        throw(DomainError(king_factor, "Rayleigh King factor must be at least one"))

    depolarization = oftype(λ, 6) * (king_factor - one(king_factor)) /
                     (oftype(λ, 7) * king_factor + oftype(λ, 3))
    # With n = 1 + x, use (n² - 1)/(n² + 2) = x(2 + x)/(3 + x(2 + x)).
    # This avoids cancellation in n² - 1 when n is stored as Float32.
    n²_minus_one = refractivity_fraction * (oftype(λ, 2) + refractivity_fraction)
    lorentz_lorenz = n²_minus_one / (oftype(λ, 3) + n²_minus_one)
    # Evaluate ν⁴/N² as (ν²/N)².  Squaring the standard number
    # density first overflows Float32 (~6.5e38), even though the final cross
    # section is well within its representable range.
    ν²_over_N = ν^2 / number_density
    cross_section_cm2 = oftype(λ, 24) * oftype(λ, π)^3 * ν²_over_N^2 *
                        lorentz_lorenz^2 * king_factor

    return (
        species = gas,
        wavelength_nm = λ,
        refractive_index = refractive_index,
        king_factor = king_factor,
        depolarization = depolarization,
        cross_section_cm2 = cross_section_cm2,
    )
end

molecular_rayleigh_properties(species, wavelength_nm::Real) =
    molecular_rayleigh_properties(RefractiveIndexRayleigh(), species,
                                  wavelength_nm)

"""
    molecular_rayleigh_cross_section_ratio(::RefractiveIndexRayleigh,
                                            species, wavelength_nm,
                                            reference_nm)
    molecular_rayleigh_cross_section_ratio(species, wavelength_nm, reference_nm)

Return ``sigma_R(wavelength_nm) / sigma_R(reference_nm)`` for a pure gas.
This ratio is the appropriate scaling when a comparison fixes the Rayleigh
optical depth at one reference wavelength and varies molecular identity.
"""
function molecular_rayleigh_cross_section_ratio(method::RefractiveIndexRayleigh,
                                                 species,
                                                 wavelength_nm::Real,
                                                 reference_nm::Real)
    σ = molecular_rayleigh_properties(method, species,
                                      wavelength_nm).cross_section_cm2
    σ_ref = molecular_rayleigh_properties(method, species,
                                          reference_nm).cross_section_cm2
    return σ / σ_ref
end

molecular_rayleigh_cross_section_ratio(species, wavelength_nm::Real,
                                       reference_nm::Real) =
    molecular_rayleigh_cross_section_ratio(RefractiveIndexRayleigh(), species,
                                           wavelength_nm, reference_nm)

"""
    molecular_rayleigh_provenance(::RefractiveIndexRayleigh, species)

Return the coefficient provenance and separate evidence intervals for the HWO
refractive-index/King-factor additions Ar and CO₂. The primary RI measurement
range, adopted dispersion-formula envelope, King-factor support, and direct
scattering checks are distinct facts. Callers spanning a wider interval must
record and approve each extrapolation; this function does not disguise an
adopted formula envelope as measured or directly validated coverage.

Ar uses the Peck--Fisher dispersion relation and unity King factor.  CO₂ uses
the Bideau-Mehu dispersion relation (scaled to the stated reference density)
and the Sneep--Ubachs King-factor fit to Alms--Burnham--Flygare measurements.
Wilmouth and Sayres (2019) directly measured both gases at 333--363 nm;
Thalman et al. (2014) checked Ar in selected regions between 350 and 660 nm;
and He et al. (2021) provide a continuous 307--725 nm direct check for CO₂.
This provenance API is intentionally limited to Ar and CO₂ because those are
the only newly required RI cross sections in the initial HWO A--D campaign.
"""
function molecular_rayleigh_provenance(::RefractiveIndexRayleigh, species)
    gas = _canonical_rayleigh_species(species)
    gas in (:Ar, :CO2) || throw(ArgumentError(
        "HWO RI provenance is available only for Ar and CO2; got $species"))

    if gas === :Ar
        # Sneep & Ubachs quote this domain in wavenumber.
        adopted_formula_range_cm1 = (5000.0, 33000.0)
        adopted_formula_range_nm =
            (1.0e7 / adopted_formula_range_cm1[2],
             1.0e7 / adopted_formula_range_cm1[1])
    else
        # The frozen refractiveindex.info record quotes this envelope in μm.
        adopted_formula_range_nm = (180.7, 1694.5)
        adopted_formula_range_cm1 =
            (1.0e7 / adopted_formula_range_nm[2],
             1.0e7 / adopted_formula_range_nm[1])
    end
    primary_measurement_range_nm = gas === :Ar ?
        (467.9465, 2058.6909) : (180.7, 253.7)
    dispersion_reference = gas === :Ar ?
        "Peck & Fisher (1964), Dispersion of Argon" :
        "Bideau-Mehu et al. (1973), Interferometric determination of the refractive index of carbon dioxide in the ultraviolet region"
    dispersion_reference_doi = gas === :Ar ?
        "https://doi.org/10.1364/JOSA.54.001362" :
        "https://doi.org/10.1016/0030-4018(73)90289-7"
    implemented_formula_reference = gas === :Ar ?
        dispersion_reference :
        "Sneep & Ubachs (2005), corrected Bideau-Mehu formula scaled to 15 C and 1013 hPa"
    implemented_formula_reference_doi = gas === :Ar ?
        dispersion_reference_doi :
        "https://doi.org/10.1016/j.jqsrt.2004.07.025"
    king_factor_kind = gas === :Ar ? :monatomic_identity : :empirical_fit
    king_factor_source = gas === :Ar ?
        "Sneep & Ubachs (2005), F_K=1 for spherical argon" :
        "Sneep & Ubachs (2005), Eq. (14), fit to Alms-Burnham-Flygare measurements"
    king_factor_reference_doi =
        "https://doi.org/10.1016/j.jqsrt.2004.07.025"
    # Keep the older generic source names, while exposing explicit fit names
    # so audit consumers cannot confuse the fit paper with the underlying
    # anisotropy measurements.
    king_factor_fit_source = king_factor_source
    king_factor_fit_reference_doi = king_factor_reference_doi
    king_factor_measurement_source = gas === :Ar ? nothing :
        "Alms, Burnham & Flygare (1975), Measurement of the dispersion in polarizability anisotropies"
    king_factor_measurement_reference_doi = gas === :Ar ? nothing :
        "https://doi.org/10.1063/1.431821"
    king_factor_measurement_wavelengths_nm = gas === :Ar ? () :
        (457.9, 488.0, 514.5, 568.2, 647.1)
    king_factor_interpolation_wavelength_range_nm = gas === :Ar ? nothing :
        (457.9, 647.1)
    king_factor_fit_residual_bound_percent = gas === :Ar ? nothing : 0.5
    king_factor_domain_note = gas === :Ar ?
        "F_K=1 is a monatomic-gas identity, not an empirical fit with a finite measurement interval" :
        "the empirical interpolation span is 457.9-647.1 nm; use outside it is King-factor extrapolation unless separately qualified"
    king_factor_coefficients = gas === :Ar ?
        (constant = 1.0, wavenumber_squared_cm2 = 0.0) :
        (constant = 1.1364, wavenumber_squared_cm2 = 2.53e-11)
    king_factor_coefficient_uncertainties = gas === :Ar ? nothing :
        (constant = 0.0005, wavenumber_squared_cm2 = 1.5e-12)
    wilmouth_validation = (
        source = "Wilmouth & Sayres (2019)",
        doi = "https://doi.org/10.1016/j.jqsrt.2019.05.031",
        coverage_kind = :continuous,
        wavelength_range_nm = (333.0, 363.0),
        reported_outer_bounds_nm = (333.0, 363.0),
        tabulated_wavelengths_nm = (),
        reported_mean_difference_percent = gas === :Ar ? 0.2 : 1.0,
        reported_difference_standard_deviation_percent = nothing,
        reported_mean_one_sigma_uncertainty_percent =
            gas === :Ar ? 2.5 : 1.9,
        reported_absolute_accuracy_bound_percent = nothing,
        correction_doi = nothing,
    )
    direct_validation = gas === :Ar ?
        (wilmouth_validation,
         (source = "Thalman et al. (2014)",
          doi = "https://doi.org/10.1016/j.jqsrt.2014.05.030",
          coverage_kind = :selected_regions,
          wavelength_range_nm = nothing,
          reported_outer_bounds_nm = (350.0, 660.0),
          tabulated_wavelengths_nm = (370.0, 405.8, 532.2, 660.0),
          reported_mean_difference_percent = 0.2,
          reported_difference_standard_deviation_percent = 0.4,
          reported_mean_one_sigma_uncertainty_percent = nothing,
          reported_absolute_accuracy_bound_percent = 1.3,
          correction_doi =
              "https://doi.org/10.1016/j.jqsrt.2016.12.014")) :
        (wilmouth_validation,
         (source = "He et al. (2021)",
          doi = "https://doi.org/10.5194/acp-21-14927-2021",
          coverage_kind = :continuous,
          wavelength_range_nm = (307.0, 725.0),
          reported_outer_bounds_nm = (307.0, 725.0),
          tabulated_wavelengths_nm = (),
          reported_mean_difference_percent = 0.4,
          reported_difference_standard_deviation_percent = 1.2,
          reported_mean_one_sigma_uncertainty_percent = (1.2, 1.3),
          reported_absolute_accuracy_bound_percent = nothing,
          correction_doi = nothing))
    precision_refractivity_validation = gas === :Ar ?
        ((doi = "https://doi.org/10.1007/s10765-023-03291-2",
          wavelengths_nm = (632.9908,)),
         (doi = "https://doi.org/10.1007/s10765-024-03412-5",
          wavelengths_nm = (1542.383,))) :
        ((doi = "https://doi.org/10.1088/1361-6501/ae6e94",
          wavelengths_nm = (632.9908, 1542.383)),)

    return (
        method = :refractive_index_king_factor,
        species = gas,
        coefficient_compilation =
            "Wilmouth & Sayres (2019), JQSRT 234, 32-39, Table 1",
        coefficient_compilation_doi =
            "https://doi.org/10.1016/j.jqsrt.2019.05.031",
        corrected_compilation_references =
            ("https://doi.org/10.1016/j.jqsrt.2004.07.025",
             "https://doi.org/10.1016/j.jqsrt.2019.05.031"),
        dispersion_reference = dispersion_reference,
        dispersion_reference_doi = dispersion_reference_doi,
        implemented_formula_reference = implemented_formula_reference,
        implemented_formula_reference_doi =
            implemented_formula_reference_doi,
        primary_dispersion_measurement_wavelength_range_nm =
            primary_measurement_range_nm,
        adopted_formula_wavenumber_range_cm1 = adopted_formula_range_cm1,
        adopted_formula_wavelength_range_nm = adopted_formula_range_nm,
        adopted_formula_range_basis = gas === :Ar ?
            "Sneep & Ubachs (2005) adopted Peck-Fisher interval" :
            "frozen refractiveindex.info Bideau-Mehu formula record; wider than the original UV measurements",
        adopted_formula_range_reference = gas === :Ar ?
            "https://doi.org/10.1016/j.jqsrt.2004.07.025" :
            "https://doi.org/10.1038/s41597-023-02898-2",
        adopted_formula_range_reference_doi = gas === :Ar ?
            "https://doi.org/10.1016/j.jqsrt.2004.07.025" :
            "https://doi.org/10.1038/s41597-023-02898-2",
        adopted_formula_range_source_repository = gas === :Ar ? nothing :
            "https://github.com/polyanskiy/refractiveindex.info-database",
        adopted_formula_range_source_commit = gas === :Ar ? nothing :
            "37773dce8128ccabcd9d5108d9958d6041244661",
        adopted_formula_range_source_path = gas === :Ar ? nothing :
            "database/data/main/CO2/nk/Bideau-Mehu.yml",
        adopted_formula_range_source_url = gas === :Ar ? nothing :
            "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/37773dce8128ccabcd9d5108d9958d6041244661/database/data/main/CO2/nk/Bideau-Mehu.yml",
        adopted_formula_range_source_sha256 = gas === :Ar ? nothing :
            "d3b52d2a0b667a2e95aec110f35a7fced0c2ab3f02583d2e604aeb58c204f55c",
        primary_measurement_reference_temperature_k =
            gas === :Ar ? 288.15 : 273.15,
        primary_measurement_reference_pressure_pa = 101325.0,
        king_factor_kind = king_factor_kind,
        king_factor_source = king_factor_source,
        king_factor_reference_doi = king_factor_reference_doi,
        king_factor_fit_source = king_factor_fit_source,
        king_factor_fit_reference_doi = king_factor_fit_reference_doi,
        king_factor_measurement_source = king_factor_measurement_source,
        king_factor_measurement_reference_doi =
            king_factor_measurement_reference_doi,
        king_factor_measurement_wavelengths_nm =
            king_factor_measurement_wavelengths_nm,
        king_factor_interpolation_wavelength_range_nm =
            king_factor_interpolation_wavelength_range_nm,
        king_factor_fit_residual_bound_percent =
            king_factor_fit_residual_bound_percent,
        king_factor_domain_note = king_factor_domain_note,
        king_factor_coefficients = king_factor_coefficients,
        king_factor_coefficient_uncertainties =
            king_factor_coefficient_uncertainties,
        direct_validation = direct_validation,
        precision_refractivity_validation = precision_refractivity_validation,
        reference_number_density_cm3 = 2.546899e19,
        reference_temperature_k = 288.15,
        reference_pressure_pa = 101325.0,
        wavelength_medium = :vacuum,
        cross_section_convention =
            :total_molecular_rayleigh_with_effective_king_factor,
        scope = :no_frequency_redistribution,
        raman_scattering = :disabled,
    )
end

function _canonical_rayleigh_species(species)
    key = species isa Symbol ? species : Symbol(species)
    key in (:He, :HE, :helium, :Helium) && return :He
    key in (:Ar, :AR, :argon, :Argon) && return :Ar
    key in (:N2, :N₂, :nitrogen, :Nitrogen) && return :N2
    key in (:O2, :O₂, :oxygen, :Oxygen) && return :O2
    key in (:CO2, :CO₂, :carbon_dioxide, :CarbonDioxide) && return :CO2
    throw(ArgumentError(
        "unsupported Rayleigh species $species; choose He, Ar, N2, O2, or CO2"))
end

# Refractivities are tabulated as (n - 1) * 1e8. The standard number density
# is 2.546899e19 molecule cm^-3 (288.15 K, 1013.25 hPa), except for the O2 fit
# marked with an asterisk by He et al. (2021), which uses 2.68678e19 molecule
# cm^-3 (273.15 K, 1013.25 hPa). Keeping the matching N with each fit is
# essential for an absolute cross section; it cancels in same-species ratios.
function _rayleigh_refractivity(gas::Symbol, ν)
    ν² = ν^2
    if gas === :He
        r = oftype(ν, 2283) + oftype(ν, 1.8102e13) /
            (oftype(ν, 1.5342e10) - ν²)
        N = oftype(ν, 2.546899e19)
    elseif gas === :Ar
        r = oftype(ν, 6432.135) + oftype(ν, 286.06021e12) /
            (oftype(ν, 1.44e10) - ν²)
        N = oftype(ν, 2.546899e19)
    elseif gas === :N2
        r = oftype(ν, 5677.465) + oftype(ν, 318.81874e12) /
            (oftype(ν, 1.44e10) - ν²)
        N = oftype(ν, 2.546899e19)
    elseif gas === :O2
        r = oftype(ν, 20564.8) + oftype(ν, 2.480899e13) /
            (oftype(ν, 4.09e9) - ν²)
        N = oftype(ν, 2.68678e19)
    elseif gas === :CO2
        r = oftype(ν, 1.1427e11) * (
            oftype(ν, 5799.25) / (oftype(ν, 128908.9)^2 - ν²) +
            oftype(ν, 120.05) / (oftype(ν, 89223.8)^2 - ν²) +
            oftype(ν, 5.3334) / (oftype(ν, 75037.5)^2 - ν²) +
            oftype(ν, 4.3244) / (oftype(ν, 67837.7)^2 - ν²) +
            oftype(ν, 1.218145e-5) /
                (oftype(ν, 2418.136)^2 - ν²))
        N = oftype(ν, 2.546899e19)
    else
        error("unreachable Rayleigh species $gas")
    end
    return r, N
end

function _rayleigh_king_factor(gas::Symbol, ν)
    ν² = ν^2
    if gas === :He || gas === :Ar
        return one(ν)
    elseif gas === :N2
        return oftype(ν, 1.034) + oftype(ν, 3.17e-12) * ν²
    elseif gas === :O2
        # Preserve the pre-HWO helper exactly. The A--D addition uses this
        # explicit RI method only for Ar and CO2; existing N2/O2 production
        # scattering paths and their numerical values are intentionally
        # outside this branch's scope.
        return oftype(ν, 1.09) + oftype(ν, 1.385e-11) * ν² +
               oftype(ν, 1.448e-20) * ν²^2
    elseif gas === :CO2
        return oftype(ν, 1.1364) + oftype(ν, 2.53e-11) * ν²
    end
    error("unreachable Rayleigh species $gas")
end
