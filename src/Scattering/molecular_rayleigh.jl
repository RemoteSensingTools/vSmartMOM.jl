"""
    AbstractMolecularRayleighMethod

Dispatch category for independently implemented molecular-Rayleigh
cross-section methods with no explicit frequency redistribution.  The method
is an explicit argument so a
refractive-index calculation can never be confused with a polarizability-
tensor calculation in provenance records or calling code.

Concrete methods:

- [`UsingRefInd`](@ref): refractive index plus effective King factor
  (`RefractiveIndexRayleigh` is a compatibility alias);
- [`UsingPolTen`](@ref): mean polarizability and anisotropy of the molecular
  polarizability tensor, with the wavelength dependence selected by an
  [`AbstractPolarizabilityDispersion`](@ref) (default
  [`DispersiveAnisotropy`](@ref); [`ConstantAnisotropy`](@ref) is the
  Raman-module tensor with a wavelength-independent anisotropy).

This abstraction is part of `Scattering`, not `Inelastic`. Only
`UsingRefInd` is independent of the Raman-module (`InelasticScattering`)
constants. `UsingPolTen` deliberately uses them: it builds fresh Raman-module
molecular constants on each call and never mutates shared Raman state.
"""
abstract type AbstractMolecularRayleighMethod end

"""
    AbstractPolarizabilityDispersion

Wavelength dependence of the N₂ and O₂ polarizability-tensor invariants, the
mean polarizability ``\\bar\\alpha`` and the anisotropy ``\\gamma``. Two
formulations are kept and selected by dispatch:

- [`DispersiveAnisotropy`](@ref) (default for [`UsingPolTen`](@ref) and for
  the automatic `depol: -1` Rayleigh depolarization of the RT model):
  ``\\bar\\alpha(\\nu)`` from the [`UsingRefInd`](@ref) refractivity through
  the Lorentz--Lorenz relation and ``\\gamma(\\nu)`` from its King factor;
- [`ConstantAnisotropy`](@ref) (default for the Raman/Cabannes cross
  sections): the published Raman-module tensor, Buldakov-type
  ``\\bar\\alpha(\\nu, T)`` and a wavelength-independent ``\\gamma``.

See `docs/dev_notes/molecular_rayleigh_dispersion.md` for the sources and the
size of the difference.
"""
abstract type AbstractPolarizabilityDispersion end

"""
    ConstantAnisotropy()

The polarizability tensor of the Raman module (`getMolecularConstants` in
`src/Inelastic/src/molecular_constructors.jl`): mean polarizability with the
Buldakov et al. single-pole dispersion and temperature dependence,
``\\bar\\alpha(\\nu,T) = \\bar\\alpha_{00}(1+\\alpha_b T+\\alpha_c T^2)
/(1-(2\\pi c\\nu/\\omega_0)^2)``, and the anisotropy ``\\gamma = \\bar\\gamma_{00}``,
independent of wavelength and temperature (0.71×10⁻²⁴ cm³ for N₂,
1.080×10⁻²⁴ cm³ for O₂). These constants produced the published Raman tables
(Sanghavi 2022, JQSRT 291, 108328; the Sanghavi 2026 Table 1 checked in
`test/test_ray_cab_raman.jl`), and they remain the default of the Raman and
Cabannes cross sections. For N₂ and O₂ the constant anisotropy equals the
King-factor anisotropy of [`DispersiveAnisotropy`](@ref) near 488 nm
(within 0.03 % for N₂ and 0.32 % for O₂); it is 2 % (N₂) and 4 % (O₂) low at 400 nm and 4 %
(N₂) and 8 % (O₂) high at 2060 nm.
"""
struct ConstantAnisotropy <: AbstractPolarizabilityDispersion end

"""
    DispersiveAnisotropy()

Polarizability-tensor invariants of N₂ and O₂ that follow the refractive-index
data of [`UsingRefInd`](@ref):

```math
\\bar\\alpha(\\nu, T) = \\frac{3}{4\\pi N}\\,\\frac{n^2-1}{n^2+2}\\,
                        \\frac{f(T)}{f(T_\\mathrm{fit})},\\qquad
\\gamma(\\nu) = \\frac{3}{4\\pi N}\\,\\frac{n^2-1}{n^2+2}\\,
                \\sqrt{\\tfrac{9}{2}\\left(F_K(\\nu)-1\\right)},
```

where ``n`` and ``N`` are the refractive index and number density of the
dispersion fit used at ``\\nu`` (reference temperature ``T_\\mathrm{fit}``),
``F_K`` is the `UsingRefInd` King factor (Bates 1984 for N₂; He et al. 2021,
Table 1, for O₂), and ``f(T) = 1 + \\alpha_b T + \\alpha_c T^2`` is the
Buldakov et al. temperature factor of the Raman-module mean polarizability.
The anisotropy is that implied by the measured depolarization ratios behind the
King-factor fits (Bates 1984, as in Bodhaine et al. 1999 and Tomasi et al.
2005), so the depolarization ``\\rho = 6\\gamma^2/(45\\bar\\alpha^2+7\\gamma^2)``
equals the `UsingRefInd` value ``6(F_K-1)/(7F_K+3)`` at ``T = T_\\mathrm{fit}``.
Validity follows the `UsingRefInd` fits: N₂ 254--2059 nm, O₂ 288--1542 nm
measured (longer wavelengths are dispersion extrapolations, which are smooth
because ``\\nu^2 \\to 0``); the King factors are fitted to visible
depolarization data and extrapolated into the UV and NIR.

H₂ has no refractive-index fit in `UsingRefInd`; for H₂ this formulation uses
the Raman-module tensor of [`ConstantAnisotropy`](@ref) unchanged.
"""
struct DispersiveAnisotropy <: AbstractPolarizabilityDispersion end

"""
    UsingRefInd()
    RefractiveIndexRayleigh()

Select the conventional refractive-index plus effective-King-factor Rayleigh
formulation without constructing frequency-shifted Raman branches.
The Habitable Worlds Observatory (HWO) exo-Earth studies use this explicit
method for the Ar and CO₂ cross sections. `RefractiveIndexRayleigh` is a compatibility alias
(`RefractiveIndexRayleigh === UsingRefInd`) kept so existing calling code,
e.g. ExoOptics, continues to work unchanged. The polarizability-tensor
implementation uses the distinct dispatch type [`UsingPolTen`](@ref); it never
silently replaces this method.
"""
struct UsingRefInd <: AbstractMolecularRayleighMethod end

# Compatibility alias: `RefractiveIndexRayleigh()` constructs `UsingRefInd()`.
const RefractiveIndexRayleigh = UsingRefInd

"""
    molecular_rayleigh_properties(::UsingRefInd, species, wavelength_nm)
    molecular_rayleigh_properties(species, wavelength_nm)

Return no-frequency-redistribution Rayleigh properties for a pure gas at
vacuum wavelength `wavelength_nm` (nm). Supported species are `:He`, `:Ar`,
`:N2`, `:O2`, and `:CO2` (Unicode aliases such as `:N₂` are accepted).

The two-argument form is equivalent to passing `UsingRefInd()`. The
machine-readable provenance record ([`molecular_rayleigh_provenance`](@ref))
currently covers Ar and CO₂ only; the He, N₂ and O₂ coefficient sets have no
such record yet (their sources are listed below).

The returned named tuple contains

- `refractive_index`: refractive index at the species reference number
  density: 2.546899×10¹⁹ cm⁻³ (288.15 K, 1013.25 hPa) for He, Ar, N₂ and CO₂
  and 2.68678×10¹⁹ cm⁻³ (273.15 K, 1013.25 hPa) for O₂. A fit given at
  another density is converted through the Lorentz--Lorenz relation
  (`(n²-1)/(n²+2)` proportional to `N`). It is stored as ``n`` itself, so a
  `Float32` result resolves ``n-1`` only to about `eps(1f0)/(n-1)` relative
  (`eps(1f0)` ≈ 1.2e-7): a few times 1e-4 for N₂, O₂, Ar and CO₂, and a few
  times 1e-3 for He. The cross section is computed from the refractivity
  directly and is not affected;
- `king_factor`: molecular-anisotropy correction ``F_K``;
- `depolarization`: depolarization ratio ``\\rho`` in the convention consumed by
  [`get_greek_rayleigh`](@ref), related by
  ``F_K = (6 + 3\\rho)/(6 - 7\\rho)``;
- `cross_section_cm2`: Rayleigh scattering cross section in cm² molecule⁻¹.

The individual dispersion and King-factor sources are returned by
[`molecular_rayleigh_provenance`](@ref) for Ar and CO₂.
Wilmouth and Sayres (2019) provide a direct 333--363 nm check of both gases.
The cross-section expression is

```math
\\sigma_R = \\frac{24 \\pi^3 \\nu^4}{N^2}
          \\left(\\frac{n^2-1}{n^2+2}\\right)^2 F_K,
```

where ``\\nu`` is in cm⁻¹ and `N` in molecule cm⁻³, the number density at
which the dispersion fit gives ``n``. Every fit is used with the ideal-gas
number density of its stated reference state, as in Bates (1984), Bodhaine et
al. (1999), Sneep and Ubachs (2005) and He et al. (2021); gas
compressibility (5--10×10⁻⁴ in ``n-1`` per molecule at 1 atm, 0--20 °C) is
neglected. The implementation can evaluate any positive wavelength, but
callers must use the distinct dispersion, King-factor, and validation evidence
returned by `molecular_rayleigh_provenance` and mark extrapolation. The
returned conventional total molecular-Rayleigh cross section uses an effective
King factor; it is not a separately computed Cabannes-only cross section. No
species-specific rotational or vibrational Raman redistribution is built.

# N₂ and O₂ dispersion and King factors

Each gas joins two dispersion fits at the wavenumber where they agree, so
``n-1`` and ``\\sigma`` are continuous (there is a kink in their slope):

- N₂, ``\\nu \\ge`` 22 919.551 cm⁻¹ (``\\lambda \\le`` 436.309 nm): Bates
  (1984), ``(n-1)\\times10^8 = 5677.465 + 318.81874\\times10^{12}/
  (1.44\\times10^{10}-\\nu^2)`` at 288.15 K (He et al. 2021, Table 1; fitted
  over 254--468 nm). Longer wavelengths: Peck and Khanna (1966),
  ``(n-1)\\times10^8 = 6855.200 + 324.3157\\times10^{12}/(1.44\\times10^{10}-\\nu^2)``
  at 0 °C and 760 Torr (measured 467.9--2058.7 nm; checked against the
  refractiveindex.info `N2/nk/Peck-0C.yml` record). The two fits differ by
  1.4×10⁻⁴ in ``n-1`` at 468 nm; joining them there would leave a 0.03 % jump
  in ``\\sigma``.
- O₂, ``\\nu \\ge`` 21 428.071 cm⁻¹ (``\\lambda \\le`` 466.678 nm): Bates
  (1984), ``(n-1)\\times10^8 = 20564.8 + 2.480899\\times10^{13}/
  (4.09\\times10^{9}-\\nu^2)`` at 273.15 K (He et al. 2021, Table 1; fitted
  over 288--546 nm). Longer wavelengths: Egan (2026), Eq. (8),
  ``(n-1)\\times10^8 = 5627886.8/(281.0-s^2) + 197258.6/(45.4-s^2)``,
  ``s = 1/\\lambda`` in µm⁻¹, at 20 °C and 100 kPa, constrained by refractivity
  at 0.633 and 1.542 µm (uncertainty 7×10⁻⁶ in ``n-1``), visible data and a
  0.254/0.436 µm ratio. At 546 nm the two fits differ by 0.10 % in ``n-1``
  (0.21 % in ``\\sigma``); at 1.6 µm the extrapolated Bates fit is 0.43 % high
  in ``n-1``.
- King factors: N₂ ``F_K = 1.034 + 3.17\\times10^{-12}\\nu^2`` (Bates 1984);
  O₂ ``F_K = 1.09 + 1.385\\times10^{-11}\\nu^2 + 1.448\\times10^{-20}\\nu^4``
  (He et al. 2021, Table 1). The O₂ constant 1.09 is He et al.'s; Bates (1984)
  and Bodhaine et al. (1999) use 1.096, which gives a 0.54 % larger O₂ cross
  section (0.54 % at 400 nm, 0.55 % at 1.6 µm). Both King factors are fitted to
  visible depolarization data and extrapolated elsewhere.

# References

- Bates (1984), *Planet. Space Sci.* 32, 785--790,
  https://doi.org/10.1016/0032-0633(84)90102-8.
- Bodhaine et al. (1999), *J. Atmos. Oceanic Technol.* 16, 1854--1861,
  https://doi.org/10.1175/1520-0426(1999)016<1854:ORODC>2.0.CO;2.
- Egan (2026), *Meas. Sci. Technol.* 37, 255003 (CC BY 4.0),
  https://doi.org/10.1088/1361-6501/ae6e94.
- He et al. (2021), *Atmospheric Chemistry and Physics* 21, 14927--14940,
  https://doi.org/10.5194/acp-21-14927-2021.
- Peck and Khanna (1966), *J. Opt. Soc. Am.* 56, 1059--1063,
  https://doi.org/10.1364/JOSA.56.001059.
- Sneep and Ubachs (2005), *JQSRT* 92, 293--310,
  https://doi.org/10.1016/j.jqsrt.2004.07.025.
- Wilmouth and Sayres (2019), *JQSRT* 234, 32--39,
  https://doi.org/10.1016/j.jqsrt.2019.05.031.
"""
function molecular_rayleigh_properties(::UsingRefInd,
                                        species, wavelength_nm::Real)
    gas = _canonical_rayleigh_species(species)
    λ = float(wavelength_nm)
    isfinite(λ) && λ > zero(λ) ||
        throw(ArgumentError("wavelength_nm must be finite and positive; got $wavelength_nm"))

    ν = oftype(λ, 1.0e7) / λ
    fit = _rayleigh_refractivity(gas, ν)
    number_density = fit.number_density
    refractivity_fraction = fit.refractivity * oftype(λ, 1.0e-8)
    king_factor = _rayleigh_king_factor(gas, ν)
    king_factor >= one(king_factor) ||
        throw(DomainError(king_factor, "Rayleigh King factor must be at least one"))

    depolarization = oftype(λ, 6) * (king_factor - one(king_factor)) /
                     (oftype(λ, 7) * king_factor + oftype(λ, 3))
    # With n = 1 + x, use (n² - 1)/(n² + 2) = x(2 + x)/(3 + x(2 + x)).
    # This avoids cancellation in n² - 1 when n is stored as Float32.
    n²_minus_one = refractivity_fraction * (oftype(λ, 2) + refractivity_fraction)
    lorentz_lorenz = n²_minus_one / (oftype(λ, 3) + n²_minus_one)
    # Report n at the species reference density; a fit given at another
    # density is rescaled through Lorentz–Lorenz, L ∝ N.
    refractive_index = number_density == fit.report_density ?
        one(λ) + refractivity_fraction :
        _lorentz_lorenz_refractive_index(
            lorentz_lorenz * (fit.report_density / number_density))
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

# Invert the Lorentz–Lorenz relation L = (n² - 1)/(n² + 2) without forming
# n² - 1 by subtraction: n² - 1 = 3L/(1 - L), n - 1 = (n² - 1)/(n + 1).
function _lorentz_lorenz_refractive_index(L)
    n²_minus_one = 3L / (one(L) - L)
    return one(L) + n²_minus_one / (sqrt(one(L) + n²_minus_one) + one(L))
end

molecular_rayleigh_properties(species, wavelength_nm::Real) =
    molecular_rayleigh_properties(UsingRefInd(), species, wavelength_nm)

"""
    molecular_rayleigh_cross_section_ratio(method::AbstractMolecularRayleighMethod,
                                            species, wavelength_nm,
                                            reference_nm)
    molecular_rayleigh_cross_section_ratio(::UsingPolTen, species,
                                            wavelength_nm, reference_nm;
                                            temperature_K = 288.15)
    molecular_rayleigh_cross_section_ratio(species, wavelength_nm, reference_nm)

Return ``\\sigma_R(\\lambda) / \\sigma_R(\\lambda_\\mathrm{ref})`` for a pure gas,
with ``\\lambda`` = `wavelength_nm` and ``\\lambda_\\mathrm{ref}`` = `reference_nm`.
This ratio is the appropriate scaling when a comparison fixes the Rayleigh
optical depth at one reference wavelength and varies molecular identity.
For [`UsingPolTen`](@ref) both cross sections are evaluated at the same
`temperature_K`; the number density cancels. The three-argument form uses
[`UsingRefInd`](@ref).
"""
function molecular_rayleigh_cross_section_ratio(method::AbstractMolecularRayleighMethod,
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
    molecular_rayleigh_cross_section_ratio(UsingRefInd(), species,
                                           wavelength_nm, reference_nm)

"""
    molecular_rayleigh_provenance(::UsingRefInd, species)

Return the coefficient provenance and separate evidence intervals for the
refractive-index/King-factor coefficients of Ar and CO₂. The primary RI
measurement range, adopted dispersion-formula envelope, King-factor support,
and direct scattering checks are reported as distinct fields, so callers can
tell measured or directly validated coverage from formula extrapolation.

Ar uses the Peck--Fisher dispersion relation and unity King factor.  CO₂ uses
the Bideau-Mehu dispersion relation (scaled to the stated reference density)
and the Sneep--Ubachs King-factor fit to Alms--Burnham--Flygare measurements.
Wilmouth and Sayres (2019) directly measured both gases at 333--363 nm;
Thalman et al. (2014) checked Ar in selected regions between 350 and 660 nm;
and He et al. (2021) provide a continuous 307--725 nm direct check for CO₂.
This provenance API is currently limited to Ar and CO₂.
"""
function molecular_rayleigh_provenance(::UsingRefInd, species)
    gas = _canonical_rayleigh_species(species)
    gas in (:Ar, :CO2) || throw(ArgumentError(
        "RI provenance is available only for Ar and CO2; got $species"))

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

# Refractivities are tabulated as (n - 1) * 1e8 at the number density N of
# each fit's reference state (ideal gas): 2.546899e19 molecule cm^-3 at
# 288.15 K and 1013.25 hPa; 2.68678e19 at 273.15 K and 1013.25 hPa (the O2 fit
# marked with an asterisk by He et al. 2021, Table 1, and the Peck & Khanna N2
# fit); and 100 kPa/(k_B 293.15 K) for the Egan (2026) O2 fit. Keeping the
# matching N with each fit is essential for an absolute cross section; it
# cancels in same-fit ratios. `report_density` is the density at which
# `molecular_rayleigh_properties` reports the refractive index of the species,
# and `temperature_K` the reference temperature of the fit (used by
# `DispersiveAnisotropy` for the temperature scaling of the mean
# polarizability).
#
# N2 and O2 each switch between two fits at the wavenumber where the per-
# molecule Lorentz–Lorenz volume [(n²-1)/(n²+2)]/N of the two fits is equal
# (found by bisection; the rounded joins below leave mismatches of 6e-12 and
# 8e-11). See the `molecular_rayleigh_properties(::UsingRefInd, ...)`
# docstring and docs/dev_notes/molecular_rayleigh_dispersion.md.
const _N2_RI_JOIN_CM1 = 22919.551   # 436.309 nm: Bates (1984) above, Peck & Khanna (1966) below
const _O2_RI_JOIN_CM1 = 21428.071   # 466.678 nm: Bates (1984) above, Egan (2026) below
const _RI_N_288K_CM3 = 2.546899e19
const _RI_N_273K_CM3 = 2.68678e19
# Egan (2026) quotes n at 20 °C and 100 kPa; ideal-gas density with exact k_B.
const _RI_N_EGAN_O2_CM3 = 1.0e5 / (1.380649e-23 * 293.15) * 1.0e-6

function _rayleigh_refractivity(gas::Symbol, ν)
    ν² = ν^2
    T = oftype(ν, 288.15)
    if gas === :He
        r = oftype(ν, 2283) + oftype(ν, 1.8102e13) /
            (oftype(ν, 1.5342e10) - ν²)
        N = oftype(ν, _RI_N_288K_CM3)
        N_report = N
    elseif gas === :Ar
        r = oftype(ν, 6432.135) + oftype(ν, 286.06021e12) /
            (oftype(ν, 1.44e10) - ν²)
        N = oftype(ν, _RI_N_288K_CM3)
        N_report = N
    elseif gas === :N2
        N_report = oftype(ν, _RI_N_288K_CM3)
        if ν >= oftype(ν, _N2_RI_JOIN_CM1)
            # Bates (1984), 254-468 nm, 15 °C (He et al. 2021, Table 1).
            r = oftype(ν, 5677.465) + oftype(ν, 318.81874e12) /
                (oftype(ν, 1.44e10) - ν²)
            N = N_report
        else
            # Peck & Khanna (1966), 467.9-2058.7 nm, 0 °C and 760 Torr:
            # 1e8(n - 1) = 6855.200 + 3243157.0/(144 - s²), s in µm⁻¹.
            r = oftype(ν, 6855.2) + oftype(ν, 324.3157e12) /
                (oftype(ν, 1.44e10) - ν²)
            N = oftype(ν, _RI_N_273K_CM3)
            T = oftype(ν, 273.15)
        end
    elseif gas === :O2
        N_report = oftype(ν, _RI_N_273K_CM3)
        if ν >= oftype(ν, _O2_RI_JOIN_CM1)
            # Bates (1984), 288-546 nm, 0 °C (He et al. 2021, Table 1).
            r = oftype(ν, 20564.8) + oftype(ν, 2.480899e13) /
                (oftype(ν, 4.09e9) - ν²)
            N = N_report
            T = oftype(ν, 273.15)
        else
            # Egan (2026), Eq. (8), 20 °C and 100 kPa:
            # 1e8(n - 1) = 5627886.8/(281.0 - s²) + 197258.6/(45.4 - s²).
            r = oftype(ν, 5.6278868e14) / (oftype(ν, 2.81e10) - ν²) +
                oftype(ν, 1.972586e13) / (oftype(ν, 4.54e9) - ν²)
            N = oftype(ν, _RI_N_EGAN_O2_CM3)
            T = oftype(ν, 293.15)
        end
    elseif gas === :CO2
        r = oftype(ν, 1.1427e11) * (
            oftype(ν, 5799.25) / (oftype(ν, 128908.9)^2 - ν²) +
            oftype(ν, 120.05) / (oftype(ν, 89223.8)^2 - ν²) +
            oftype(ν, 5.3334) / (oftype(ν, 75037.5)^2 - ν²) +
            oftype(ν, 4.3244) / (oftype(ν, 67837.7)^2 - ν²) +
            oftype(ν, 1.218145e-5) /
                (oftype(ν, 2418.136)^2 - ν²))
        N = oftype(ν, _RI_N_288K_CM3)
        N_report = N
    else
        error("unreachable Rayleigh species $gas")
    end
    return (refractivity = r, number_density = N, report_density = N_report,
            temperature_K = T)
end

function _rayleigh_king_factor(gas::Symbol, ν)
    ν² = ν^2
    if gas === :He || gas === :Ar
        return one(ν)
    elseif gas === :N2
        return oftype(ν, 1.034) + oftype(ν, 3.17e-12) * ν²
    elseif gas === :O2
        # He et al. (2021), Table 1 (kept by decision 13a, 2026-10-19). Bates
        # (1984) and Bodhaine et al. (1999) use the constant 1.096, which
        # raises the O2 cross section by 0.54-0.55 %.
        return oftype(ν, 1.09) + oftype(ν, 1.385e-11) * ν² +
               oftype(ν, 1.448e-20) * ν²^2
    elseif gas === :CO2
        return oftype(ν, 1.1364) + oftype(ν, 2.53e-11) * ν²
    end
    error("unreachable Rayleigh species $gas")
end

# ---------------------------------------------------------------------------
# Polarizability-tensor method
# ---------------------------------------------------------------------------

"""
    UsingPolTen(anisotropy::AbstractPolarizabilityDispersion = DispersiveAnisotropy())

Select the polarizability-tensor Rayleigh formulation, built from the mean
polarizability ``\\bar\\alpha`` and the anisotropy ``\\gamma`` of N₂, O₂ or H₂
(Unicode aliases such as `:N₂` are accepted). `anisotropy` selects their
wavelength dependence:

- [`DispersiveAnisotropy`](@ref) (default): ``\\bar\\alpha(\\nu)`` and
  ``\\gamma(\\nu)`` from the [`UsingRefInd`](@ref) refractivity and King factor
  of N₂ and O₂, with the Buldakov temperature factor on ``\\bar\\alpha``. This
  is the formulation of the RT model's automatic (`depol: -1`) Rayleigh
  depolarization.
- [`ConstantAnisotropy`](@ref): the Raman-module tensor (Buldakov
  ``\\bar\\alpha(\\nu, T)``, wavelength-independent ``\\gamma``), which the Raman
  and Cabannes cross sections use. `UsingPolTen(ConstantAnisotropy())`
  reproduces the `UsingPolTen()` results of the v2.2 pre-releases.

Both go through the Raman-module routines `getMolecularConstants`,
`compute_effective_coefficents!` and `compute_σ_Rayl_coeff!` and never mutate
shared state. See [`molecular_rayleigh_properties`](@ref) for the formulation.
"""
struct UsingPolTen{A<:AbstractPolarizabilityDispersion} <: AbstractMolecularRayleighMethod
    "Wavelength dependence of the tensor invariants."
    anisotropy::A
end

UsingPolTen() = UsingPolTen(DispersiveAnisotropy())

# Reference number density (molecule cm⁻³) at 288.15 K and 1013.25 hPa, the
# same standard density used by the RI fits for He, Ar, N2 and CO2.
const _POLTEN_DEFAULT_NUMBER_DENSITY_CM3 = 2.546899e19
const _POLTEN_DEFAULT_TEMPERATURE_K = 288.15

"""
    _poltensor_rayleigh_coefficients(gas::Symbol, ν_cm1::Float64, T_K::Float64,
                                     [anisotropy = ConstantAnisotropy()])

Internal hook implemented in `InelasticScattering`
(`src/Inelastic/molecular_rayleigh_poltensor.jl`). It returns the effective
Raman-module coefficients `(α̅, γ̅, γ_C, ρ, σ_coeff)` for canonical species
`gas ∈ (:N2, :O2, :H2)` at wavenumber `ν_cm1` (cm⁻¹) and temperature `T_K`.
The three-argument form is the Raman-module (constant-anisotropy) tensor.
The function is declared here so that `Scattering` does not need to load the
Raman module at definition time.
"""
function _poltensor_rayleigh_coefficients end

"""
    _dispersive_polarizability_tensor(gas, ν_cm1, T_K, α_b, α_c) -> (α̅, γ̅)

Mean polarizability and anisotropy (cm³) of [`DispersiveAnisotropy`](@ref) for
`gas ∈ (:N2, :O2)` at wavenumber `ν_cm1` (cm⁻¹) and temperature `T_K`:
``\\bar\\alpha = 3L/(4\\pi N)\\,f(T)/f(T_\\mathrm{fit})`` and
``\\gamma = 3L/(4\\pi N)\\,\\sqrt{9(F_K-1)/2}``, with ``L = (n^2-1)/(n^2+2)``,
``N`` and ``T_\\mathrm{fit}`` of the `UsingRefInd` dispersion fit used at
`ν_cm1`, its King factor ``F_K``, and the Buldakov temperature factor
``f(T) = 1 + \\alpha_b T + \\alpha_c T^2`` of the Raman-module constants.
Evaluated in `Float64`.
"""
function _dispersive_polarizability_tensor(gas::Symbol, ν_cm1::Float64,
                                           T_K::Float64, α_b::Float64,
                                           α_c::Float64)
    gas in (:N2, :O2) || throw(ArgumentError(
        "DispersiveAnisotropy() has refractive-index data for N2 and O2 only; got $gas"))
    fit = _rayleigh_refractivity(gas, ν_cm1)
    x = fit.refractivity * 1.0e-8
    n²_minus_one = x * (2 + x)
    L = n²_minus_one / (3 + n²_minus_one)
    α̅_fit = 3L / (4π * fit.number_density)
    f(t) = 1 + α_b * t + α_c * t^2
    F = _rayleigh_king_factor(gas, ν_cm1)
    return (α̅ = α̅_fit * f(T_K) / f(fit.temperature_K),
            γ̅ = α̅_fit * sqrt(9 * (F - 1) / 2))
end

"""
    molecular_rayleigh_properties(method::UsingPolTen, species, wavelength_nm;
                                  temperature_K = 288.15,
                                  number_density_cm3 = 2.546899e19)

Polarizability-tensor variant. Species: `:N2`, `:O2`, `:H2`. For any other
species an `ArgumentError` is thrown that names the species and suggests
`UsingRefInd()`.

Tensor invariants, selected by `method.anisotropy`:

- [`DispersiveAnisotropy`](@ref) (`UsingPolTen()`), N₂ and O₂:
  ``\\bar\\alpha(\\nu,T) = 3L/(4\\pi N)\\,f(T)/f(T_\\mathrm{fit})`` and
  ``\\gamma(\\nu) = 3L/(4\\pi N)\\sqrt{9(F_K-1)/2}`` from the `UsingRefInd`
  Lorentz--Lorenz factor ``L``, density ``N`` and King factor ``F_K`` at
  ``\\nu = 10^7/\\lambda`` cm⁻¹, with the Buldakov temperature factor
  ``f(T) = 1 + \\alpha_b T + \\alpha_c T^2``. At ``T = T_\\mathrm{fit}`` (288.15 K
  for N₂ below 436.3 nm, 273.15 K above; 273.15 K for O₂ below 466.7 nm,
  293.15 K above) the cross section and depolarization equal the `UsingRefInd`
  values; at other temperatures only ``\\bar\\alpha`` changes (by 1.8×10⁻⁶ K⁻¹
  for N₂ and about 2.6×10⁻⁶ K⁻¹ for O₂ near 288 K). H₂ uses the
  constant-anisotropy tensor below.
- [`ConstantAnisotropy`](@ref) (Raman-module constants of
  `getMolecularConstants(::N₂/::O₂/::H₂, vmr)`): mean polarizability with
  the Buldakov et al. single-pole dispersion and temperature dependence,
  ``\\bar\\alpha(\\nu, T) = \\bar\\alpha_{00}(1 + \\alpha_b T + \\alpha_c T^2)
  / (1 - (2\\pi c \\nu/\\omega_0)^2)``, and the anisotropy
  ``\\gamma = \\bar\\gamma_{00}``, **independent of wavelength and temperature**.

Derived quantities (no further inputs):

- `depolarization` ``\\rho = 6\\gamma^2/(45\\bar\\alpha^2 + 7\\gamma^2)``
  (the module's `rho_depol_Rayl`, the same convention as the RI method and
  [`get_greek_rayleigh`](@ref));
- `king_factor` ``F_K = 1 + 2\\gamma^2/(9\\bar\\alpha^2)
  = (6+3\\rho)/(6-7\\rho)``;
- `cross_section_cm2`
  ``\\sigma = \\tfrac{128}{3}\\pi^5\\nu^4\\bar\\alpha^2 F_K``, i.e. the
  module's `σ_Rayl_coeff · ν⁴` (cm² molecule⁻¹, ``\\bar\\alpha`` in cm³,
  ``\\nu`` in cm⁻¹). This is the total (Cabannes plus rotational Raman)
  molecular Rayleigh cross section, as for the RI method;
- `refractive_index`: the Lorentz--Lorenz equivalent
  ``(n^2-1)/(n^2+2) = 4\\pi N\\bar\\alpha/3`` at `number_density_cm3`
  (default 2.546899e19 molecule cm⁻³, 288.15 K and 1013.25 hPa). Pass
  `number_density_cm3 = nothing` to return `refractive_index = nothing`.
  The cross section, depolarization and King factor are per molecule and do
  not depend on `number_density_cm3`.

The default `temperature_K = 288.15` K (kept by decision 13d, 2026-10-19) is
the temperature of the default `number_density_cm3`, the standard state of the
He, Ar, N₂ (λ < 436 nm) and CO₂ RI fits and of the Bates (1984) and Bodhaine et
al. (1999) cross sections. It therefore compares the tensor with those data at
their own state; the O₂ RI refractive index refers to 273.15 K and
2.68678e19 molecule cm⁻³, so compare O₂ refractive indices at that density.
The temperature enters only through the small Buldakov factor on
``\\bar\\alpha``: between 200 and 300 K the cross sections change by less
than 4×10⁻⁴ (N₂, O₂) and 2×10⁻³ (H₂). The RT model's own
molecular-depolarization path does not use this default: it evaluates the
tensor at the column-weighted mean temperature of the profile and the
band-centre wavenumber. Calling this function does not change any model
default.

The internal evaluation is always in `Float64` (``\\bar\\alpha^2 \\sim 10^{-48}``
underflows `Float32`); results are converted to the floating-point type of
`float(wavelength_nm)`. The method is therefore not ForwardDiff-aware (dual
numbers are not propagated) and does not gain precision from `BigFloat`
inputs. As for the RI method, a `Float32` `refractive_index` resolves ``n-1``
only to about `eps(1f0)/(n-1)` relative (up to ~1e-3, for H₂). The
Raman-module polarizability constants and their temperature coefficients come
from Buldakov et al.; the O₂ and H₂ anisotropies from Asawaroengchai &
Rosenblatt (1980), per the source comments. The H₂ tensor is checked against
the Peck & Huang (1977) refractivity and the Raj, Hamaguchi & Witek (2018)
anisotropy in `test/test_molecular_rayleigh.jl`.
"""
function molecular_rayleigh_properties(method::UsingPolTen, species,
                                        wavelength_nm::Real;
                                        temperature_K::Real =
                                            _POLTEN_DEFAULT_TEMPERATURE_K,
                                        number_density_cm3::Union{Nothing,Real} =
                                            _POLTEN_DEFAULT_NUMBER_DENSITY_CM3)
    gas = _canonical_poltensor_species(species)
    λ = float(wavelength_nm)
    isfinite(λ) && λ > zero(λ) ||
        throw(ArgumentError("wavelength_nm must be finite and positive; got $wavelength_nm"))
    T = Float64(temperature_K)
    isfinite(T) && T > 0 ||
        throw(ArgumentError("temperature_K must be finite and positive; got $temperature_K"))
    if number_density_cm3 !== nothing
        N = Float64(number_density_cm3)
        isfinite(N) && N > 0 ||
            throw(ArgumentError("number_density_cm3 must be finite and positive or nothing; got $number_density_cm3"))
    end

    FT = typeof(λ)
    ν = 1.0e7 / Float64(λ)
    c = _poltensor_rayleigh_coefficients(gas, ν, T, method.anisotropy)
    c.α̅ > 0 || throw(DomainError(c.α̅,
        "mean polarizability must be positive; wavelength $(λ) nm is at or beyond the dispersion pole of the $gas tensor constants"))

    king_factor = 3 * (1 + 2 * c.γ_C) / (3 - 4 * c.γ_C)
    cross_section_cm2 = c.σ_coeff * ν^4

    refractive_index = if number_density_cm3 === nothing
        nothing
    else
        # Lorentz–Lorenz: L = (n² - 1)/(n² + 2) = 4πNα̅/3, so
        # n² - 1 = 3L/(1 - L) and n - 1 = (n² - 1)/(n + 1).
        L = 4π * Float64(number_density_cm3) * c.α̅ / 3
        L < 1 || throw(DomainError(L, "Lorentz–Lorenz factor must be below one"))
        n²_minus_one = 3L / (1 - L)
        FT(1 + n²_minus_one / (sqrt(1 + n²_minus_one) + 1))
    end

    return (
        species = gas,
        wavelength_nm = λ,
        refractive_index = refractive_index,
        king_factor = FT(king_factor),
        depolarization = FT(c.ρ),
        cross_section_cm2 = FT(cross_section_cm2),
    )
end

function molecular_rayleigh_cross_section_ratio(method::UsingPolTen, species,
                                                 wavelength_nm::Real,
                                                 reference_nm::Real;
                                                 temperature_K::Real =
                                                     _POLTEN_DEFAULT_TEMPERATURE_K)
    σ = molecular_rayleigh_properties(method, species, wavelength_nm;
                                      temperature_K = temperature_K,
                                      number_density_cm3 = nothing).cross_section_cm2
    σ_ref = molecular_rayleigh_properties(method, species, reference_nm;
                                          temperature_K = temperature_K,
                                          number_density_cm3 = nothing).cross_section_cm2
    return σ / σ_ref
end

function _canonical_poltensor_species(species)
    key = species isa Symbol ? species : Symbol(species)
    key in (:N2, :N₂, :nitrogen, :Nitrogen) && return :N2
    key in (:O2, :O₂, :oxygen, :Oxygen) && return :O2
    key in (:H2, :H₂, :hydrogen, :Hydrogen) && return :H2
    throw(ArgumentError(
        "UsingPolTen(): no polarizability tensor is available for species " *
        "$species; tensors exist only for N2, O2 and H2. Use UsingRefInd() " *
        "(refractive index + King factor) for He, Ar, N2, O2 or CO2."))
end
