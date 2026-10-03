"""
    molecular_rayleigh_properties(species, wavelength_nm)

Return elastic Rayleigh properties for a pure gas at vacuum wavelength
`wavelength_nm` (nm). Supported species are `:He`, `:Ar`, `:N2`, `:O2`, and
`:CO2` (Unicode aliases such as `:N₂` are accepted).

The returned named tuple contains

- `refractive_index`: refractive index at the source reference number density;
- `king_factor`: molecular-anisotropy correction ``F_K``;
- `depolarization`: depolarization ratio ``rho`` in the convention consumed by
  [`get_greek_rayleigh`](@ref), related by
  ``F_K = (6 + 3rho)/(6 - 7rho)``;
- `cross_section_cm2`: Rayleigh scattering cross section in cm² molecule⁻¹.

The dispersion and King-factor fits are those collected and experimentally
checked over the UV-visible by He et al. (2021), with the Ar fit from Wilmouth
and Sayres (2019). Their cross-section expression is

```math
sigma_R = \\frac{24 pi^3 nu^4}{N^2}
          \\left(\\frac{n^2-1}{n^2+2}\\right)^2 F_K,
```

where `nu` is in cm⁻¹ and `N` in molecule cm⁻³. The implementation is
intended for the 360--400 nm comparison and has been validated only over that
interval in vSmartMOM. It describes elastic Rayleigh scattering; it does not
add species-specific rotational or vibrational Raman redistribution.

# References

- He et al. (2021), *Atmospheric Chemistry and Physics* 21, 14927--14940,
  https://doi.org/10.5194/acp-21-14927-2021.
- Wilmouth and Sayres (2019), *Atmospheric Measurement Techniques* 12,
  1277--1293, https://doi.org/10.5194/amt-12-1277-2019.
- Sneep and Ubachs (2005), *JQSRT* 92, 293--310,
  https://doi.org/10.1016/j.jqsrt.2004.07.025.
"""
function molecular_rayleigh_properties(species, wavelength_nm::Real)
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

"""
    molecular_rayleigh_cross_section_ratio(species, wavelength_nm, reference_nm)

Return ``sigma_R(wavelength_nm) / sigma_R(reference_nm)`` for a pure gas.
This ratio is the appropriate scaling when a comparison fixes the Rayleigh
optical depth at one reference wavelength and varies molecular identity.
"""
function molecular_rayleigh_cross_section_ratio(species, wavelength_nm::Real,
                                                 reference_nm::Real)
    σ = molecular_rayleigh_properties(species, wavelength_nm).cross_section_cm2
    σ_ref = molecular_rayleigh_properties(species, reference_nm).cross_section_cm2
    return σ / σ_ref
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
        return oftype(ν, 1.09) + oftype(ν, 1.385e-11) * ν² +
               oftype(ν, 1.448e-20) * ν²^2
    elseif gas === :CO2
        return oftype(ν, 1.1364) + oftype(ν, 2.53e-11) * ν²
    end
    error("unreachable Rayleigh species $gas")
end
