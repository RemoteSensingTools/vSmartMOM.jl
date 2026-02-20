"""
Shared physical constants for Raman/inelastic scattering calculations.
Used by inelastic_helper.jl, apply_lineshape.jl, and related modules.
"""
const c₂                 = 1.4387769
const cSqrtLn2divSqrtPi  = 0.469718639319144059835  # √(ln2/π)
const cLn2               = 0.6931471805599          # ln2
const cSqrtLn2           = 0.8325546111577          # √(ln2)
const cSqrt2Ln2          = 1.1774100225             # √(2ln2)
const cc_                = 2.99792458e8             # speed of light [m/s]
const cBolts_            = 1.3806503e-23            # Boltzmann const. [J/K]
const p_ref              = 1013.25                  # reference pressure [hPa]
const t_ref              = 296.0                    # reference temperature [K]
const nm_per_m           = 1.0e7                    # conversion: nm to cm⁻¹ (wavenumber)
# Molecular mass constants for unit molecular mass (1 amu)
const cMassMol           = 1.66053873e-27           # kg per molecule (used in inelastic_helper apply_lineshape!)
const cMassMolIE         = 1.66053873e-30           # kg per molecule (used in apply_lineshape.jl apply_lineshape_!)
