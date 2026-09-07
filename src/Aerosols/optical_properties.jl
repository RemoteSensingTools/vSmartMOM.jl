"""
    compute_optical_properties(data::AerosolData, wavelengths::AbstractVector,
                               ri_database::RefractiveIndexDatabase)

Reserved entry point for converting ingested aerosol schemes to physical optics.
Currently throws `ArgumentError`: the TOMAS and two-moment ingestion framework
has no validated conversion from its concentration units and vertical geometry
to extinction, single-scattering albedo, and phase coefficients.

Earlier versions returned placeholder efficiencies, number densities, and a
constant asymmetry parameter. Those values were not suitable for scientific
use and are no longer returned. Use the production `Scattering` Mie APIs and
CoreRT aerosol configuration for radiative-transfer calculations. A future
scheme adapter must supply explicit units, meteorology, size distributions,
and a validated Mie/phase calculation before implementing this interface.
"""
function compute_optical_properties(data::AerosolData, wavelengths::AbstractVector,
                                    ri_database::RefractiveIndexDatabase)
    throw(ArgumentError("Aerosols.compute_optical_properties is not implemented: " *
        "scheme-to-optics conversion requires validated units, meteorology, and Mie optics. " *
        "Use vSmartMOM.Scattering and CoreRT aerosol configuration."))
end

"""
    compute_mie_efficiencies(x::Real, n::Complex)

Unsupported legacy helper. Throws `ArgumentError` instead of returning the
former interpolation between unvalidated efficiency approximations. Use the
Mie implementation in `vSmartMOM.Scattering`.
"""
function compute_mie_efficiencies(x::Real, n::Complex)
    throw(ArgumentError("Aerosols.compute_mie_efficiencies was a placeholder; " *
                        "use the validated vSmartMOM.Scattering Mie APIs."))
end

"""
    integrate_phase_function(data::AerosolData, wavelength::Real,
                             ri_database::RefractiveIndexDatabase,
                             scattering_angles::AbstractVector)

Reserved scheme-to-phase-function adapter. Throws `ArgumentError`; no phase
integration for ingested aerosol schemes is implemented. The former fixed
Henyey–Greenstein curve did not depend on the input aerosol state.
"""
function integrate_phase_function(data::AerosolData, wavelength::Real,
                                   ri_database::RefractiveIndexDatabase,
                                   scattering_angles::AbstractVector)
    throw(ArgumentError("Aerosols.integrate_phase_function is not implemented; " *
                        "use vSmartMOM.Scattering phase calculations."))
end
