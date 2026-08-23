#=

Convenience alias: `absorption_cross_section` forwards to `compute_absorption_cross_section`.

=#

"""
    absorption_cross_section(model, grid, pressure, temperature;
                             wavelength_flag=false)

Calculate absorption cross-section. Convenience wrapper around `compute_absorption_cross_section`.

# Arguments
- `model::AbstractCrossSectionModel`: HitranModel or InterpolationModel
- `grid`: Wavelength [nm] or wavenumber [cm⁻¹] grid (see wavelength_flag)
- `pressure::Real`: Pressure (hPa)
- `temperature::Real`: Temperature (K)
- `wavelength_flag::Bool=false`: If true, grid is wavelength in nm; else wavenumber in cm⁻¹
- `vmr=model.vmr`: Additional keyword accepted by the
  `AtmosphericAbsorption.LineByLineModel` dispatch. This is the absorber's
  moist-air mole fraction used to mix foreign and self broadening.

# Returns
- `Vector`: Absorption cross-sections (cm²/molecule) at each grid point
"""
absorption_cross_section(model::AbstractCrossSectionModel,
                         grid::Union{AbstractRange{<:Real}, AbstractArray},
                         pressure::Real,
                         temperature::Real;
                         wavelength_flag::Bool=false) =
    compute_absorption_cross_section(model, grid, pressure, temperature; wavelength_flag=wavelength_flag)

# AA dispatch: `wavelength_flag=true` takes a grid in nm and returns σ on that
# nm grid in input order (the nm→cm⁻¹ conversion and wing-cutoff windowing
# happen inside AA, in wavenumber space). Mirrors the legacy method above.
# Both `wavelength_flag` and the per-call broadener `vmr` are forwarded to AA.
absorption_cross_section(model::AtmosphericAbsorption.LineByLineModel,
                         grid::Union{AbstractRange{<:Real}, AbstractArray},
                         pressure::Real,
                         temperature::Real;
                         wavelength_flag::Bool=false,
                         vmr::Real=model.vmr) =
    wavelength_flag ?
        AtmosphericAbsorption.compute_cross_section(model, grid, pressure, temperature;
                                                    wavelength_flag=true, vmr=vmr) :
        AtmosphericAbsorption.compute_cross_section(model, grid, pressure, temperature;
                                                    vmr=vmr)

# Native ABSCO tables are linear in spectral position, pressure, temperature, and their tabulated
# H₂O-broadener axis. `vmr` here is the H₂O broadener VMR, not the absorber abundance used later to
# turn σ into optical depth. Keep this public standalone route aligned with the profile path below.
absorption_cross_section(model::AtmosphericAbsorption.AbscoLUT,
                         grid::Union{AbstractRange{<:Real}, AbstractArray},
                         pressure::Real,
                         temperature::Real;
                         wavelength_flag::Bool=false,
                         vmr::Real=0) =
    AtmosphericAbsorption.compute_cross_section(
        model,
        wavelength_flag ? 1e7 ./ grid : grid,
        pressure,
        temperature;
        vmr,
        interp=:linear,
    )
