#=

Convenience alias: `absorption_cross_section` forwards to `compute_absorption_cross_section`.

=#

"""
    absorption_cross_section(model, grid, pressure, temperature; wavelength_flag=false)

Calculate absorption cross-section. Convenience wrapper around `compute_absorption_cross_section`.

# Arguments
- `model::AbstractCrossSectionModel`: HitranModel or InterpolationModel
- `grid`: Wavelength [nm] or wavenumber [cm⁻¹] grid (see wavelength_flag)
- `pressure::Real`: Pressure (hPa)
- `temperature::Real`: Temperature (K)
- `wavelength_flag::Bool=false`: If true, grid is wavelength in nm; else wavenumber in cm⁻¹

# Returns
- `Vector`: Absorption cross-sections (cm²/molecule) at each grid point
"""
absorption_cross_section(model::AbstractCrossSectionModel,
                         grid::Union{AbstractRange{<:Real}, AbstractArray},
                         pressure::Real,
                         temperature::Real;
                         wavelength_flag::Bool=false) =
    compute_absorption_cross_section(model, grid, pressure, temperature; wavelength_flag=wavelength_flag)

absorption_cross_section(model::AtmosphericAbsorption.LineByLineModel,
                         grid::Union{AbstractRange{<:Real}, AbstractArray},
                         pressure::Real,
                         temperature::Real) =
    AtmosphericAbsorption.compute_cross_section(model, grid, pressure, temperature)
