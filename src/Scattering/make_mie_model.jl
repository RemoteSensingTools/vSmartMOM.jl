#=
 
This file provides helper functions to create a MieModel object

=#

"""
    make_mie_model(::NAI2, aerosol::Aerosol, λ, polarization, truncation_type, r_max, nquad_radius) -> MieModel

Construct a [`MieModel`](@ref) configured for the Siewert NAI-2 workflow.

# Arguments
- `aerosol`: aerosol size-distribution and refractive-index specification.
- `λ`: wavelength (must use the same length units as `r_max` and the aerosol radius scale).
- `polarization`: one of [`Stokes_I`](@ref), [`Stokes_IQ`](@ref),
  [`Stokes_IQU`](@ref), [`Stokes_IQUV`](@ref).
- `truncation_type`: typically [`δBGE`](@ref).
- `r_max`: upper radius used in size quadrature.
- `nquad_radius`: number of radius quadrature points.

# Returns
- `MieModel` ready for [`compute_aerosol_optical_properties`](@ref).
"""
function make_mie_model(computation_type::NAI2, 
                        aerosol::Aerosol, 
                        λ::Real,
                        polarization::AbstractPolarizationType, 
                        truncation_type::AbstractTruncationType,
                        r_max::Real,
                        nquad_radius::Integer)
    return MieModel(computation_type, aerosol, λ, polarization, truncation_type, r_max, nquad_radius, zeros(1, 1, 1), zeros(1, 1, 1))
end

"""
    make_mie_model(::PCW, aerosol::Aerosol, λ, polarization, truncation_type, r_max, nquad_radius, wigner_filepath::String) -> MieModel

Construct a [`MieModel`](@ref) configured for the Domke PCW workflow, loading
Wigner tables from `wigner_filepath` via [`load_wigner_values`](@ref).
"""
function make_mie_model(computation_type::PCW, 
                        aerosol::Aerosol, 
                        λ::Real,
                        polarization::AbstractPolarizationType, 
                        truncation_type::AbstractTruncationType, 
                        r_max::Real,
                        nquad_radius::Integer,
                        wigner_filepath::String)
    wigner_A, wigner_B = Scattering.load_wigner_values(wigner_filepath)
    return MieModel(computation_type, aerosol, λ, polarization, truncation_type, r_max, nquad_radius, wigner_A, wigner_B)
end

"""
    make_mie_model(::PCW, aerosol::Aerosol, λ, polarization, truncation_type, r_max, nquad_radius, wigner_A, wigner_B) -> MieModel

Construct a [`MieModel`](@ref) configured for the Domke PCW workflow, using
precomputed Wigner tables `wigner_A` and `wigner_B`.

Use this overload when Wigner tensors are already in memory.
"""
function make_mie_model(computation_type::PCW, 
                        aerosol::Aerosol, 
                        λ::Real,
                        polarization::AbstractPolarizationType, 
                        truncation_type::AbstractTruncationType, 
                        r_max::Real,
                        nquad_radius::Integer,
                        wigner_A, wigner_B)
return MieModel(computation_type, aerosol, λ, polarization, truncation_type, r_max, nquad_radius, wigner_A, wigner_B)
end
