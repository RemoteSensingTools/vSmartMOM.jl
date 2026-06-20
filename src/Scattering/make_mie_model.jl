#=
 
This file provides helper functions to create a MieModel object

=#

"""
    make_mie_model(::NAI2, aerosol::Aerosol, λ, polarization, truncation_type, r_max, nquad_radius;
                   architecture = Architectures.CPU(), precision_policy = nothing) -> MieModel

Construct a [`MieModel`](@ref) configured for the Siewert NAI-2 workflow.

# Arguments
- `aerosol`: aerosol size-distribution and refractive-index specification.
- `λ`: wavelength (must use the same length units as `r_max` and the aerosol radius scale).
- `polarization`: one of [`Stokes_I`](@ref), [`Stokes_IQ`](@ref),
  [`Stokes_IQU`](@ref), [`Stokes_IQUV`](@ref).
- `truncation_type`: typically [`δBGE`](@ref).
- `r_max`: upper radius used in size quadrature.
- `nquad_radius`: number of radius quadrature points.

# Keyword arguments
- `architecture`: `Architectures.CPU()` (default) or `Architectures.GPU()`.
  Selects the CPU vs GPU Mie path at
  [`compute_aerosol_optical_properties`](@ref) time.
- `precision_policy`: GPU `Dₙ`-recursion precision (`NativeFloat64`/`DSEmulated`).
  `nothing` (default) auto-selects on the GPU and is ignored on the CPU.

# Returns
- `MieModel` ready for [`compute_aerosol_optical_properties`](@ref).
"""
function make_mie_model(computation_type::NAI2,
                        aerosol::Aerosol{FT},
                        λ::Real,
                        polarization::AbstractPolarizationType,
                        truncation_type::AbstractTruncationType,
                        r_max::Real,
                        nquad_radius::Integer;
                        architecture::Architectures.AbstractArchitecture = Architectures.CPU(),
                        precision_policy = nothing) where {FT}
    # Promote λ and r_max to the aerosol's element type so the MieModel is
    # FT-consistent: aerosol::Aerosol{FT}, λ::FT, r_max::FT all agree.
    # Trivial placeholder wigner arrays are allocated in FT (NAI2 never reads them).
    FT_r = typeof(FT(r_max))
    return MieModel(computation_type, aerosol, FT_r(λ), polarization, truncation_type,
                    FT_r(r_max), nquad_radius,
                    zeros(FT_r, 1, 1, 1), zeros(FT_r, 1, 1, 1),
                    architecture, precision_policy)
end

"""
    make_mie_model(::PCW, aerosol::Aerosol, λ, polarization, truncation_type, r_max, nquad_radius, wigner_filepath::String) -> MieModel

Construct a [`MieModel`](@ref) configured for the Domke PCW workflow, loading
Wigner tables from `wigner_filepath` via [`load_wigner_values`](@ref).
"""
function make_mie_model(computation_type::PCW,
                        aerosol::Aerosol{FT},
                        λ::Real,
                        polarization::AbstractPolarizationType,
                        truncation_type::AbstractTruncationType,
                        r_max::Real,
                        nquad_radius::Integer,
                        wigner_filepath::String;
                        architecture::Architectures.AbstractArchitecture = Architectures.CPU(),
                        precision_policy = nothing) where {FT}
    wigner_A_raw, wigner_B_raw = Scattering.load_wigner_values(wigner_filepath)
    # Convert Wigner tables to the model's element type for FT-consistency.
    FT_r = typeof(FT(r_max))
    wigner_A = convert(Array{FT_r,3}, wigner_A_raw)
    wigner_B = convert(Array{FT_r,3}, wigner_B_raw)
    return MieModel(computation_type, aerosol, FT_r(λ), polarization, truncation_type,
                    FT_r(r_max), nquad_radius, wigner_A, wigner_B, architecture, precision_policy)
end

"""
    make_mie_model(::PCW, aerosol::Aerosol, λ, polarization, truncation_type, r_max, nquad_radius, wigner_A, wigner_B) -> MieModel

Construct a [`MieModel`](@ref) configured for the Domke PCW workflow, using
precomputed Wigner tables `wigner_A` and `wigner_B`.

Use this overload when Wigner tensors are already in memory.
"""
function make_mie_model(computation_type::PCW,
                        aerosol::Aerosol{FT},
                        λ::Real,
                        polarization::AbstractPolarizationType,
                        truncation_type::AbstractTruncationType,
                        r_max::Real,
                        nquad_radius::Integer,
                        wigner_A, wigner_B;
                        architecture::Architectures.AbstractArchitecture = Architectures.CPU(),
                        precision_policy = nothing) where {FT}
    # Convert Wigner tables to the model's element type for FT-consistency.
    FT_r = typeof(FT(r_max))
    return MieModel(computation_type, aerosol, FT_r(λ), polarization, truncation_type,
                    FT_r(r_max), nquad_radius,
                    convert(Array{FT_r,3}, wigner_A),
                    convert(Array{FT_r,3}, wigner_B),
                    architecture, precision_policy)
end
