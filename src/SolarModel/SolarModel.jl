module SolarModel

using DocStringExtensions       # For simplifying docstring
using DelimitedFiles            # For easily reading in solar spectrum 
using Interpolations            # For interpolating solar spectrum

"""
    $(FUNCTIONNAME)(T::Real, grid::Vector; wavelength_flag=false)

Produce the black-body planck spectrum (mW/m²-sr-cm⁻¹), given the temperature (K) 
and calculation grid (ν in cm⁻¹; if wavelength_flag=true, λ in nm)

"""
function planck_spectrum(T::Real, grid::Vector; wavelength_flag=false)

    c1 = 1.1910427 * 10^(-5)    # mW/m²-sr-cm⁻¹
    c2 = 1.4387752              # K⋅cm

    # Convert to wavenumbers if given in wavelengths
    #    λ?                 λ -> ν            ν
    νs = wavelength_flag ? (1e7 ./ (grid)) : grid

    # L(ν, T) = c1⋅ν³/(exp(c2⋅ν/T) - 1)
    radiance = c1 .* (νs.^3) ./ (exp.(c2 * νs / T) .- 1)

    return radiance
end

"""
    $(FUNCTIONNAME)(T::Real; stride_length::Integer = 100)

Produce the black-body planck spectrum (mW/m²-sr-cm⁻¹), given the temperature (K). 
Use a unit calculation grid and check for convergence every `stride_length` cm⁻¹ until the 
spectrum dies off. 

"""
function planck_spectrum(T::Real; stride_length::Integer = 100)

    # νs, starting with ν0 = 1.0 cm⁻¹
    νs = [1.0]

    # radiances corresponding with νs
    radiances = planck_spectrum(T, νs)

    # Loop until convergence
    while true 
        
        # Add the next ν
        νs = vcat(νs, collect(νs[end] + 1 : νs[end] + stride_length))

        # Compute the next radiance
        radiances = vcat(radiances, planck_spectrum(T, νs[(end - stride_length + 1) : end]))

        # Exit if spectrum has died off
        (radiances[end] < radiances[1]) && break 

    end

    return [νs[1:(end-1)] radiances[1:(end-1)]]
end

"""
    $(FUNCTIONNAME)(file_name::String)

Get the solar transmission from the specified file
"""
solar_transmission_from_file(file_name::String) = readdlm(file_name)

"""
    $(FUNCTIONNAME)(file_name::String, ν_grid::Union{AbstractRange{<:Real}, AbstractArray})

Get the solar transmission from the specified file, and interpolate to wavenumber grid
"""
function solar_transmission_from_file(file_name::String, 
                                    ν_grid::Union{AbstractRange{<:Real}, AbstractArray})

    solar = solar_transmission_from_file(file_name)

    solar_idx_start = argmin(abs.(solar[:, 1] .- minimum(ν_grid)))
    solar_idx_end   = argmin(abs.(solar[:, 1] .- maximum(ν_grid)))

    solar_subset = solar[(solar_idx_start-10):(solar_idx_end+10), :]

    itp = LinearInterpolation(solar_subset[:, 1], 
                              solar_subset[:, 2])

    return itp.(ν_grid)
end

export planck_spectrum, solar_transmission_from_file

end