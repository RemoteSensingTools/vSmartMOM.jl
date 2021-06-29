module SolarModel

using DocStringExtensions       # For simplifying docstring
using DelimitedFiles            # For easily reading in solar spectrum 
using Interpolations            # For interpolating solar spectrum

"""
    $(FUNCTIONNAME)(T::Real, grid::Vector; wavelength_flag=false)

Produce the black-body planck spectrum, given the temperature and calculation grid

"""
function planck_spectrum(T::Real, grid::Vector; wavelength_flag=false)

    h = 6.626070e-34    # J⋅Hz−1
    c = 3e8             # m / s
    k = 1.380649e-23    # J⋅K−1
    
    #    ν?            ν -> λ                 ν: cm⁻¹ to m⁻¹
    νs = wavelength_flag ? (1 ./ (1e-9 .* grid)) : grid * 100 

    L = map(ν-> (2*h*c^2) * (1/((1/ν)^5*(exp(h*c/((1/ν)*k*T)) - 1))) , νs)

    return [grid L]

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