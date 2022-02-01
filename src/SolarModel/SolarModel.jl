module SolarModel

using ..vSmartMOM               # For locating default solar T
using DocStringExtensions       # For simplifying docstring
using DelimitedFiles            # For easily reading in solar spectrum 
using Interpolations            # For interpolating solar spectrum
using Pkg.Artifacts             # For default solar spectrum

"""
    $(FUNCTIONNAME)(T::Real, ν_grid::Vector)

Produce the black-body planck spectrum (mW/m²-sr-cm⁻¹), given the temperature (K) 
and calculation grid (ν in cm⁻¹)

"""
function planck_spectrum_wn(T::Real, ν_grid::Vector)

    c1 = 1.1910427 * 10^(-5)    # mW/m²-sr-cm⁻¹
    c2 = 1.4387752              # K⋅cm

    # L(ν, T) = c1⋅ν³/(exp(c2⋅ν/T) - 1)
    radiance = c1 .* (ν_grid.^3) ./ (exp.(c2 * ν_grid / T) .- 1)

    return radiance
end

"""
    $(FUNCTIONNAME)(T::Real, λ_grid::Vector)

Produce the black-body planck spectrum (W/m²-sr-μm), given the temperature (K) 
and calculation grid (λ in μm)

"""
function planck_spectrum_wl(T::Real, λ_grid::Vector)

    c1 = 1.1910427 * 10^8    # W/m²-sr-μm
    c2 = 1.4387752 * 10^4    # K⋅μm

    # L(ν, T) = c1⋅ν³/(exp(c2⋅ν/T) - 1)
    radiance = c1 ./ (λ_grid.^5 .* (exp.(c2 ./ (λ_grid * T)) .- 1))

    return radiance
end

# W/m²-sr-μm to Ph/s-m²-sr-um
# λ_grid in micron
function watts_to_photons(λ_grid::Vector, radiance::Vector)

    h = 6.62607015e-34 # J⋅Hz−1
    c = 299792458 # m/s

    E_per_λ = h * c ./ (λ_grid / 1e6)
    photons = radiance ./ E_per_λ

    return photons
end

"""
    $(FUNCTIONNAME)(T::Real; stride_length::Integer = 100)

Produce the black-body planck spectrum (mW/m²-sr-cm⁻¹), given the temperature (K). 
Use a unit calculation grid and check for convergence every `stride_length` cm⁻¹ until the 
spectrum dies off. 

"""
function planck_spectrum_wn(T::Real; stride_length::Integer = 100)

    # νs, starting with ν0 = 1.0 cm⁻¹
    νs = [1.0]

    # radiances corresponding with νs
    radiances = planck_spectrum_wn(T, νs)

    # Loop until convergence
    while true 
        
        # Add the next ν
        νs = vcat(νs, collect(νs[end] + 1 : νs[end] + stride_length))

        # Compute the next radiance
        radiances = vcat(radiances, planck_spectrum_wn(T, νs[(end - stride_length + 1) : end]))

        # Exit if spectrum has died off
        (radiances[end] < radiances[1]) && break 

    end

    return [νs[1:(end-1)] radiances[1:(end-1)]]
end

"""
    $(FUNCTIONNAME)(solar, ν_grid)

Interpolate a solar linelist to the ν_grid
"""
function itp_solar_to_ν_grid(solar, ν_grid)

    solar_idx_start = maximum((argmin(abs.(solar[:, 1] .- minimum(ν_grid))) - 10, 1))
    solar_idx_end   = minimum((argmin(abs.(solar[:, 1] .- maximum(ν_grid))) + 10, length(solar[:,1])))

    solar_subset = solar[solar_idx_start:solar_idx_end, :]

    itp = LinearInterpolation(solar_subset[:, 1], 
                              solar_subset[:, 2])

    return itp.(ν_grid)
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
    return itp_solar_to_ν_grid(solar, ν_grid)
end

"""
    $(FUNCTIONNAME)(ν_grid::Union{AbstractRange{<:Real}, AbstractArray} = 600.0:0.01:26316.0)

Get the default solar transmission and interpolate to wavenumber grid (entire grid if not specified)
"""
function default_solar_transmission(ν_grid::Union{AbstractRange{<:Real}, AbstractArray} = 600.0:0.01:26316.0)
    
    @info "Using line-list from:\nToon, G. C., Solar line list for GGG2014, TCCON data archive, hosted by the Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, Tennessee, U.S.A., doi:10. 14291/tccon.ggg2014.solar.R0/1221658, 2014."
    @info "Found at: https://mark4sun.jpl.nasa.gov/toon/solar/solar_spectrum.html"

    filename = joinpath(dirname(pathof(RadiativeTransfer)), "SolarModel", "solar.out")
    !isfile(filename) && download("http://web.gps.caltech.edu/~cfranken/hitran_2016/solar_merged_20160127_600_26316_100.out", filename)

    return hcat(ν_grid, solar_transmission_from_file(filename, ν_grid))
end

"""
    $(FUNCTIONNAME)(ν_grid::Union{AbstractRange{<:Real}, AbstractArray} = 600.0:0.01:26316.0)

Get the default solar spectrum and interpolate to wavenumber grid (entire grid if not specified)
"""
function default_solar_spectrum_at_earth(ν_grid::Union{AbstractRange{<:Real}, AbstractArray} = 600.0:0.01:26316.0)

    T = 5777 # K
    λ_grid = reverse(1e4 ./ ν_grid) # Wavenumber to micron
    black_body = reverse(SolarModel.watts_to_photons(λ_grid, planck_spectrum_wl(T, λ_grid) * 2.1629e-05 * pi))
    solar_transmission = default_solar_transmission(ν_grid)[:,2]

    return hcat(ν_grid, black_body .* solar_transmission)
end

export planck_spectrum_wn, planck_spectrum_wl, solar_transmission_from_file, default_solar_transmission

end