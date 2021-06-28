module SolarModel


function planck_spectrum(T::Real, grid::Vector; wavelength_flag=false)

    h = 6.626070e-34    # J⋅Hz−1
    c = 3e8             # m / s
    k = 1.380649e-23    # J⋅K−1
    
    #    ν?            ν -> λ                 ν: cm⁻¹ to m⁻¹
    νs = wavelength_flag ? (1 ./ (1e-9 .* grid)) : grid * 100 

    L = map(ν-> (2*h*c^2) * (1/((1/ν)^5*(exp(h*c/((1/ν)*k*T)) - 1))) , νs)

    return [grid L]

end

function output_irradiance(L_solar, R)

    
end

export planck_spectrum

end