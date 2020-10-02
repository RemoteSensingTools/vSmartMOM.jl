
""" Convenience function to create Mie Model with NAI2 computation type """
function make_mie_model(computation_type::NAI2, 
                        aerosol::AbstractAerosolType, 
                        λ::Real,
                        polarization::AbstractPolarizationType, 
                        truncation_type::AbstractTruncationType)
    return MieModel(computation_type, aerosol, λ, polarization, truncation_type, zeros(1, 1, 1), zeros(1, 1, 1))
end

""" Convenience function to load Wigner matrices and create Mie Model with PCW computation type """
function make_mie_model(computation_type::PCW, 
                        aerosol::AbstractAerosolType, 
                        λ::Real,
                        polarization::AbstractPolarizationType, 
                        truncation_type::AbstractTruncationType, 
                        wigner_filepath::String)
    wigner_A, wigner_B = PhaseFunction.load_wigner_values(wigner_filepath)
    return MieModel(computation_type, aerosol, λ, polarization, truncation_type, wigner_A, wigner_B)
end