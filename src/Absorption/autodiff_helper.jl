#=
 
This file wraps `compute_absorption_cross_section` so that autodiff users and non-autodiff 
users can call the same function with just a keyword argument change. 

(Tried doing this with multiple dispatch but ran into keyword argument issues)

=#

"""
    $(FUNCTIONNAME)(model::HitranModel, grid::AbstractRange{<:Real}, pressure::Real, temperature::Real; wavelength_flag::Bool=false)

Calculate absorption cross-section at the given pressure, temperature, and grid of wavelengths 
(or wavenumbers), and have the option to perform auto-differentiation

"""
function absorption_cross_section(model::AbstractCrossSectionModel,          # Model to use 
                                  grid::Union{AbstractRange{<:Real}, AbstractArray}, # Wavelength [nm] or wavenumber [cm-1] grid 
                                  pressure::Real,              # actual pressure [hPa]
                                  temperature::Real;
                                  autodiff::Bool=false,           # actual temperature [K]    
                                  wavelength_flag::Bool=false)
                                  

    # This function takes in the "x-vector" along with the input model so that ForwardDiff will work
    function absorption_cross_section_autodiff(x ; model::AbstractCrossSectionModel = model, grid=grid, pressure=pressure, 
                                               temperature=temperature, wavelength_flag=wavelength_flag)

        if length(x) !== 2
            @error "Must receive two cross-section parameters for auto-differentiation (p, T)" x
        end
    
        # Make sure that 𝐱 and parameter match
        @assert (pressure == x[1])
        @assert (temperature == x[2])
    
        return compute_absorption_cross_section(model, grid, x[1], x[2], wavelength_flag=wavelength_flag);
    end

    if (autodiff)

        x = [pressure, temperature]

        result = DiffResults.JacobianResult(zeros(length(collect(grid))), x);
        ForwardDiff.jacobian!(result, absorption_cross_section_autodiff, x);

        return (result.value, result.derivs[1]);

    else 
        return compute_absorption_cross_section(model, grid, pressure, temperature, wavelength_flag=wavelength_flag)
    end

end