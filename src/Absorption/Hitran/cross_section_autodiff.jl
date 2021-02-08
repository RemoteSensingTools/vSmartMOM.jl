
"""
    $(FUNCTIONNAME)(model::HitranModel, grid::AbstractRange{<:Real}, pressure::Real, temperature::Real; wavelength_flag::Bool=false)

Calculate absorption cross-section at the given pressure, temperature, and grid of wavelengths 
(or wavenumbers), and have the option to perform auto-differentiation

"""
function absorption_cross_section(model::AbstractCrossSectionModel,          # Model to use 
                                  grid::AbstractRange{<:Real}, # Wavelength [nm] or wavenumber [cm-1] grid 
                                  pressure::Real,              # actual pressure [hPa]
                                  temperature::Real,           # actual temperature [K]    
                                  wavelength_flag::Bool=false ; autodiff=false)

    # This function takes in the "x-vector" along with the input model so that ForwardDiff will work
    function absorption_cross_section_autodiff(x ; model::AbstractCrossSectionModel = model, grid=grid, pressure=pressure, 
                                               temperature=temperature, wavelength_flag=wavelength_flag)

        if length(x) !== 2
            @error "Must receive two cross-section parameters for auto-differentiation (p, T)" x
        end
    
        # Make sure that 𝐱 and parameter match
        @assert (pressure== x[1])
        @assert (temperature == x[2])
    
        cross_section = absorption_cross_section(model, grid, x[1], x[2]);
    
        return cross_section
    end

    if (autodiff)

        x = [pressure, temperature]

        result = DiffResults.JacobianResult(zeros(length(collect(grid))), x);
        ForwardDiff.jacobian!(result, absorption_cross_section_autodiff, x);

        return (result.value, result.derivs[1]);

    else 
        return absorption_cross_section(model, grid, pressure, temperature, wavelength_flag)
    end

end