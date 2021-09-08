#####
##### Functions to deal with a cross-section interpolator
#####

using ..Architectures: GPU

"""
    $(FUNCTIONNAME)(hitran::HitranTable, 
                    broadening::AbstractBroadeningFunction, 
                    wave_grid::AbstractRange{<:Real}, 
                    p_grid::AbstractRange{<:Real},
                    t_grid::AbstractRange{<:Real}; 
                    wavelength_flag::Bool=false,
                    wing_cutoff::Real=40, 
                    vmr::Real=0, 
                    CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(),
                    architecture::AbstractArchitecture=default_architecture )

Using a HitranModel, create an InterpolationModel by interpolating the lineshape function 
at the given pressure and temperature grids

"""
function make_interpolation_model(
                                  # Required
                                  hitran::HitranTable, 
                                  broadening::AbstractBroadeningFunction, 
                                  wave_grid::AbstractRange{<:Real}, 
                                  p_grid::AbstractRange{<:Real},
                                  t_grid::AbstractRange{<:Real}; 
                                  # Optionals
                                  wavelength_flag::Bool=false,
                                  wing_cutoff::Integer=40, 
                                  vmr::Real=0, 
                                  CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(),
                                  architecture::AbstractArchitecture=default_architecture     # Computer `Architecture` on which `Model` is run
                                )

    # Warn user if using incompatible/untested CEF
    if architecture isa GPU && !(CEF isa HumlicekWeidemann32SDErrorFunction)
        @warn "Cross-section calculations on GPU may or may not work with this CEF (use HumlicekWeidemann32SDErrorFunction if you encounter issues)"
    end

    # Convert from wavelength to wavenumber if necessary
    ν_grid = wavelength_flag ? reverse(nm_per_m ./ wave_grid) : wave_grid

    # Empty matrix to store the calculated cross-sections
    cs_matrix = zeros(length(p_grid), length(t_grid), length(ν_grid));

    # Calculate all the cross-sections at the pressure and temperature grids
    @showprogress 1 "Computing Cross Sections for Interpolation..." for i in 1:length(p_grid)
        for j in 1:length(t_grid)
            # make_hitran_model(hitran_data, Voigt(), wing_cutoff = 40, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=CPU())
            model = make_hitran_model(hitran, broadening, wing_cutoff=wing_cutoff, CEF=CEF, architecture=architecture)
            cs_matrix[i,j,:] = Array(compute_absorption_cross_section(model, collect(ν_grid), p_grid[i], t_grid[j], wavelength_flag=wavelength_flag))
        end
    end
    
    # Perform the interpolation
    itp = interpolate(cs_matrix, BSpline(Cubic(Line(OnGrid()))))

    # Get the molecule and isotope numbers from the HitranTable
    mol = hitran.mol[1]
    iso = all(x->x==hitran.iso[1], hitran.iso) ? hitran.iso[1] : -1

    # Return the interpolation model with all the proper parameters
    return InterpolationModel(itp, mol, iso, broadening, ν_grid, p_grid, t_grid, wing_cutoff, vmr, CEF, architecture)
end

""" Convenience function to save an InterpolationModel at a specified filepath (Using JLD2) """
function save_interpolation_model(itp_model::InterpolationModel, filepath::String)
    @save filepath itp_model
end

""" Convenience function to load an InterpolationModel from a specified filepath (Using JLD2) """
function load_interpolation_model(filepath::String)
    itp_model::InterpolationModel = @load filepath itp_model
    return itp_model
end
