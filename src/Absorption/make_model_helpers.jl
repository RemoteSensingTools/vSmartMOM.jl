#=
 
This file contains helper functions for creating calculation models 

`make_hitran_model` creates a HitranModel with all specified parameters. 

`make_interpolation_model` creates a cross-section matrix based on the pressure and 
temperature grids, and fits an interpolation. 

There are also some save/load convenience functions for the interpolation model. 

=#

"""
    make_hitran_model(hitran::HitranTable, 
                      broadening::AbstractBroadeningFunction; 
                      wing_cutoff::Real=40, 
                      vmr::Real=0, 
                      CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(), 
                      architecture = default_architecture)

Convenience function to make a HitranModel out of the parameters (Matches make_interpolation_model)

"""
function make_hitran_model(hitran::HitranTable, 
                           broadening::AbstractBroadeningFunction; 
                           wing_cutoff::Integer=40, 
                           vmr::Union{Real, Vector}=0, 
                           CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(), 
                           architecture = default_architecture)

    if architecture isa GPU && !(CEF isa HumlicekWeidemann32SDErrorFunction)
        @warn "Cross-section calculations on GPU may or may not work with this CEF (use HumlicekWeidemann32SDErrorFunction if you encounter issues)"
    end

    return HitranModel(hitran=hitran, broadening=broadening , wing_cutoff=wing_cutoff , vmr=vmr, CEF=CEF, architecture=architecture)
end

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
    if architecture isa Architectures.GPU && !(CEF isa HumlicekWeidemann32SDErrorFunction)
        @warn "Cross-section calculations on GPU may or may not work with this CEF (use HumlicekWeidemann32SDErrorFunction if you encounter issues)"
    end

    # Convert from wavelength to wavenumber if necessary
    ν_grid = wavelength_flag ? reverse(nm_per_m ./ wave_grid) : wave_grid

    # Empty matrix to store the calculated cross-sections
    cs_matrix = zeros(length(ν_grid),length(p_grid), length(t_grid));

    # Calculate all the cross-sections at the pressure and temperature grids
    @showprogress 1 "Computing Cross Sections for Interpolation..." for i in 1:length(p_grid)
        for j in 1:length(t_grid)
            # make_hitran_model(hitran_data, Voigt(), wing_cutoff = 40, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=CPU())
            model = make_hitran_model(hitran, broadening, wing_cutoff=wing_cutoff, CEF=CEF, architecture=architecture)
            cs_matrix[:,i,j] = Array(compute_absorption_cross_section(model, collect(ν_grid), p_grid[i], t_grid[j], wavelength_flag=wavelength_flag))
        end
    end
    
    # Perform the interpolation
    itp = interpolate(cs_matrix, BSpline(Cubic(Line(OnGrid()))))

    # Get the molecule and isotope numbers from the HitranTable
    mol = hitran.mol[1]
    iso = all(x->x==hitran.iso[1], hitran.iso) ? hitran.iso[1] : -1

    # Return the interpolation model with all the proper parameters
    return InterpolationModel(itp, mol, iso,  ν_grid, p_grid, t_grid)
end

""" Convenience function to save an InterpolationModel at a specified filepath (Using JLD2) """
function save_interpolation_model(itp_model::InterpolationModel, filepath::String)
    @save filepath itp_model
end

""" Convenience function to load an InterpolationModel from a specified filepath (Using JLD2) """
function load_interpolation_model(filepath::String)
    @load filepath itp_model
    return itp_model
end

function make_interpolation_model(
        # Required
        absco::AbscoTable,  
        wave_grid::AbstractRange{<:Real}, 
        p_grid::AbstractRange{<:Real},
        t_grid::AbstractRange{<:Real}; 
        # Optionals
        wavelength_flag::Bool=false,
        architecture::AbstractArchitecture=default_architecture     # Computer `Architecture` on which `Model` is run
    )

    # Convert from wavelength to wavenumber if necessary
    ν_grid = wavelength_flag ? reverse(nm_per_m ./ wave_grid) : wave_grid

    # Empty matrix to store the calculated cross-sections
    cs_matrix = zeros(length(ν_grid),length(p_grid), length(t_grid));
    
    # Define indices in Absco:
    T = absco.T
    p = absco.p
    xs = absco.σ
    index_p = 1:length(p)
    index_t = 1:size(T,1)
    inter_p = LinearInterpolation(p, collect(index_p), extrapolation_bc = Flat())
    # Calculate all the cross-sections at the pressure and temperature grids
    @showprogress 1 "Computing Cross Sections for Interpolation..." for i in 1:length(p_grid)
        p_ref = p_grid[i]
        fractional_index_p = inter_p(p_ref)
        inter_T_top    = LinearInterpolation(T[:,Int(ceil(fractional_index_p))] , collect(index_t), extrapolation_bc = Flat())
        inter_T_bottom = LinearInterpolation(T[:,Int(floor(fractional_index_p))], collect(index_t), extrapolation_bc = Flat())
      
        for j in 1:length(t_grid)
            
            t_ref = t_grid[j]

            # Cumbersome linear interpolation of non-equidistant grid in ABSCO!
            fractional_index_T_top    = inter_T_top(t_ref)
            fractional_index_T_bottom = inter_T_bottom(t_ref)
            a1 = ceil(fractional_index_T_top) - fractional_index_T_top
            xs_interp_p1 = (1-a1)*xs[:,1, Int(ceil(fractional_index_T_top)), Int(ceil(fractional_index_p))]    + a1 * xs[:,1,Int(floor(fractional_index_T_top)), Int(ceil(fractional_index_p))]
            
            a1 = ceil(fractional_index_T_bottom)-fractional_index_T_bottom
            xs_interp_p2 = (1-a1)*xs[:,1, Int(ceil(fractional_index_T_bottom)), Int(floor(fractional_index_p))] + a1 * xs[:,1,Int(floor(fractional_index_T_bottom)), Int(floor(fractional_index_p))]
            a1 = ceil(fractional_index_p)-fractional_index_p
            xs_interp = a1*xs_interp_p2+(1-a1)*xs_interp_p1 
            interp_xs = LinearInterpolation(absco.ν, xs_interp, extrapolation_bc = Flat())
            cs_matrix[:,i,j] = Array(interp_xs(ν_grid))
        end
    end

    # Perform the interpolation
    # BSpline(Cubic(Line(OnGrid())))
    # BSpline(Quadratic(Line(OnGrid())))

    itp = interpolate(cs_matrix,  (BSpline(Linear()),BSpline(Quadratic(Line(OnGrid()))),BSpline(Quadratic(Line(OnGrid())))))

    # Get the molecule and isotope numbers from the HitranTable
    mol = absco.mol
    iso = absco.iso

    # Return the interpolation model with all the proper parameters
    return InterpolationModel(itp, mol, iso,  ν_grid, p_grid, t_grid)
end

function make_interpolation_model_test(
        # Required
        absco::AbscoTable,  
        wave_grid::AbstractRange{<:Real}, 
        p_grid::AbstractRange{<:Real},
        t_grid::AbstractRange{<:Real}; 
        # Optionals
        wavelength_flag::Bool=false,
        architecture::AbstractArchitecture=default_architecture     # Computer `Architecture` on which `Model` is run
    )

    # Convert from wavelength to wavenumber if necessary
    ν_grid = wavelength_flag ? reverse(nm_per_m ./ wave_grid) : wave_grid

    # Empty matrix to store the calculated cross-sections
    cs_matrix = zeros(length(ν_grid),length(p_grid), length(t_grid));
    
    # Define indices in Absco:
    T = absco.T
    dT = mean(diff(T, dims=1))
    @assert std(diff(T, dims=1))<1e-4 "T-grid too variable in ABSCO!"
    absco_grid = absco.ν[1]:mean(diff(absco.ν)):absco.ν[end]
    p = absco.p
    xs = absco.σ

    index_p = 1:length(p)
    index_t = 1:size(T,1)

    # Create an array of cubic splines in T (per pressure level in ABSCO)
    # no H2O broadening for now, can be added later:
    cs_matrix_temp = zeros(length(absco.ν),length(p), length(t_grid));
    t_interpolators = [CubicSplineInterpolation((absco_grid,T[1,i]:dT:T[end,i]+0.01),xs[:,1,:,i], extrapolation_bc = Flat()) for i in eachindex(p)]
    for i in eachindex(p)
        #@show i, p[i]
        cs_matrix_temp[:,i,:] = t_interpolators[i](absco.ν,t_grid);
    end
    # final interpolation routine on equidistant grid
    #full_interp = CubicSplineInterpolation((absco_grid,p, t_grid),cs_matrix_temp, extrapolation_bc = Flat());
    # This only works in Linear interpolation so far for irregular grids, need to update at a certain point:
    full_interp = LinearInterpolation((absco_grid,p, t_grid),cs_matrix_temp, extrapolation_bc = Flat());
    cs_matrix[:,:,:] = full_interp(ν_grid, p_grid, t_grid);
    
    # Perform the interpolation (can add a switch whether quadratic or cubic later, ideally in the setup file)
    itp = interpolate(cs_matrix, BSpline(Cubic(Line(OnGrid()))));
    #itp = interpolate(cs_matrix, BSpline(Quadratic(Line(OnGrid()))))
    # interpolate(table, (BSpline(Constant()),BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid())))));

    # Get the molecule and isotope numbers from the HitranTable
    mol = absco.mol;
    iso = absco.iso;

    # Return the interpolation model with all the proper parameters
    return InterpolationModel(itp, mol, iso,  ν_grid, p_grid, t_grid)
end

