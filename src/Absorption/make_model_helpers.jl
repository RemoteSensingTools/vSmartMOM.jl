#=
 
This file contains helper functions for creating calculation models 

`make_hitran_model` creates a HitranModel with all specified parameters. 

`make_interpolation_model` creates a cross-section matrix based on the pressure and 
temperature grids, and fits an interpolation. 

There are also some save/load convenience functions for the interpolation model. 

=#

# The standalone line-by-line engine is deprecated in favour of AtmosphericAbsorption.jl,
# which the RT pipeline already uses. Kept through vSmartMOM 2.x; removal planned for 3.0.
const _LEGACY_LINE_ENGINE_MIGRATION = """
vSmartMOM's standalone line-by-line model (`make_hitran_model`, `HitranModel`) is
deprecated and will be removed in vSmartMOM 3.0. The RT pipeline already computes
line absorption with AtmosphericAbsorption.jl (add it with
`Pkg.add(name="AtmosphericAbsorption", version="0.1.2")`). The corresponding calls:

    import AtmosphericAbsorption as AA
    lines = AA.load_lines(AA.HitranPort(artifact("CO2")); mol=2, iso=1,
                          ν_min=6000, ν_max=6400)
    model = AA.LineByLineModel(lines; profile=AA.Voigt(), wing_cutoff=40, vmr=0,
                               cpf=AA.HumlicekWeideman32(), architecture=AA.CPU())
    σ     = AA.compute_cross_section(model, ν_grid, p_hPa, T_K)

The two engines agree closely for vmr = 0. For vmr > 0 the differences grow with
vmr and can reach the percent level: AtmosphericAbsorption mixes the pressure shift,
(1 - vmr) δ_air + vmr δ_self, whereas this engine applies δ_air unchanged. Partition
sums are TIPS-2021 there and TIPS-2017 here. Validate migrated results. Tables for YAML `LUTfiles` are still built
with `make_interpolation_model` in vSmartMOM 2.x. If you publish results computed
from HITRAN data, please cite the HITRAN edition you used (see the vSmartMOM
References page).
"""

"Emit the legacy line-engine deprecation notice once per session."
_warn_legacy_line_engine() = @warn _LEGACY_LINE_ENGINE_MIGRATION maxlog = 1

"""
    make_hitran_model(hitran::HitranTable, 
                      broadening::AbstractBroadeningFunction; 
                      wing_cutoff::Real=40, 
                      vmr::Real=0, 
                      CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(), 
                      architecture::AbstractArchitecture = default_architecture())

Convenience function to make a HitranModel out of the parameters (Matches make_interpolation_model)

!!! warning "Deprecated"
    The standalone line-by-line model is deprecated and will be removed in vSmartMOM 3.0;
    use `AtmosphericAbsorption.LineByLineModel` (see `load_lines`, `compute_cross_section`).
    A one-time notice with the corresponding calls and their convention differences is
    logged on first use.
"""
function make_hitran_model(hitran::HitranTable, 
                           broadening::AbstractBroadeningFunction; 
                           wing_cutoff::Real=40, 
                           vmr::Union{Real, Vector}=0, 
                           CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(), 
                           architecture::AbstractArchitecture = default_architecture())

    _warn_legacy_line_engine()
    if architecture isa GPU && !(CEF isa HumlicekWeidemann32SDErrorFunction)
        @warn "Cross-section calculations on GPU may or may not work with this CEF (use HumlicekWeidemann32SDErrorFunction if you encounter issues)"
    end

    return HitranModel(hitran=hitran, broadening=broadening , wing_cutoff=wing_cutoff , vmr=vmr, CEF=CEF, architecture=architecture)
end

# The interpolation table stores its spectral axis as an `AbstractRange`, which only a
# uniform wavenumber grid satisfies; wavelength input is converted at query time instead.
function _require_wavenumber_lut_grid(wavelength_flag::Bool)
    wavelength_flag && throw(ArgumentError(
        "make_interpolation_model builds tables on a uniform wavenumber grid (cm⁻¹); " *
        "wavelength_flag=true is not supported. Pass a uniform wavenumber grid whose ends " *
        "lie outside the band, e.g. range(floor(1e7/λ_max_nm/Δν)*Δν, " *
        "ceil(1e7/λ_min_nm/Δν)*Δν; step=Δν), and query the table in wavelength " *
        "with absorption_cross_section(model, λ_nm, p, T; wavelength_flag=true)."))
    return nothing
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
                    architecture::AbstractArchitecture = default_architecture())

Using a HitranModel, create an InterpolationModel by interpolating the lineshape function 
at the given pressure and temperature grids. `vmr` is the absorber's volume mixing
ratio passed to the underlying line model; in this engine it mixes only the Lorentz
widths (self vs. air), while the pressure shift and temperature exponent stay air-only.
The resulting table is specific to that value.

The table is built on a uniform **wavenumber** grid (cm⁻¹). `wavelength_flag=true` is
rejected: a uniform wavelength grid maps to a non-uniform wavenumber grid, which the
interpolation table cannot store. Build the table on a uniform wavenumber grid whose ends
lie outside the band, e.g.
`wave_grid = range(floor(1e7/λ_max_nm/Δν)*Δν, ceil(1e7/λ_min_nm/Δν)*Δν; step=Δν)`, and
query it in wavelength with `absorption_cross_section(model, λ_nm, p, T;
wavelength_flag=true)`, which remains supported.

!!! note "Line engine"
    The tables are computed with the deprecated standalone line-by-line model. They remain
    a supported format for YAML `LUTfiles` throughout vSmartMOM 2.x. Tables saved with
    AtmosphericAbsorption's `save_interpolation_model` use a different format that this
    loader does not read; AtmosphericAbsorption `.absco` tables (`save_absco_lut`) are read.
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
                                  wing_cutoff::Real=40,
                                  vmr::Real=0,
                                  CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(),
                                  architecture::AbstractArchitecture = default_architecture(),
                                  # Linear avoids cubic spline overshoot/undershoot on under-sampled
                                  # narrow absorption lines. Default stays :cubic for back-compat.
                                  interp_type::Symbol = :cubic,
                                  # Storage type of the LUT coefficient array. Float32 halves disk/VRAM
                                  # with negligible accuracy loss for σ values that span ~10 orders
                                  # of magnitude logarithmically but are interpolated linearly.
                                  storage_type::Type{<:AbstractFloat} = Float64,
                                )

    # Warn user if using incompatible/untested CEF
    if architecture isa Architectures.GPU && !(CEF isa HumlicekWeidemann32SDErrorFunction)
        @warn "Cross-section calculations on GPU may or may not work with this CEF (use HumlicekWeidemann32SDErrorFunction if you encounter issues)"
    end

    _require_wavenumber_lut_grid(wavelength_flag)
    ν_grid = wave_grid

    # Empty matrix to store the calculated cross-sections
    cs_matrix = zeros(storage_type, length(ν_grid), length(p_grid), length(t_grid))

    # One line model (with the requested self-broadening vmr) for all (p, T) nodes
    model = HitranModel(hitran=hitran, broadening=broadening, wing_cutoff=wing_cutoff,
                        vmr=vmr, CEF=CEF, architecture=architecture)
    @showprogress 1 "Computing Cross Sections for Interpolation..." for i in 1:length(p_grid)
        for j in 1:length(t_grid)
            cs_matrix[:,i,j] = storage_type.(collect(compute_absorption_cross_section(model, collect(ν_grid), p_grid[i], t_grid[j])))
        end
    end

    # Perform the interpolation
    itp = if interp_type === :linear
        interpolate(cs_matrix, BSpline(Linear()))
    elseif interp_type === :cubic
        interpolate(cs_matrix, BSpline(Cubic(Line(OnGrid()))))
    else
        throw(ArgumentError("interp_type must be :linear or :cubic, got $(interp_type)"))
    end

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

"""
    make_interpolation_model(absco::AbscoTable, wave_grid, p_grid, t_grid; 
                             wavelength_flag=false, architecture=default_architecture())

Create an InterpolationModel from ABSCO tabulated cross-section data.

Linearly interpolates the pre-computed cross-sections onto the specified pressure and
temperature grids. Uses the ABSCO table layout σ(ν, broadener VMR, T, p).
Only the first broadener slice (`σ[:, 1, :, :]`, H₂O VMR = 0 in OCO-2 ABSCO
tables) is used, so the result has no H₂O broadening. Listing the original
`.hdf` file in `absorption.LUTfiles` keeps the H₂O-broadener axis instead.

# Arguments
- `absco::AbscoTable`: ABSCO table with cross-sections σ(ν, broadener VMR, T, p)
- `wave_grid`: uniform wavenumber grid in cm⁻¹; `wavelength_flag=true` is rejected
  (see the HITRAN method)
- `p_grid`: Pressure grid (hPa)
- `t_grid`: Temperature grid (K)

# Returns
- `InterpolationModel`: Interpolator for absorption cross-sections
"""
function make_interpolation_model(
        # Required
        absco::AbscoTable,  
        wave_grid::AbstractRange{<:Real}, 
        p_grid::AbstractRange{<:Real},
        t_grid::AbstractRange{<:Real}; 
        # Optionals
        wavelength_flag::Bool=false,
        architecture::AbstractArchitecture = default_architecture()     # Computer `Architecture` on which `Model` is run
    )

    _require_wavenumber_lut_grid(wavelength_flag)
    ν_grid = wave_grid

    # Empty matrix to store the calculated cross-sections
    cs_matrix = zeros(length(ν_grid),length(p_grid), length(t_grid));
    
    # Define indices in Absco:
    T = absco.T
    p = absco.p
    xs = absco.σ
    index_p = 1:length(p)
    index_t = 1:size(T,1)
    inter_p = linear_interpolation(p, collect(index_p), extrapolation_bc = Flat())
    # Calculate all the cross-sections at the pressure and temperature grids
    @showprogress 1 "Computing Cross Sections for Interpolation..." for i in 1:length(p_grid)
        p_ref = p_grid[i]
        fractional_index_p = inter_p(p_ref)
        inter_T_top    = linear_interpolation(T[:,Int(ceil(fractional_index_p))] , collect(index_t), extrapolation_bc = Flat())
        inter_T_bottom = linear_interpolation(T[:,Int(floor(fractional_index_p))], collect(index_t), extrapolation_bc = Flat())
      
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
            interp_xs = linear_interpolation(absco.ν, xs_interp, extrapolation_bc = Flat())
            cs_matrix[:,i,j] = collect(interp_xs(ν_grid))
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

"""
    make_interpolation_model_test(absco::AbscoTable, wave_grid, p_grid, t_grid; 
                                  wavelength_flag=false, architecture=default_architecture())

Alternative InterpolationModel from ABSCO using cubic spline interpolation in temperature.

Uses cubic splines per pressure level for smoother T-interpolation, then linear
interpolation on the combined (ν, p, T) grid. Assumes equidistant T-grid in ABSCO.

# Arguments
- `absco::AbscoTable`: ABSCO table with cross-sections
- `wave_grid`, `p_grid`, `t_grid`: Output grids (same as make_interpolation_model)

# Returns
- `InterpolationModel`: Interpolator for absorption cross-sections
"""
function make_interpolation_model_test(
        # Required
        absco::AbscoTable,  
        wave_grid::AbstractRange{<:Real}, 
        p_grid::AbstractRange{<:Real},
        t_grid::AbstractRange{<:Real}; 
        # Optionals
        wavelength_flag::Bool=false,
        architecture::AbstractArchitecture = default_architecture()     # Computer `Architecture` on which `Model` is run
    )

    _require_wavenumber_lut_grid(wavelength_flag)
    ν_grid = wave_grid

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
    full_interp = linear_interpolation((absco_grid,p, t_grid),cs_matrix_temp, extrapolation_bc = Flat());
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

