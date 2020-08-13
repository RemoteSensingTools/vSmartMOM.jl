#####
##### Functions to deal with a cross-section interpolator
#####

function make_interpolation_model(model::HitranModel, 
                                 wave_grid::Array{<:Real,1}, 
                                 wavelength_flag::Bool,
                                 p_grid::Array{<:Real,1},
                                 t_grid::Array{<:Real,1}
                                )

    ν_grid = wavelength_flag ? nm_per_m ./ wave_grid : wave_grid
    cs_matrix = zeros(length(p_grid), length(t_grid), length(ν_grid));

    @showprogress 1 "Computing Cross Sections..." for i in 1:length(p_grid)
        for j in 1:length(t_grid)
            println(t_grid[j])
            @time cs_matrix[i,j,:] = absorption_cross_section(model, ν_grid, false, p_grid[i], t_grid[j])
        end
    end
    
    itp = interpolate(cs_matrix, BSpline(Cubic(Line(OnGrid()))))

    mol = model.hitran.mol[1]
    iso = all(x->x==model.hitran.iso[1], model.hitran.iso) ? model.hitran.iso[1] : -1

    if model.broadening isa Doppler
        broadening = "Doppler"
    elseif model.broadening isa Lorentz
        broadening = "Lorentz"
    elseif model.broadening isa Voigt
        broadening = "Voigt"
    end

    wing_cutoff = model.wing_cutoff
    vmr = model.vmr

    return InterpolationModel(itp, mol, iso, broadening, ν_grid, p_grid, t_grid, wing_cutoff, vmr, model.CEF)
end

function save_interpolation_model(itp_model::InterpolationModel, filepath::String)
    @save filepath itp_model
end

function load_interpolation_model(filepath::String)
    itp_model::InterpolationModel = @load filepath itp
    return itp_model
end