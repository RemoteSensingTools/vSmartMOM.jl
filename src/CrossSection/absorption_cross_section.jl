using KernelAbstractions
using CUDA
using Parameters
using ..Architectures: devi

"""
Given the hitran data and necessary parameters, calculate an absorption cross-section at the given pressure, 
temperature, and grid of wavelengths (or wavenumbers)
"""
function compute_absorption_cross_section(
                #Required
                hitran::HitranTable,          # Model to use in this cross section calculation 
                                              # (Calculation from Hitran data vs. using Interpolator)
                broadening::AbstractBroadeningFunction, # Broadening function to use
                grid::Array{<:Real,1},        # Wavelength [nm] or wavenumber [cm-1] grid
                pressure::Real,               # actual pressure [hPa]
                temperature::Real;            # actual temperature [K] 
                # Optionals
                wavelength_flag::Bool=false,  # Use wavelength in nm (true) or wavenumber cm-1 units (false)
                wing_cutoff::Real=40,         # Wing cutoff [cm -1]
                vmr::Real=0,                  # VMR of gas [0-1]
                CEF::AbstractComplexErrorFunction=ErfcErrorFunction(), # Error function to use (if Voigt broadening)
                architecture::AbstractArchitecture=Architectures.CPU()
                )

    # Store results here to return
    gridC = array_type(architecture)(grid);

    result = similar(gridC);
    fill!(result,0);

    # Convert to wavenumber from [nm] space if necessary
    grid = wavelength_flag ? reverse(nm_per_m ./ grid) : grid

    # Calculate the minimum and maximum grid bounds, including the wing cutoff
    grid_max = maximum(grid) + wing_cutoff
    grid_min = minimum(grid) - wing_cutoff

    # Interpolators from grid bounds to index values
    grid_idx_interp_low  = LinearInterpolation(grid, 1:1:length(grid),extrapolation_bc = 1)
    grid_idx_interp_high = LinearInterpolation(grid, 1:1:length(grid),extrapolation_bc = length(grid))

    # Temporary storage array for output of qoft!. Compiler/speed issues when returning value in qoft
    rate = zeros(eltype(AbstractFloat(temperature)),1)

    # Declare the device being used
    device = devi(architecture)

    # Loop through all transition lines:
    for j in eachindex(hitran.Sᵢ)

        # Test that this ν lies within the grid
        # (Aside: I ♥︎ chained comparisons)
        if grid_min < hitran.νᵢ[j] < grid_max

            # Apply pressure shift
            ν   = hitran.νᵢ[j] + pressure/p_ref*hitran.δ_air[j]

            # Compute Lorentzian HWHM
            γ_l = (hitran.γ_air[j] *
                  (1-vmr)*pressure/p_ref+hitran.γ_self[j] *
                  vmr*pressure/p_ref) *
                  (t_ref/temperature)^hitran.n_air[j]

            # Compute Doppler HWHM
            γ_d = ((cSqrt2Ln2/cc_)*sqrt(cBolts_/cMassMol)*sqrt(temperature) * 
            hitran.νᵢ[j]/sqrt(mol_weight(hitran.mol[j],hitran.iso[j])))

            # Ratio of widths
            y = sqrt(cLn2) * γ_l/γ_d

            # Apply line intensity temperature corrections
            S = hitran.Sᵢ[j]
            if hitran.E″[j] != -1
                qoft!(hitran.mol[j],hitran.iso[j],temperature,t_ref, rate)
                S = S * rate[1] *
                        exp(c₂*hitran.E″[j]*(1/t_ref-1/temperature)) *
                        (1-exp(-c₂*hitran.νᵢ[j]/temperature))/(1-exp(-c₂*hitran.νᵢ[j]/t_ref));

            end

            # Calculate index range that this ν impacts
            ind_start = Int64(round(grid_idx_interp_low(ν - wing_cutoff)))
            ind_stop  = Int64(round(grid_idx_interp_high(ν + wing_cutoff)))
            
            # Create views from the result and grid arrays
            result_view   = view(result ,ind_start:ind_stop);
            grid_view     = view(gridC   ,ind_start:ind_stop);

            # Kernel for performing the lineshape calculation
            kernel! = line_shape!(device)

            # Run the event on the kernel 
            # That this, this function adds to each element in result, the contribution from this transition
            event = kernel!(result_view, grid_view, ν, γ_d, γ_l, y, S, broadening, CEF, ndrange=length(grid_view))
            wait(device, event)
        end
    end

    # Return the resulting lineshape
    return result
end

"""
    $(FUNCTIONNAME)(model::HitranModel, grid::AbstractRange{<:Real}, pressure::Real, temperature::Real; wavelength_flag::Bool=false)

Given a HitranModel, return the calculated absorption cross-section at the given pressure, 
temperature, and grid of wavelengths (or wavenumbers)

"""
function absorption_cross_section(
                # Required
                model::HitranModel,          # Model to use in this cross section calculation 
                                             # (Calculation from Hitran data vs. using Interpolator)
                grid::AbstractRange{<:Real}, # Wavelength [nm] or wavenumber [cm-1] grid (modify using wavelength_flag)
                pressure::Real,              # actual pressure [hPa]
                temperature::Real,           # actual temperature [K]    
                # Optionals
                wavelength_flag::Bool=false  # Use wavelength in nm (true) or wavenumber cm-1 units (false)       
                )
    FT = typeof(AbstractFloat(model.vmr))
    return compute_absorption_cross_section(model.hitran, model.broadening, collect(grid), pressure, temperature, wavelength_flag=wavelength_flag, wing_cutoff=model.wing_cutoff, vmr=model.vmr, CEF=model.CEF, architecture=model.architecture)

end

"""
    $(FUNCTIONNAME)(model::InterpolationModel, grid::AbstractRange{<:Real}, pressure::Real, temperature::Real; wavelength_flag::Bool=false)

Given an Interpolation Model, return the interpolated absorption cross-section at the given pressure, 
temperature, and grid of wavelengths (or wavenumbers)

"""
function absorption_cross_section(
    # Required
    model::InterpolationModel,      # Model to use in this cross section calculation 
                                    # (Calculation from Interpolator vs. Hitran Data)
    grid::AbstractRange{<:Real},    # Wavelength [nm] or wavenumber [cm-1] grid
    pressure::Real,                 # actual pressure [hPa]
    temperature::Real,              # actual temperature [K]  
    #Optionals 
    wavelength_flag::Bool=false,    # Use wavelength in nm (true) or wavenumber cm-1 units (false)            
    )

    # Convert to wavenumber from [nm] space if necessary
    grid = wavelength_flag ? reverse(nm_per_m ./ collect(grid)) : collect(grid)

    # Scale the interpolation to match the model grids
    sitp = scale(model.itp, model.p_grid, model.t_grid, model.ν_grid)

    # Perform the interpolation and return the resulting grid
    return sitp(pressure, temperature, grid)
end

#####
##### Lineshape functions that are called by absorption_cross_section
#####

@kernel function line_shape!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Doppler, CEF)
    FT = eltype(ν)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(cSqrtLn2divSqrtPi) * exp(-FT(cLn2) * ((FT(grid[I]) - FT(ν)) / FT(γ_d)) ^2) / FT(γ_d)
end

@kernel function line_shape!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Lorentz, CEF)
    FT = eltype(ν)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(γ_l) / (FT(pi) * (FT(γ_l) ^2 + (FT(grid[I]) - FT(ν)) ^ 2))
end

@kernel function line_shape!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Voigt, CEF)
    FT = eltype(ν)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(cSqrtLn2divSqrtPi)/FT(γ_d) * real(w(CEF, FT(cSqrtLn2) / FT(γ_d) * (FT(grid[I]) - FT(ν)) + im * FT(y)))
end

@kernel function line_shape32!(A, @Const(grid), ν, γ_d, γ_l, y, S, broadening, CEF)
    line_shape!(A, grid, Float32(ν), Float32(γ_d), Float32(γ_l), Float32(y), Float32(S), broadening, CEF)
end