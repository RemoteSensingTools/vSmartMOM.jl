module CrossSection

using Parameters                # For constructing HitranTable with keywords
using DocStringExtensions       # For simplifying docstring
using Interpolations            # For interpolating in both lookup tables and qoft!
using JLD2                      # For saving and loading the interpolator
using ProgressMeter             # For showing progress, especially in creating the interpolator
using KernelAbstractions
using CUDA

include("constants.jl")         # Scientific and mathematical constants
include("types.jl")             # All types used in this module
include("hitran.jl")            # HITRAN file-related functions
include("complex_error_functions.jl")   # Complex error functions used in line broadening
include("mol_weights.jl")       # Molecular weights (TODO: replace with netCDF file)
include("TIPS_2017.jl")         # Partition sums data (TODO: replace with netCDF file)
include("partition_sums.jl")    # Partition sums interpolator (TODO: replace with LinearInterpolation)
include("cross_section_interpolator.jl") # Cross-section interpolator functions

# Export the Cross Section models
export HitranModel, InterpolationModel

# Export the absorption_cross_section functions
export compute_absorption_cross_section, absorption_cross_section

# Export the hitran functions from hitran.jl
export read_hitran, make_hitran_model

# Export the broadening function types
export Doppler, Lorentz, Voigt

# Export the complex error functions
export HumlicekErrorFunction, HumlicekWeidemann32VoigtErrorFunction, HumlicekWeidemann32SDErrorFunction, CPF12ErrorFunction, ErfcHumliErrorFunctionVoigt, ErfcHumliErrorFunctionSD, ErfcErrorFunction

# Export the hitran table struct type
export HitranTable

# Export the interpolator functions
export make_interpolation_model, save_interpolation_model, load_interpolation_model

"""
    $(FUNCTIONNAME)(hitran::HitranTable,broadening::AbstractBroadeningFunction,grid::Array{<:Real,1},wavelength_flag::Bool,pressure::Real,temperature::Real,wing_cutoff::Real;vmr::Real=0,CEF::AbstractComplexErrorFunction=ErfcErrorFunction())

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
                CEF::AbstractComplexErrorFunction=ErfcErrorFunction() # Error function to use (if Voigt broadening)
                )

    # Store results here to return
    gridC = CuArray(grid);
    result = similar(grid);
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
    rate = zeros(1)

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

            # Pre factor sqrt(ln(2)/pi)/γ_d
            pre_fac = sqrt(cLn2/pi)/γ_d

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
            
            grid_     = view(grid   ,ind_start:ind_stop);
            result_   = view(result ,ind_start:ind_stop);
            
            device = CPU()
            kernelVoigt! = voigt_shape32!(device,1)
            #kernelCPU = doppler_shape!(CUDADevice(), 1024)
            #γ_d,ν,y,S,
            event = kernelVoigt!(result_, grid_,Float32(γ_d),Float32(ν),Float32(y),Float32(S), CEF, ndrange=length(grid_));wait(device,event)
            #testerCPU!(result_, grid_,γ_d,ν,y,S);
            # Loop over every ν in input grid that this transition impacts
            #= for i in ind_start:ind_stop

                # Depending on the type of input broadening specified, apply that transformation
                # and add to the result array

                # Doppler
                if broadening isa Doppler
                    lineshape_val = cSqrtLn2divSqrtPi*exp(-cLn2*((grid[i] - ν) /γ_d) ^2) /γ_d
                    result[i] += S * lineshape_val

                # Lorentz
                elseif broadening isa Lorentz
                    lineshape_val = γ_l/(pi*(γ_l^2 + (grid[i] - ν) ^ 2))
                    result[i] += S * lineshape_val

                # Voigt
                elseif broadening isa Voigt
                    compl_error = w(CEF, sqrt(cLn2)/γ_d * (grid[i] - ν) + 1im * y);
                    #compl_error = weideman32a(sqrt(cLn2)/γ_d * (grid[i] - ν) + 1im * y);
                    result[i] += pre_fac * S * real(compl_error)
                end
            end =#
        end
    end

    # Return the resulting lineshape
    return Array(result)
end

"""
    $(FUNCTIONNAME)(model::HitranModel, grid::Array{<:Real,1}, wavelength_flag::Bool, pressure::Real, temperature::Real)

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

    return compute_absorption_cross_section(model.hitran, model.broadening, collect(grid), pressure, temperature, wavelength_flag=wavelength_flag, wing_cutoff=model.wing_cutoff, vmr=model.vmr, CEF=model.CEF)

end

"""
    $(FUNCTIONNAME)(model::InterpolationModel, grid::Array{<:Real,1}, wavelength_flag::Bool, pressure::Real, temperature::Real)

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

@kernel function voigt_shape32!(A, @Const(grid),γ_d,ν,y,S, modShape  )
    cSqrtLn2divSqrtPi  = Float32(0.469718639319144059835)
    cLn2               = Float32(0.6931471805599)
    #sqrt(cLn2)/γ_d * (grid[I] - ν) + im * y;
    #γ_d = Float32(0.08);
    #ν = Float32(6050);
    I = @index(Global, Linear)
    @inbounds A[I] += S * sqrt(cLn2/pi)/γ_d * real(w(modShape,sqrt(cLn2)/γ_d * (grid[I] - ν) + im * y))
end

function testerCPU!(A, grid,γ_d,ν,y,S)
    cSqrtLn2divSqrtPi  = Float32(0.469718639319144059835)
    cLn2               = Float32(0.6931471805599)
    for I in eachindex(grid)
        @inbounds A[I] += S * sqrt(cLn2/pi)/γ_d * real(weideman32a(sqrt(cLn2)/γ_d * (grid[I] - ν) + im * y))
    end
end


end
