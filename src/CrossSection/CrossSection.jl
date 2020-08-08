module CrossSection

using Parameters                # For constructing HitranTable with keywords
using DocStringExtensions       # For simplifying docstring
using Interpolations            # For interpolating in both lookup tables and qoft!

include("constants.jl")         # Scientific and mathematicsl constants
include("types.jl")             # All types used in this module
include("hitran.jl")            # HITRAN file-related functions
include("complex_error_functions.jl")   # Complex error functions used in line broadening
include("mol_weights.jl")       # Molecular weights (TODO: replace with netCDF file)
include("TIPS_2017.jl")         # Partition sums data (TODO: replace with netCDF file)
include("partition_sums.jl")    # Partition sums interpolator (TODO: replace with LinearInterpolation)

# Export the absorption_cross_section functions
export absorption_cross_section

# Export the read_hitran function from hitran.jl
export read_hitran

# Export the broadening function types
export Doppler, Lorentz, Voigt

# Export the complex error functions
export HumlicekErrorFunction, HumlicekWeidemann32VoigtErrorFunction, HumlicekWeidemann32SDErrorFunction, CPF12ErrorFunction, ErfcHumliErrorFunctionVoigt, ErfcHumliErrorFunctionSD, ErfcErrorFunction

# Export the hitran table struct type
export HitranTable, AbstractCrossSection

"""
    $(FUNCTIONNAME)(broadening::AbstractBroadeningFunction, hitran::HitranTable, grid::Array{<:Real,1}, wavelength_flag::Bool, pressure::Real, temperature::Real, wingCutoff::Real; vmr=0)

Read/parse a HITRAN data file and return the data in [`HitranTable`](@ref) format

"""
function absorption_cross_section(
                broadening::AbstractBroadeningFunction,    # Broadening function (AbstractBroadeningFunction)
                hitran::HitranTable,                       # Struct with hitran data
                grid::Array{<:Real,1},                       # Wavelength [nm] or wavenumber [cm-1] grid
                wavelength_flag::Bool,                     # Use wavelength in nm (true) or wavenumber cm-1 units (false)
                pressure::Real,                            # actual pressure [hPa]
                temperature::Real,                         # actual temperature [K]
                wingCutoff::Real;                          # wing cutoff [cm-1]
                vmr=0                                      # VMR of gas itself [0-1]
                )

    # store results here to return
    result = similar(grid);
    result .= 0.0;

    # convert to wavenumber from [nm] space if necessary
    grid = wavelength_flag ? nm_per_m ./ grid : grid

    # Calculate the minimum and maximum grid bounds, including the wing cutoff
    grid_max = maximum(grid) + wingCutoff
    grid_min = minimum(grid) - wingCutoff

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
                qoft!(2,1,temperature,t_ref, rate)
                # println(rate)
                S = S * rate[1] *
                        exp(c₂*hitran.E″[j]*(1/t_ref-1/temperature)) *
                        (1-exp(-c₂*hitran.νᵢ[j]/temperature))/(1-exp(-c₂*hitran.νᵢ[j]/t_ref));

            end

            # Calculate index range that this ν impacts
            ind_start = Int(round(grid_idx_interp_low(ν - wingCutoff)))
            ind_stop  = Int(round(grid_idx_interp_high(ν + wingCutoff)))

            # Loop over every ν in input grid that this transition impacts
            for i=ind_start:ind_stop

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
                    compl_error = w(broadening.CEF, sqrt(cLn2)/γ_d * (grid[i] - ν) + 1im * y);
                    result[i] += pre_fac * S * real(compl_error)
                end
            end
        end
    end

    # Return the resulting lineshape
    return result
end

end
