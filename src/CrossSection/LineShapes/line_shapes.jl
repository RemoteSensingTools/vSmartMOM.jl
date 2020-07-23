using Parameters
using DocStringExtensions
using Interpolations

"""
    type AbstractCrossSection
Abstract Cross Section type for generic cross section calculations
"""
abstract type AbstractCrossSection end
"""
    struct HitranCrossSection{FT}
An [`AbstractCrossSection`](@ref) type struct, which provides all HITRAN line parameters needed to compute absorption cross sections 
see https://hitran.org/docs/definitions-and-units/ for details
# Fields
$(DocStringExtensions.FIELDS)
"""
struct HitranTable{FT<:AbstractFloat} <: AbstractCrossSection
    "The molecular species identification (ID) number"
    mol::Array{Int,1};
    "The isotopologue ID number"
    iso::Array{Int,1};
    "The wavenumber of the spectral line transition (cm-1) in vacuum"    
    νᵢ::Array{FT,1};
    "The spectral line intensity (cm−1/(molecule·cm−2)) at Tref=296K"
    Sᵢ::Array{FT,1};
    "The Einstein-A coefficient (s-1) of a transition"
    Aᵢ::Array{FT,1};
    "The air-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm"
    γ_air::Array{FT,1};
    "The self-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm"
    γ_self::Array{FT,1};
    "The lower-state energy of the transition (cm-1)"
    E″::Array{FT,1};
    "The coefficient of the temperature dependence of the air-broadened half width"
    n_air::Array{FT,1};
    "The pressure shift (cm−1/atm) at Tref=296K and pref=1atm of the line position with respect to the vacuum transition wavenumber νij"
    δ_air::Array{FT,1};
end


function line_shape(
                mod,                   # Line shape model (Voigt here)
                modCEF,                # chosen model for complex error function
                hitran::HitranTable,   # struct  with hitran data  
                grid,                  # wavelength [nm] or wavenumber [cm-1] grid
                wavelength::Bool,      # use wavelength in nm (true) or wavenumber cm-1 units (false)  
                pressure,              # actual pressure [hPa]               
                temperature,           # actual temperature [K]              
                wingCutoff,            # wing cutoff [cm-1]                  
                vmr                    # VMR of gas itself [0-1]
                )
    # Get Type
    FT = eltype(grid);
    CGS_SPEED_OF_LIGHT = FT(2.99792458e10);
    CGS_BOLTZMANN      = FT(1.3806513e-16);
    Nₐ                 = FT(6.02214076e23);
    c₂                 = FT(1.4387769);

    # store results here (or create ! function later)
    result = similar(grid);
    result .= 0.0;

    # need to add partition functions here:...

    p_ref = FT(1013.0);   # reference pressure [hPa]
    t_ref = FT(296.0);    # reference temperature [K]
    wuPi  = FT(sqrt(1/pi));

    # convert to wavenumber from [nm] space if necessary
    if(wavelength) 
        grid = 1.0e7./grid;
    end

    gridMax = maximum(grid)+wingCutoff
    gridMin = minimum(grid)-wingCutoff

    interp_linear_low  = LinearInterpolation(grid, 1:1:length(grid),extrapolation_bc = 1)
    interp_linear_high = LinearInterpolation(grid, 1:1:length(grid),extrapolation_bc = length(grid))
    # Loop through all lines:
    for j in eachindex(hitran.Sᵢ)
        if hitran.νᵢ[j] < gridMax && hitran.νᵢ[j] > gridMin
            # Apply pressure shift
            ν   = hitran.νᵢ[j] + pressure/p_ref*hitran.δ_air[j]
            # Lorentzian HWHM
            γ_l = (hitran.γ_air[j] * 
                  (1-vmr)*pressure/p_ref+hitran.γ_self[j] *
                  vmr*pressure/p_ref) *
                  (t_ref/temperature)^hitran.n_air[j]
            
            # Gaussian HWHM
            #γ_d = hitran.νᵢ[j]/CGS_SPEED_OF_LIGHT*sqrt(2*CGS_BOLTZMANN*temperature*pi/(qmolWeight(MOLEC(j),ISOTOP(j))/Nₐ));
            molWeight = FT(43.98983) # Hardcoded for CO2 now, need qmolWeight(molecID,isoID)
            # Doppler HWHM
            γ_d = hitran.νᵢ[j]/CGS_SPEED_OF_LIGHT *
                  sqrt(2*CGS_BOLTZMANN*temperature*FT(pi)/(molWeight/Nₐ));
            #println(γ_d)
            # Ratio of widths
            y = sqrt(log(FT(2))) * γ_l/γ_d
            # pre factor sqrt(ln(2)/pi)/γ_d
            pre_fac = sqrt(log(FT(2))/pi)/γ_d
            # Line intensity (temperature corrections)
            S = hitran.Sᵢ[j]
            if hitran.E″[j] != -1
                #still needs partition sum correction:
                #S = S * qoft(MOLEC(j),ISOTOP(j),t_ref)/qoft(MOLEC(j),ISOTOP(j),temperature) * 
                S = S * 1 * 
                        exp(c₂*hitran.E″[j]*(1/t_ref-1/temperature)) * 
                        (1-exp(-c₂*hitran.νᵢ[j]/temperature))/(1-exp(-c₂*hitran.νᵢ[j]/t_ref));
            end
            ind_start = Int(round(interp_linear_low(ν-wingCutoff)))
            ind_stop  = Int(round(interp_linear_high(ν+wingCutoff)))
            for i=ind_start:ind_stop
                compl_error = w(modCEF, sqrt(log(FT(2)))/γ_d * (grid[i] - ν) + 1im*y);
                result[i]  += pre_fac * S * real(compl_error)
            end

            #ind = findall(x->x>-wingCutoff && x<wingCutoff, grid .- hitran.νᵢ[j])
            #x = sqrt(log(FT(2)))/γ_d * (grid[ind] .- hitran.νᵢ[j])
            # Complex Error Function w
            #compl_error = w.([modCEF],x .+ 1im*y);
            #plot(real(compl_error))
            
        end
    end
    return result
end
