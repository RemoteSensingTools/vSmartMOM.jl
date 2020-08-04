module CrossSection

using Parameters
using DocStringExtensions
using Interpolations

include("types.jl")
include("hitran.jl")
include("complex_error_functions.jl")
include("iso_properties.jl")

include("TIPS_2017.jl")
include("partition_sums.jl")

export readHITRAN, line_shape, doppler, lorentz, voigt, HumlicekErrorFunction, HumlicekWeidemann32VoigtErrorFunction, HumlicekWeidemann32SDErrorFunction, CPF12ErrorFunction, ErfcHumliErrorFunctionVoigt, ErfcHumliErrorFunctionSD, ErfcErrorFunction, qoft

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

    cMassMol = 1.66053873e-27
    cSqrtLn2divSqrtPi = 0.469718639319144059835
    cLn2 = 0.6931471805599
    cSqrtLn2 = 0.8325546111577
    cSqrt2Ln2 = 1.1774100225
    fSqrtMass = sqrt(43.98983)
    cc_ = 2.99792458e8
    cBolts_ = 1.3806503e-23
    cSqrt2Ln2 = 1.1774100225

    # store results here (or create ! function later)
    # result = SharedArray{Float64}(length(grid))
    result = similar(grid);
    result .= 0.0;

    # need to add partition functions here:...

    p_ref = FT(1013.25);  # reference pressure [hPa]
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

    times = []

    # rate=0.0
    # rate = Vector

    # println("runnning latest")

    rate = zeros(1)
    # println(rate)

    # Loop through all lines:
    for j in eachindex(hitran.Sᵢ)
        # println(j)
        if hitran.νᵢ[j] < gridMax && hitran.νᵢ[j] > gridMin
            # Apply pressure shift
            ν   = hitran.νᵢ[j] + pressure/p_ref*hitran.δ_air[j]
            # Lorentzian HWHM
            γ_l = (hitran.γ_air[j] *
                  (1-vmr)*pressure/p_ref+hitran.γ_self[j] *
                  vmr*pressure/p_ref) *
                  (t_ref/temperature)^hitran.n_air[j]

            # Gaussian HWHM


            # Doppler HWHM
            molWeight = FT(mol_weight(hitran.mol[j],hitran.iso[j]))
            γ_d = (cSqrt2Ln2/cc_)*sqrt(cBolts_/cMassMol)*sqrt(temperature) * hitran.νᵢ[j]/fSqrtMass

            # Ratio of widths
            y = sqrt(log(FT(2))) * γ_l/γ_d

            # pre factor sqrt(ln(2)/pi)/γ_d
            pre_fac = sqrt(log(FT(2))/pi)/γ_d

            # Line intensity (temperature corrections)
            S = hitran.Sᵢ[j]
            if hitran.E″[j] != -1
                #still needs partition sum correction:

                # S = S * 1 *
                # CURRT = typeof()

                qoft!(2,1,temperature,t_ref, rate)


                # append!(times, time)
                # println(rate)
                # println(typeof(S))
                # println(typeof(rate))
                # println(typeof(c₂))
                # println(typeof(hitran.E″[j]))
                # println(typeof(1/t_ref-1/temperature))
                # println(typeof((1-exp(-c₂*hitran.νᵢ[j]/temperature))))
                # println(typeof(1-exp(-c₂*hitran.νᵢ[j]/t_ref)))
                #
                # println((S))
                # println((rate))
                # println((c₂))
                # println((hitran.E″[j]))
                # println((1/t_ref-1/temperature))
                # println(((1-exp(-c₂*hitran.νᵢ[j]/temperature))))
                # println((1-exp(-c₂*hitran.νᵢ[j]/t_ref)))

                S = S * rate[1] *
                        exp(c₂*hitran.E″[j]*(1/t_ref-1/temperature)) *
                        (1-exp(-c₂*hitran.νᵢ[j]/temperature))/(1-exp(-c₂*hitran.νᵢ[j]/t_ref));

                # break

            end
            ind_start = Int(round(interp_linear_low(ν-wingCutoff)))
            ind_stop  = Int(round(interp_linear_high(ν+wingCutoff)))
            for i=ind_start:ind_stop

                if mod === doppler
                    lineshape_val = cSqrtLn2divSqrtPi*exp(-cLn2*((grid[i] - ν) /γ_d) ^2) /γ_d
                    result[i] += S * lineshape_val

                elseif mod === lorentz

                    lineshape_val = γ_l/(pi*(γ_l^2+(grid[i] - ν) ^ 2))
                    result[i] += S * lineshape_val

                elseif mod === voigt
                    compl_error = w(modCEF, sqrt(log(FT(2)))/γ_d * (grid[i] - ν) + 1im*y);
                    result[i] += pre_fac * S * real(compl_error)
                end


            end

            #ind = findall(x->x>-wingCutoff && x<wingCutoff, grid .- hitran.νᵢ[j])
            #x = sqrt(log(FT(2)))/γ_d * (grid[ind] .- hitran.νᵢ[j])
            # Complex Error Function w
            #compl_error = w.([modCEF],x .+ 1im*y);
            #plot(real(compl_error))

        end
    end

    # println(sum(times))

    return result # (result, times)
end

end
