module VoigtShape

# using ..HITRAN
using Pandas

export computeVoigt

function computeVoigt(hitranDF, ν_min, ν_max, ν_step, pressure, temperature, wingCutoff, vmr)

    p_ref = 1013.0
    t_ref = 296.0

    wuPi = sqrt(1/π)

    # grid = 1.0e7 ./ grid

    # Loop over every transition wavelength
    for ν_curr in collect(ν_min:ν_step:ν_max)

        # r = iloc(hitranDF)[i,:]

        # linePos = r[:nu]+pressure/p_ref*r[:delta_air]

        println(ν_curr)
    end

end

end
