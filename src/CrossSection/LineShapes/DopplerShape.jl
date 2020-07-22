module DopplerShape

# using ..HITRAN
using Pandas

export computeDoppler

function computeDoppler(hitranDF, ν_min, ν_max, ν_step, pressure, temperature, vmr, wingCutoff)

    println(hitranDF)

    for i in 0:length(hitranDF)-1

        r = iloc(hitranDF)[i,:]

        println(r[:nu])
    end

end

end
