# using Plots
# using PyCall
# using TimerOutputs
# using Statistics
#
# # const to = TimerOutput();
#
# # ht = readHITRAN("/net/fluo/data1/projects/RTM/Mie/hitran/par/02_hit04.par", 2, 1, 6000, 6400)

# using RadiativeTransfer
# include("/home/rjeyaram/RadiativeTransfer/src/RadiativeTransfer.jl")
using Revise
using RadiativeTransfer.CrossSection

# using Test

ht = CrossSection.readHITRAN("/home/rjeyaram/RadiativeTransfer/test/CO2.data", 2, 1, 6000, 6400)
modCEF = CrossSection.ErfcHumliErrorFunctionVoigt();
grid = collect(6000:0.01:6400);

# @time cs = line_shape(voigt,modCEF,ht,grid,false,1013.0,296,10,0)

# display(plot(grid, cs, ylims=(0, 3e-17)))

# @time cs = line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)

# @time plot(grid, line_shape(doppler,modCEF,ht,grid,false,1013.15,293,40,0), ylims=(0,1e-22))
# @time plot!(grid, line_shape(lorentz,modCEF,ht,grid,false,1013.15,293,40,0), ylims=(0,1e-22))
# @time display(plot!(grid, line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0), ylims=(0,1e-22)))

# @time doppler_shift = line_shape(doppler,modCEF,ht,grid,false,1013.15,293,40,0)
# @time lorentz_shift = line_shape(lorentz,modCEF,ht,grid,false,1013.15,293,40,0)

# voigt_times = []
#
# append!(voigt_times, 1)
#
# for i in 1:100
#     # append!(voigt_times, 1)
#     time = @elapsed res = line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)
#     append!(voigt_times, time)
#     println(time)
# end
#
# println("Mean and stddev")
# println(mean(voigt_times))
# println(std(voigt_times))
#
# # time = @elapsed res = line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)
#
using DelimitedFiles

@time cs = line_shape(voigt,modCEF,ht,grid,false,1013.25,296.0,40,0)

using Plots

maximum(cs)
plot(grid, cs, ylims=(0, 8e-23))
plot!(grid, cs, ylims=(0, 8e-23))

plot

using Interpolations

T = 200:10:380;
p = 1:25:1050;

cs_matrix = zeros(length(T), length(p), length(grid));

for iP in eachindex(p)
    println("iP: ", iP)
    for iT in eachindex(T)
        println("iT: ", iT)
       cs_matrix[iT,iP,:] = line_shape(voigt,modCEF,ht,grid,false,p[iP],T[iT],10.0,0.0);
    end
    # println(iP)
end
itp = interpolate( cs_matrix, BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp,T,p,6000:0.01:6400)


plot(grid, )

using Interact

writedlm( "/home/rjeyaram/RadiativeTransfer/test/voigt_test_julia.csv",  cs, ',')

using Plots

b = @manipulate for curr_T = 200:10:380, curr_P = 1:25:1050

    display(plot(grid, sitp(curr_T, curr_P, grid), ylims=(0, 8e-23)))

end

using JLD2
@save "/home/rjeyaram/RadiativeTransfer/test/savedInterpolation.jld" itp



function app(req)

    grid= -10:1:10

    a = @manipulate for m = 0.1:0.1:2.0#curr_T = 200:10:380, curr_P = 1:25:1050
        plot(grid, m .* grid)
    end

    return a

end

webio_serve(page("/", app), port=5055)

using Mux

ui = button()
WebIO.webio_serve(page("/", req -> b), 5055)

@testset "RadiativeTransfer.jl" begin
    # Write your tests here.




end
