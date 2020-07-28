using RadiativeTransfer.CrossSection
using Test

using Plots
using PyCall
using TimerOutputs
using Statistics

# const to = TimerOutput();

ht = readHITRAN("/net/fluo/data1/projects/RTM/Mie/hitran/par/02_hit04.par", 2, 1, 6000, 6400)

modCEF = HumlicekWeidemann32VoigtErrorFunction();
grid = collect(6000:0.001:6400);

# @time cs = line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)

# @time plot(grid, line_shape(doppler,modCEF,ht,grid,false,1013.15,293,40,0), ylims=(0,1e-22))
# @time plot!(grid, line_shape(lorentz,modCEF,ht,grid,false,1013.15,293,40,0), ylims=(0,1e-22))
# @time display(plot!(grid, line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0), ylims=(0,1e-22)))

# @time doppler_shift = line_shape(doppler,modCEF,ht,grid,false,1013.15,293,40,0)
# @time lorentz_shift = line_shape(lorentz,modCEF,ht,grid,false,1013.15,293,40,0)

voigt_times = []

# append!(voigt_times, 1)

for i in 1:100
    # append!(voigt_times, 1)
    time = @elapsed res = line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)
    append!(voigt_times, time)
    println(time)
end

println("Mean and stddev")
println(mean(voigt_times))
println(std(voigt_times))

# time = @elapsed res = line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)



# tic()

# @timeit to "a" line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)

# println(toc() - tic())

# push!(pyimport("sys")["path"], pwd());
# hapi = pyimport("hapi")
# #
# hapi.absorptionCoefficient_Voigt()
# #
# hapi.fetch("CO2",2,1,6000,6400)
# nu_l, cs_co2_Voigt220 = hapi.absorptionCoefficient_Doppler(SourceTables='CO2',WavenumberStep=0.001, WavenumberRange=[6000,6400],Environment={'p':1.,'T':220},IntensityThreshold=1e-27)

using DelimitedFiles

cs = line_shape(voigt,modCEF,ht,grid,false,1013.15,293,40,0)

writedlm( "/home/rjeyaram/RadiativeTransfer/test/voigt_test_julia.csv",  cs, ',')

@testset "RadiativeTransfer.jl" begin
    # Write your tests here.





end
