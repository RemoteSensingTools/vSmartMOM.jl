using Revise
using Plots
using Pkg
Pkg.activate(".")
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
# using Distributions
# using BenchmarkTools
# using Test
# using CUDA

# Sets all the "specific" parameters
parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters
model = default_model(parameters);
#model.params.architecture = CPU()
@time R_GPU, T_GPU, R_SFI, T_SFI = vSmartMOM.rt_run(model);
#@show R_GPU[1,1,:]
#@show R_SFI[1,1,:]
#plot(R_GPU[1,1,:])
#plot!(R_SFI[1,1,:])
#plot(R_GPU[1,1,:]./R_SFI[1,1,:])
# curr_parameters.

