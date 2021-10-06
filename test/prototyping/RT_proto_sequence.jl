using Revise
using Plots, Plots.PlotMeasures
using Pkg
# Pkg.activate(".")
using BenchmarkTools
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using RadiativeTransfer.SolarModel

I_modeled_all = zeros(7, 16);
I_deltas_all  = zeros(7, 16);
Q_modeled_all = zeros(7, 16);
Q_deltas_all  = zeros(7, 16);
U_modeled_all = zeros(7, 16);
U_deltas_all  = zeros(7, 16);

μ = [0.02, 0.06, 0.10, 0.16, 0.20, 0.28, 0.32, 0.40, 0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.00]
ϕs = collect(0.0:30.0:180.0)
τ = 0.5

include("/home/rjeyaram/RadiativeTransfer/test/benchmarks/natraj_trues.jl")
parameters = parameters_from_yaml("/home/rjeyaram/RadiativeTransfer/test/benchmarks/natraj.yaml");

parameters.spec_bands = [1e7/360.0 (1e7/360.0 + 1)]
parameters.vza = acosd.(μ)
parameters.sza = acosd.(0.2)

for ϕ_i in 1:length(ϕs)
    parameters.vaz = repeat([ϕs[ϕ_i]], 16)
    model = model_from_parameters(parameters);
    model.τ_rayl[1] .= τ

    R = vSmartMOM.rt_run(model, i_band=1)
    @show size(R)

    I_modeled_all[ϕ_i,:] = R[:,1,1]
    Q_modeled_all[ϕ_i,:] = R[:,2,1]
    U_modeled_all[ϕ_i,:] = R[:,3,1]

    I_deltas_all[ϕ_i,:] = abs.(I_trues[:,ϕ_i] - I_modeled_all[ϕ_i,:]) ./ I_trues[:,ϕ_i]
    Q_deltas_all[ϕ_i,:] = abs.(Q_trues[:,ϕ_i] - Q_modeled_all[ϕ_i,:]) ./ Q_trues[:,ϕ_i]
    U_deltas_all[ϕ_i,:] = abs.(U_trues[:,ϕ_i] - U_modeled_all[ϕ_i,:]) ./ U_trues[:,ϕ_i]
end


plot(parameters.vza, I_modeled_all[1,:])
plot!(parameters.vza, I_trues[:,1])

## 

plot_phis = [1, 4, 7]
plot_list = []
gr(display_type=:inline)
for i in plot_phis
    p = plot(parameters.vza, Q_modeled_all[i,:], xlabel = "VZA", title = "Q, VAZ = $(ϕs[i])°", legend=:bottomleft, label="model", left_margin = 10mm)
    p = plot!(parameters.vza, Q_trues[:,i], label="natraj")
    push!(plot_list, p)

    p = plot(parameters.vza, Q_modeled_all[i,:] - Q_trues[:,i], xlabel = "VZA", title = "Model – Natraj", legend=:bottomright, label="diff", ylims=(-0.0005, 0.0005))

    push!(plot_list, p)

    p = plot(parameters.vza, Q_deltas_all[i,:] * 100, xlabel = "VZA", title = "% Difference", legend=false, label="diff")

    push!(plot_list, p)
end

plot(plot_list..., layout=(3,3), size=(1200,1200))
