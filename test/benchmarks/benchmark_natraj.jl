##

using RadiativeTransfer, RadiativeTransfer.vSmartMOM

## 

I_modeled_all = zeros(7, 16);
I_deltas_all = zeros(7, 16);

Q_modeled_all = zeros(7, 16);
Q_deltas_all = zeros(7, 16);

U_modeled_all = zeros(7, 16);
U_deltas_all = zeros(7, 16);

## 

μ = [0.02, 0.06, 0.10, 0.16, 0.20, 0.28, 0.32, 0.40, 0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.00]
ϕs = collect(0.0:30.0:180.0)

acosd.(μ)

include("test/benchmarks/natraj_trues.jl")
parameters = parameters_from_yaml("test/benchmarks/natraj.yaml");

λ = 360.0
parameters.spec_bands = [1e7/λ (1e7/λ + 1)]
parameters.vza = acosd.(μ)
parameters.sza = acosd.(0.2)

τ = 0.5

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

@show maximum(I_deltas_all)
println(maximum(Q_deltas_all[findall(i -> i >= 0.01, Q_modeled_all)]))
@show maximum(filter(!isnan, U_deltas_all[findall(i -> i >= 0.01, U_modeled_all)])) 