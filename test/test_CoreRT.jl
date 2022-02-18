
# Test the pure Rayleigh case
@testset "pure_rayleigh" begin

    println("Testing pure rayleigh RT...")

    ϵ = 1e-4

    parameters = CoreRT.parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
    model = model_from_parameters(parameters)

    R_model = rt_run(model)[1]

    R_test = [1.6376193549450170E-002  2.1561609277437792E-003  5.5387198966105483E-019 0.0 ; 
              1.2818490637667997E-002  3.7825131902305691E-004  2.0685375352086151E-019 0.0 ;
              1.0839231562441865E-002 -3.0496438891765027E-005  4.2614463687432547E-021 0.0 ;
              9.4104223210835589E-003  2.9434545102941607E-004 -1.4472678739285088E-019 0.0 ;
              8.2598796504639114E-003  1.1181829243926966E-003  0.0000000000000000      0.0 ;
              7.3749583918955438E-003  2.3298093802174317E-003  0.0000000000000000      0.0 ;
              6.9126577200309845E-003  3.8960774035191144E-003  0.0000000000000000      0.0 ;
              7.2827593896516067E-003  5.9139825670394474E-003  0.0000000000000000      0.0 ;
              9.6437995452080226E-003  8.8885549319859258E-003  0.0000000000000000      0.0 ]

    @test isapprox(R_model[:,:,1], R_test, atol=ϵ)

end

@testset "compare against 6SV1" begin 

    include("benchmarks/6SV1_R_trues.jl")
    R_modeled_all = zeros(6, 3, 3, 16);
    R_deltas_all = zeros(6, 3, 3, 16);
    parameters = parameters_from_yaml("benchmarks/6SV1_1.yaml");

    function test_against_6SV1(case_i, azs, szas, λ, τ, ρ)
        for sza_i in 1:length(szas)
            for az_i in 1:length(azs)
                R_true = R_trues[case_i]
                R_modeled = view(R_modeled_all, case_i, :, :, :)
                R_deltas = view(R_deltas_all, case_i, :, :, :)
                parameters.spec_bands = [1e7/λ (1e7/λ + 1)]
                parameters.vaz = repeat([azs[az_i]], 16)
                parameters.sza = szas[sza_i]
                parameters.brdf = [vSmartMOM.CoreRT.LambertianSurfaceScalar(ρ * π)]
                model = model_from_parameters(parameters);
                model.τ_rayl[1] .= τ
                R_modeled[sza_i, az_i, :] = CoreRT.rt_run(model, i_band=1)[1][:,1,1] / model.quad_points.μ₀
                R_deltas[sza_i, az_i, :] = abs.(R_true[sza_i][az_i] - R_modeled[sza_i, az_i, :]) ./ R_true[sza_i][az_i]
            end
        end
        return maximum(R_deltas_all[case_i,:,:,:])
    end

    ϵ = 0.006

    @test test_against_6SV1(1, [180, 90, 0], [23.0739, 53.1301, 78.4630], 530, 0.1, 0.0) < ϵ
    @test test_against_6SV1(2, [180, 90, 0], [0.0001, 36.8699, 66.4218], 530, 0.1, 0.25) < ϵ
    @test test_against_6SV1(3, [180, 90, 0], [0.0001, 36.8699, 66.4218], 440, 0.25, 0.0) < ϵ
    @test test_against_6SV1(4, [180, 90, 0], [23.0739, 53.1301, 78.4630], 440, 0.25, 0.25) < ϵ
    @test test_against_6SV1(5, [180, 90, 0], [23.0739, 53.1301, 78.4630], 360, 0.50, 0.0) < ϵ
    @test test_against_6SV1(6, [180, 90, 0], [0.0001, 36.8699, 66.4218], 360, 0.50, 0.25) < ϵ

end

@testset "compare against natraj paper" begin 

    I_modeled_all = zeros(7, 16);
    I_deltas_all  = zeros(7, 16);
    Q_modeled_all = zeros(7, 16);
    Q_deltas_all  = zeros(7, 16);
    U_modeled_all = zeros(7, 16);
    U_deltas_all  = zeros(7, 16);

    μ = [0.02, 0.06, 0.10, 0.16, 0.20, 0.28, 0.32, 0.40, 0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.00]
    ϕs = collect(0.0:30.0:180.0)
    τ = 0.5

    include("benchmarks/natraj_trues.jl")
    parameters = parameters_from_yaml("benchmarks/natraj.yaml");

    parameters.spec_bands = [1e7/360.0 (1e7/360.0 + 1)]
    parameters.vza = acosd.(μ)
    parameters.sza = acosd.(0.2)

    for ϕ_i in 1:length(ϕs)
        parameters.vaz = repeat([ϕs[ϕ_i]], 16)
        model = model_from_parameters(parameters);
        model.τ_rayl[1] .= τ

        R = CoreRT.rt_run(model, i_band=1)[1]
        @show size(R)

        I_modeled_all[ϕ_i,:] = R[:,1,1]
        Q_modeled_all[ϕ_i,:] = R[:,2,1]
        U_modeled_all[ϕ_i,:] = R[:,3,1]

        I_deltas_all[ϕ_i,:] = abs.(I_trues[:,ϕ_i] - I_modeled_all[ϕ_i,:]) ./ I_trues[:,ϕ_i]
        Q_deltas_all[ϕ_i,:] = abs.(Q_trues[:,ϕ_i] - Q_modeled_all[ϕ_i,:]) ./ Q_trues[:,ϕ_i]
        U_deltas_all[ϕ_i,:] = abs.(U_trues[:,ϕ_i] - U_modeled_all[ϕ_i,:]) ./ U_trues[:,ϕ_i]
    end

    ϵ = 0.008

    @test maximum(I_deltas_all) < 0.002
    @test maximum(Q_deltas_all[findall(i -> i >= 0.01, Q_modeled_all)]) < ϵ
    @test maximum(filter(!isnan, U_deltas_all[findall(i -> i >= 0.01, U_modeled_all)])) < ϵ

end

# # Test the O2 case
# @testset "O2 RT" begin

#     println("Testing RT with only O₂...")

#     parameters = vSmartMOM.parameters_from_yaml("helper/O2Parameters.yaml")

# end

# # Test 3 bands, 3 profiles case (O₂, CO₂, H₂O)
# @testset "pure_rayleigh" begin

#     println("Testing 3 bands, 3 profiles case (O₂, CO₂, H₂O)")

#     parameters = vSmartMOM.parameters_from_yaml("helper/ThreeBandsParameters.yaml")

# end