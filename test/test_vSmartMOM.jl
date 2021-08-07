
# Test the pure Rayleigh case
@testset "pure_rayleigh" begin

    println("Testing pure rayleigh RT...")

    ϵ = 1e-4

    parameters = vSmartMOM.parameters_from_yaml("helper/PureRayleighParameters.yaml")
    model = model_from_parameters(parameters)

    R_model = rt_run(model)

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