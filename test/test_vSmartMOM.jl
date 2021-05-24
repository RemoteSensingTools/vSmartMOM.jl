
# Test the pure Rayleigh case
@testset "pure_rayleigh" begin

    println("Testing pure rayleigh RT...")

    parameters = vSmartMOM.parameters_from_yaml("helper/PureRayleighParameters.yaml")

end

# Test the O2 case
@testset "O2 RT" begin

    println("Testing RT with only O₂...")

    parameters = vSmartMOM.parameters_from_yaml("helper/O2Parameters.yaml")

end

# Test 3 bands, 3 profiles case (O₂, CO₂, H₂O)
@testset "pure_rayleigh" begin

    println("Testing 3 bands, 3 profiles case (O₂, CO₂, H₂O)")

    parameters = vSmartMOM.parameters_from_yaml("helper/ThreeBandsParameters.yaml")

end