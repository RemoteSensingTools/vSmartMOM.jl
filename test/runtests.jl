using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel
using Test
using DelimitedFiles
using Statistics
using ProgressMeter
using WignerSymbols
using Distributions
using JLD2

# Test the Absorption module
#@testset "Absorption" begin include("test_Absorption.jl") end

# Test the Scattering module
#@testset "Scattering" begin include("test_Scattering.jl") end

# Test the CoreRT module
@testset "CoreRT" begin include("test_CoreRT.jl") end

# Test the SolarModel module
#@testset "SolarModel" begin include("test_SolarModel.jl") end
