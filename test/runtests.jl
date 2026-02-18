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

# Core module tests
@testset "Absorption" begin include("test_Absorption.jl") end
@testset "Scattering" begin include("test_Scattering.jl") end
@testset "CoreRT" begin include("test_CoreRT.jl") end
@testset "SolarModel" begin include("test_SolarModel.jl") end

# Forward model tests (these require YAML parameter files + data)
@testset "Forward noRS" begin include("test_forward_noRS.jl") end

# Jacobian unit tests
@testset "Jacobian Unit" begin include("test_jacobians_unit.jl") end
