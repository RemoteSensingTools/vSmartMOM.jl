using RadiativeTransfer
using RadiativeTransfer.Architectures
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using RadiativeTransfer.SolarModel
using Test
using DelimitedFiles
using Statistics
using ProgressMeter
using WignerSymbols
using Distributions
using JLD2

# Test the Absorption module
# @testset "RadiativeTransfer.Absorption" begin include("test_Absorption.jl") end

# Test the Scattering module
# @testset "RadiativeTransfer.Scattering" begin include("test_Scattering.jl") end

# Test the vSmartMOM module
@testset "RadiativeTransfer.vSmartMOM" begin include("test_vSmartMOM.jl") end

# Test the SolarModel module
# @testset "RadiativeTransfer.SolarModel" begin include("test_SolarModel.jl") end
