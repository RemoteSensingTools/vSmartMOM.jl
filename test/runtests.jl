using RadiativeTransfer
using RadiativeTransfer.Architectures
using RadiativeTransfer.Absorption
using RadiativeTransfer.PhaseFunction
using Test
using DelimitedFiles
using Statistics
using ProgressMeter
using WignerSymbols
using Distributions
using JLD2
using CUDA

# # Determine whether GPU is available
# test_arch = 

# Test the Cross Section module
@testset "RadiativeTransfer.Absorption" begin include("test_Absorption.jl") end

# Test the Phase Function module
@testset "RadiativeTransfer.PhaseFunction" begin include("test_PhaseFunction.jl") end
