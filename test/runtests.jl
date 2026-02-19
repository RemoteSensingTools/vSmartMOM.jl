using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel
using Test
using DelimitedFiles
using Statistics
using LinearAlgebra
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

# Type stability tests (no external data)
@testset "Type Stability" begin include("test_type_stability.jl") end

# Float32 consistency tests
@testset "Float32 Consistency" begin include("test_float32.jl") end

# RAMI smoke test (no gas absorption, bypasses sandbox scripts)
@testset "RAMI Smoke" begin include("test_rami_smoke.jl") end

# Hybrid AD tests (ForwardDiff Mie + linearized RT)
@testset "Hybrid AD" begin include("test_hybrid_ad.jl") end

# Canopy surface tests (CanopySurface composable surface type)
@testset "Canopy Surface" begin include("test_canopy.jl") end

# GPU-specific tests (conditional on CUDA availability)
CUDA_AVAILABLE = try
    using CUDA
    CUDA.functional()
catch
    false
end

if CUDA_AVAILABLE
    @testset "Raman GPU" begin include("test_forward_raman_gpu.jl") end
end
