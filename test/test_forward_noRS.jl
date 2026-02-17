# =================================================================
# Forward Model Smoke Test (no Raman scattering)
# Tests basic forward RT with an EMIT-style setup
# =================================================================
using vSmartMOM, vSmartMOM.CoreRT
using Test

println("="^60)
println("Forward Model Smoke Test (noRS)")
println("="^60)

# Load EMIT-style fast configuration
params = parameters_from_yaml("test/test_parameters/ParamsEMIT_fast.yaml")
Nbands = length(params.spec_bands)
println("Configuration loaded: $(Nbands) band(s)")

# Build model (forward-only, no linearization)
println("\nBuilding forward model...")
@time model = model_from_parameters(params)
println("Model built successfully.")
println("  Layers: $(length(model.profile.T))")
println("  Quadrature points: $(model.quad_points.Nquad)")

@testset "Forward noRS Smoke Test" begin

    @testset "Model construction" begin
        @test !isnothing(model)
        @test length(model.params.spec_bands) == Nbands
        @test all(isfinite.(model.profile.T))
        @test all(model.profile.p_full .> 0)
        
        # Optical properties should be finite
        for ib in 1:Nbands
            @test all(isfinite.(model.τ_abs[ib]))
            @test all(model.τ_abs[ib] .>= 0)
            @test all(isfinite.(model.τ_rayl[ib]))
            @test all(model.τ_rayl[ib] .>= 0)
        end
    end

    @testset "Band $ib RT run" for ib in 1:Nbands
        println("\n  Running RT for band $ib...")
        @time result = rt_run(model; i_band=ib)
        # rt_run returns a tuple: (R, T, ieR, ieT, ...) — extract reflectance
        R = result isa Tuple ? result[1] : result

        nStokes = model.params.polarization_type.n
        nVza = length(model.obs_geom.vza)
        nSpec = length(model.params.spec_bands[ib])
        println("  Output shape: ($(size(R)))")
        println("  nStokes=$nStokes, nVza=$nVza, nSpec=$nSpec")

        @testset "Output dimensions" begin
            @test ndims(R) == 3
            @test size(R, 1) == nVza
            @test size(R, 2) == nStokes
            @test size(R, 3) == nSpec
        end

        @testset "Physical sanity" begin
            # All values should be finite
            @test all(isfinite.(R))
            
            # Stokes I (radiance) should be positive
            I_vals = R[:, 1, :]
            @test all(I_vals .> 0) || @warn "Some I values non-positive: min=$(minimum(I_vals))"
            
            # I should be bounded (not unreasonably large)
            @test maximum(I_vals) < 1.0  # Reflectance < 1 for realistic scenes
            
            # |Q| and |U| should be less than I (degree of polarization <= 1)
            if nStokes >= 2
                Q_vals = R[:, 2, :]
                @test all(abs.(Q_vals) .<= I_vals .+ 1e-10)
            end
            if nStokes >= 3
                U_vals = R[:, 3, :]
                @test all(abs.(U_vals) .<= I_vals .+ 1e-10)
            end
        end

        @testset "Multi-VZA consistency" begin
            if nVza > 1
                # Different viewing angles should give different results
                @test !all(R[1, 1, :] .≈ R[2, 1, :])
            end
        end

        # Print summary statistics
        I_vals = R[:, 1, :]
        println("  I: min=$(minimum(I_vals)), max=$(maximum(I_vals)), mean=$(sum(I_vals)/length(I_vals))")
        if nStokes >= 2
            Q_vals = R[:, 2, :]
            println("  Q: min=$(minimum(Q_vals)), max=$(maximum(Q_vals))")
        end
    end
end

println("\n" * "="^60)
println("Forward noRS smoke test complete.")
println("="^60)
