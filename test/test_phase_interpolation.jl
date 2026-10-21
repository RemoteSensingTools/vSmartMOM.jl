using Test, Random, vSmartMOM
using vSmartMOM.CoreRT

function check_phase_interpolation(AT)
    rng = Xoshiro(319)
    for FT in (Float32,Float64), dims in ((6,6),(4,6,6)), nn in (2,3), reverse_nodes in (false,true)
        knots = FT.([1.0,1.4,2.0][nn == 2 ? [1,3] : [1,2,3]])
        values = [randn(rng,FT,dims...) for _ in 1:nn]
        if reverse_nodes
            reverse!(knots); reverse!(values)
        end
        # Include both endpoints, an interior knot, and extrapolation. The
        # interpolation policy must agree with the established host routine.
        grid = FT.(range(0.9,2.1;length=17))
        old = CoreRT._interpolate_phase_nodes(grid,knots,values)
        new = Array(CoreRT.interpolate_phase_blocks(grid,knots,values,AT))
        @test size(new) == (dims...,length(grid))
        @test new ≈ old rtol=20eps(FT) atol=20eps(FT)
        at_knots = Array(CoreRT.interpolate_phase_blocks(knots,knots,values,AT))
        for i in 1:nn
            @test selectdim(at_knots,ndims(at_knots),i) ≈ values[i] rtol=10eps(FT) atol=10eps(FT)
        end
        if FT === Float64
            # Independent central difference of the original host evaluator
            # verifies interpolation of tangent blocks at fixed knots.
            tangent = [randn(rng,FT,dims...) for _ in 1:nn]
            h = 1e-5
            plus = CoreRT._interpolate_phase_nodes(grid,knots,values .+ h.*tangent)
            minus = CoreRT._interpolate_phase_nodes(grid,knots,values .- h.*tangent)
            derivative = Array(CoreRT.interpolate_phase_blocks(grid,knots,tangent,AT))
            @test derivative ≈ (plus-minus)/(2h) rtol=1e-9 atol=1e-9
        end
    end
    @test_throws ArgumentError CoreRT.interpolate_phase_blocks([1.0],[1.0,1.0],[ones(2,2),ones(2,2)],AT)
end

@testset "Phase interpolation CPU" begin
    check_phase_interpolation(Array)
end
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
    @testset "Phase interpolation CUDA" begin
        check_phase_interpolation(CUDA.CuArray)
    end
end
