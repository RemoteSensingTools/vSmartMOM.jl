using Test, Random
using vSmartMOM, vSmartMOM.Scattering
using vSmartMOM.Scattering: GreekCoefs, linGreekCoefs, compute_Z_moments
function check_z_jacobian_tables(FT, arr_type=Array)
    rng = MersenneTwister(529)
    old = Scattering._Z_JACOBIAN_TABLES_ENABLED[]
    try
        for pol in (Stokes_I{FT}(),Stokes_IQ{FT}(),Stokes_IQU{FT}(),Stokes_IQUV{FT}()), ns in (nothing,1,7,64)
            L = 5
            dims = ns === nothing ? (L,) : (L,ns)
            ddims = ns === nothing ? (4,L) : (4,L,ns)
            values = ntuple(_ -> randn(rng,FT,dims),6)
            tangents = ntuple(_ -> randn(rng,FT,ddims),6)
            g, dg = GreekCoefs(values...), linGreekCoefs(tangents...)
            μ = FT[0.2,0.6,0.9]
            for m in (0,1,4,5,7)
                Scattering._Z_JACOBIAN_TABLES_ENABLED[] = false
                reference = compute_Z_moments(pol,μ,g,dg,m)
                Scattering._Z_JACOBIAN_TABLES_ENABLED[] = true
                actual = compute_Z_moments(pol,μ,g,dg,m;arr_type)
                for k in 1:4
                    @test Array(actual[k]) == reference[k]
                end
                if ns === nothing && m == 1
                    # Independent directional derivative, exploiting linearity.
                    h = FT === Float32 ? FT(0.01) : FT(1e-4)
                    p = 2
                    gp = GreekCoefs(map((v,d)->v .+ h .* d[p,:],values,tangents)...)
                    gm = GreekCoefs(map((v,d)->v .- h .* d[p,:],values,tangents)...)
                    plus, minus = compute_Z_moments(pol,μ,gp,m), compute_Z_moments(pol,μ,gm,m)
                    tol = FT === Float32 ? FT(1e-4) : FT(1e-10)
                    for k in 1:2
                        @test (plus[k]-minus[k])/(2h) ≈ reference[k+2][p,:,:] rtol=tol atol=tol
                    end
                end
            end
        end
    finally
        Scattering._Z_JACOBIAN_TABLES_ENABLED[] = old
    end
end
@testset "Tabulated phase Jacobians CPU" begin
    check_z_jacobian_tables(Float64)
    check_z_jacobian_tables(Float32)
end
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.functional() || error("CUDA explicitly requested but unavailable")
    CUDA.allowscalar(false)
    @testset "Tabulated phase Jacobians CUDA output" begin
        check_z_jacobian_tables(Float64,CuArray)
        check_z_jacobian_tables(Float32,CuArray)
    end
end
