using Test, Random
using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.Scattering: Stokes_I, Stokes_IQU
using vSmartMOM.Architectures: architecture
using vSmartMOM.InelasticScattering: noRS

function check_lambertian_jacobian_batch(FT, arr_type; rtol)
    rng = MersenneTwister(582)
    for ns in (1,7), nstokes in (1,3), albedo in (FT(0),FT(0.23)), dark in (false,true), nτ in (0,3)
        n, np = 2*nstokes, 5
        pol = nstokes == 1 ? Stokes_I{FT}() : Stokes_IQU{FT}()
        μ = arr_type(FT[0.3,0.8])
        wt = arr_type(FT[0.5,0.5])
        quad = (; qp_μ=μ, wt_μ=wt, qp_μN=repeat(μ;inner=nstokes),
                  wt_μN=repeat(wt;inner=nstokes), iμ₀Nstart=nstokes+1,
                  iμ₀=2, μ₀=FT(0.8), external_solar=false)
        a, da = CoreRT.make_added_layer(LinMode(), noRS(), FT, arr_type, np, (n,n), ns)
        ar, dar = deepcopy(a), deepcopy(da)
        τ = arr_type(rand(rng,FT,ns))
        dτ = arr_type(randn(rng,FT,ns,nτ))
        F₀ = arr_type(dark ? zeros(FT,nstokes,ns) : rand(rng,FT,nstokes,ns))
        brdf = CoreRT.LambertianSurfaceScalar(albedo)
        old = CoreRT._BATCHED_JACOBIANS_ENABLED[]
        try
            # Follow m=0 with m=1 to check clearing of the source tangents.
            for m in (0,1)
                CoreRT._BATCHED_JACOBIANS_ENABLED[] = false
                CoreRT.create_surface_layer!(noRS(),brdf,ar,dar,4,true,m,pol,quad,
                                             τ,dτ,F₀,architecture(ar.r⁻⁺))
                CoreRT._BATCHED_JACOBIANS_ENABLED[] = true
                CoreRT.create_surface_layer!(noRS(),brdf,a,da,4,true,m,pol,quad,
                                             τ,dτ,F₀,architecture(a.r⁻⁺))
                for field in (:r⁻⁺,:t⁺⁺,:r⁺⁻,:t⁻⁻,:j₀⁺,:j₀⁻)
                    @test Array(getproperty(a,field)) ≈ Array(getproperty(ar,field)) rtol=rtol atol=FT(1e-8)
                end
                for field in (:ap_ṙ⁻⁺,:ap_ṫ⁺⁺,:ap_ṙ⁺⁻,:ap_ṫ⁻⁻,:ap_J̇₀⁺,:ap_J̇₀⁻)
                    @test Array(getproperty(da,field)) ≈ Array(getproperty(dar,field)) rtol=rtol atol=FT(1e-8)
                end
            end
        finally
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = old
        end
    end
end

@testset "Batched Lambertian Jacobians CPU" begin
    check_lambertian_jacobian_batch(Float64,Array;rtol=1e-11)
    check_lambertian_jacobian_batch(Float32,Array;rtol=2e-5)
end
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.functional() || error("CUDA explicitly requested but unavailable")
    CUDA.allowscalar(false)
    @testset "Batched Lambertian Jacobians CUDA" begin
        check_lambertian_jacobian_batch(Float64,CuArray;rtol=1e-10)
        check_lambertian_jacobian_batch(Float32,CuArray;rtol=2e-5)
    end
end
