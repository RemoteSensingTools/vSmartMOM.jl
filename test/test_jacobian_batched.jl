using Test, Random, LinearAlgebra
using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.Scattering: Stokes_I, Stokes_IQU
using vSmartMOM.Architectures: architecture
using vSmartMOM.InelasticScattering: noRS

function check_jacobian_propagation(FT, arr_type; rtol)
    rng = MersenneTwister(482)
    for ns in (1, 7), nstokes in (1, 3), ndoubl in (0, 3)
        n, np = 2*nstokes, 4
        a, da = CoreRT.make_added_layer(LinMode(), noRS(), FT, arr_type, np, (n,n), ns)
        c, dc = CoreRT.make_composite_layer(LinMode(), noRS(), FT, arr_type, np, (n,n), ns)
        forward_fields = (:r⁻⁺,:t⁺⁺,:r⁺⁻,:t⁻⁻,:j₀⁺,:j₀⁻)
        tangent_fields = (:ap_ṙ⁻⁺,:ap_ṫ⁺⁺,:ap_ṙ⁺⁻,:ap_ṫ⁻⁻,:ap_J̇₀⁺,:ap_J̇₀⁻)
        for field in forward_fields
            x = getproperty(a,field)
            x .= arr_type(rand(rng,FT,size(x)) .* FT(0.03))
        end
        for field in tangent_fields
            x = getproperty(da,field)
            x .= arr_type(randn(rng,FT,size(x)) .* FT(0.001))
            # The last column is inactive in this atmosphere. Keep zero
            # surface slots to exercise the active prefix without changing it.
            x[:,:,:,np] .= zero(FT)
        end
        ar, dar = deepcopy(a), deepcopy(da)
        expk = arr_type(fill(FT(0.9),ns)); er = copy(expk)
        dτ = arr_type(rand(rng,FT,ns,np) .* FT(0.001))
        ident = arr_type(repeat(Matrix{FT}(I,n,n),1,1,ns))
        pol = nstokes == 1 ? Stokes_I{FT}() : Stokes_IQU{FT}()
        old = CoreRT._BATCHED_JACOBIANS_ENABLED[]
        try
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = false
            CoreRT.doubling_allparams!(pol,true,er,ndoubl,ar,dar,ident,
                architecture(ar.r⁻⁺),dτ,FT(0.8);N_active=np-1)
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = true
            @test CoreRT._use_batched_jacobians(da)
            CoreRT.doubling_allparams!(pol,true,expk,ndoubl,a,da,ident,
                architecture(a.r⁻⁺),dτ,FT(0.8);N_active=np-1)
            @test Array(expk) ≈ Array(er) rtol=rtol
            for field in forward_fields
                @test Array(getproperty(a,field)) ≈ Array(getproperty(ar,field)) rtol=rtol atol=FT(1e-8)
            end
            for field in tangent_fields
                @test Array(getproperty(da,field)) ≈ Array(getproperty(dar,field)) rtol=rtol atol=FT(1e-9)
            end

            for field in (:R⁻⁺,:T⁺⁺,:R⁺⁻,:T⁻⁻,:J₀⁺,:J₀⁻)
                x = getproperty(c,field)
                x .= arr_type(rand(rng,FT,size(x)) .* FT(0.04))
            end
            for field in (:Ṙ⁻⁺,:Ṫ⁺⁺,:Ṙ⁺⁻,:Ṫ⁻⁻,:J̇₀⁺,:J̇₀⁻)
                x = getproperty(dc,field)
                x .= arr_type(randn(rng,FT,size(x)) .* FT(0.002))
            end
            cr, dcr = deepcopy(c), deepcopy(dc)
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = false
            CoreRT.interaction!(CoreRT.ScatteringInterface_11(),true,cr,dcr,ar,dar,ident)
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = true
            CoreRT.interaction!(CoreRT.ScatteringInterface_11(),true,c,dc,a,da,ident)
            for field in (:R⁻⁺,:T⁺⁺,:R⁺⁻,:T⁻⁻,:J₀⁺,:J₀⁻)
                @test Array(getproperty(c,field)) ≈ Array(getproperty(cr,field)) rtol=rtol atol=FT(1e-8)
            end
            for field in (:Ṙ⁻⁺,:Ṫ⁺⁺,:Ṙ⁺⁻,:Ṫ⁻⁻,:J̇₀⁺,:J̇₀⁻)
                @test Array(getproperty(dc,field)) ≈ Array(getproperty(dcr,field)) rtol=rtol atol=FT(1e-9)
            end
        finally
            CoreRT._BATCHED_JACOBIANS_ENABLED[] = old
        end
    end
end

@testset "Batched Jacobian propagation CPU" begin
    check_jacobian_propagation(Float64,Array;rtol=1e-11)
    check_jacobian_propagation(Float32,Array;rtol=2e-5)
end

if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.functional() || error("CUDA explicitly requested but unavailable")
    CUDA.allowscalar(false)
    @testset "Batched Jacobian propagation CUDA" begin
        check_jacobian_propagation(Float64,CuArray;rtol=1e-10)
        check_jacobian_propagation(Float32,CuArray;rtol=2e-5)
    end
end
