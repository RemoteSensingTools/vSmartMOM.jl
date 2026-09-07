using Test, Random, vSmartMOM
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
    function check_vendor_jacobians(FT,n)
        C=vSmartMOM.CoreRT
        rng=Xoshiro(571)
        ns,np=512,4
        A,B=map(_->CuArray(randn(rng,FT,n,n,ns)),1:2)
        dA,dB=map(_->CuArray(randn(rng,FT,n,n,ns,np)),1:2)
        out,reference=zero(dA),zero(dA)
        w=C.make_jacobian_workspace(A,np)
        @test w.products !== nothing
        entries=0
        for (iteration,active) in enumerate((np,2,np,2))
            inactive=Array(out[:,:,:,active+1:np])
            # Fresh active-prefix views must reuse their pointer metadata.
            a,b=C._jac_active(dA,active),C._jac_active(dB,active)
            dst,ref=C._jac_active(out,active),C._jac_active(reference,active)
            C._jprod!(dst,A,a,B,b)
            expected=Array(dst)
            C._jprod!(w,dst,A,a,B,b)
            @test Array(dst) ≈ expected rtol=200eps(FT) atol=200eps(FT)
            for (left,right) in ((A,b),(a,B),(a,b))
                fill!(dst,FT(0.7));fill!(ref,FT(0.7))
                C._jmul!(ref,left,right,FT(0.3))
                C._jmul!(w,dst,left,right,FT(0.3))
                @test Array(dst) ≈ Array(ref) rtol=200eps(FT) atol=200eps(FT)
            end
            @test isequal(Array(out[:,:,:,active+1:np]),inactive)
            if iteration==2
                entries=length(w.products.pointers)
            elseif iteration>2
                @test length(w.products.pointers)==entries
            end
        end
        # Pointer vectors must never be reused across deep-copied workspaces.
        copy=deepcopy(w.products)
        @test isempty(copy.pointers)
        @test length(w.products.pointers)==entries
        src,dsrc=CuArray(randn(rng,FT,n,1,ns)),CuArray(randn(rng,FT,n,1,ns,np))
        r1,r2=similar(dsrc),similar(dsrc)
        C._jprod!(r1,A,dA,src,dsrc)
        C._jprod!(w,r2,A,dA,src,dsrc)
        @test Array(r1) ≈ Array(r2) rtol=200eps(FT) atol=200eps(FT)
        @test length(w.products.pointers)==entries # source vectors bypass cuBLAS plans
    end
    @testset "CUDA Jacobian pointer plans" begin
        C=vSmartMOM.CoreRT
        old=C._VENDOR_JACOBIANS_ENABLED[]
        try
            C._VENDOR_JACOBIANS_ENABLED[]=true
            for FT in (Float64,Float32), n in (33,57)
                check_vendor_jacobians(FT,n)
            end
        finally
            C._VENDOR_JACOBIANS_ENABLED[]=old
        end
    end
end
