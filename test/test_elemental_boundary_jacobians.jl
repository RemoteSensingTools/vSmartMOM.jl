using Test, LinearAlgebra, vSmartMOM
import KernelAbstractions as KA

# Differentiate the independent forward kernels at zero albedo, zero optical
# thickness, and an exactly zero projected phase entry. These are smooth
# elemental limits; dividing a zero forward value by its parameter loses them.
function check_elemental_boundaries(FT,AT)
    C=vSmartMOM.CoreRT
    nstokes=3; n=6; ns=2; np=4
    μ=AT(FT[0.3,0.3,0.3,0.8,0.8,0.8]); μ₀=FT(0.8)
    weights=AT(fill(FT(0.25),n)); D=AT(Diagonal(FT[1,1,-1,1,1,-1]))
    F₀=AT(FT[1.2 0.7; 0 0; 0 0]); I₀=AT(FT[1,0,0])
    phase=fill(FT(0.5),n,n,1); phase[3,4,1]=0
    dphase=zeros(FT,n,n,1,np); dphase[3,4,1,3]=1
    dt=zeros(FT,ns,np); dt[:,1].=1
    dw=zeros(FT,ns,np); dw[:,2].=1
    da=zeros(FT,ns,np); da[:,4].=1
    make_arrays()=(AT(zeros(FT,n,n,ns)),AT(zeros(FT,n,n,ns)),
                   AT(zeros(FT,n,1,ns)),AT(zeros(FT,n,1,ns)))
    primal=make_arrays()
    core=map(x->AT(zeros(FT,size(x)...,3)),primal)
    tangent=map(x->AT(zeros(FT,size(x)...,np)),primal)
    reverse=ntuple(_->AT(zeros(FT,n,n,ns,np)),2)
    backend=KA.get_backend(first(primal))
    function forward(τ,ω,above,Z)
        out=make_arrays()
        C.get_elem_rt!(backend)(out[1],out[2],ω,τ,Z,Z,μ,weights;ndrange=size(out[1]))
        C.get_elem_rt_SFI!(backend)(out[3],out[4],ω,τ,above,Z,Z,F₀,μ,μ₀,
            0,FT(0.5),nstokes,2,D;ndrange=size(out[3]))
        map(Array,out)
    end
    h=FT===Float64 ? FT(1e-6) : FT(0.002)
    rtol=FT===Float64 ? 2e-7 : 6e-3
    atol=FT===Float64 ? 1e-9 : 2e-5
    for τvalue in (zero(FT),FT(0.02)), ωvalue in (zero(FT),FT(0.8))
        τ=AT(fill(τvalue,ns)); ω=AT(fill(ωvalue,ns)); above=AT(fill(FT(0.1),ns)); Z=AT(phase)
        foreach(x->fill!(x,0),primal);foreach(x->fill!(x,0),core);foreach(x->fill!(x,0),tangent)
        C.get_elem_rt_fused!(backend)(primal[1],primal[2],core[1],core[2],
            tangent[1],tangent[2],reverse...,ω,τ,Z,Z,AT(dt),AT(dw),AT(dphase),AT(dphase),
            μ,weights,np,0,nstokes;ndrange=size(primal[1]))
        C.get_elem_rt_SFI_fused!(backend)(primal[3],primal[4],core[3],core[4],
            tangent[3],tangent[4],ω,τ,above,AT(da),Z,Z,F₀,AT(dt),AT(dw),AT(dphase),AT(dphase),
            μ,0,FT(0.5),nstokes,I₀,μ₀,2,D,np;ndrange=size(primal[3]))
        analytic=map(Array,tangent)
        for (a,b) in zip(map(Array,primal),forward(τ,ω,above,Z))
            @test a ≈ b rtol=rtol atol=atol
        end
        for p in 1:np
            plus=forward(τ.+h.*AT(dt[:,p]),ω.+h.*AT(dw[:,p]),above.+h.*AT(da[:,p]),Z.+h.*AT(dphase[:,:,:,p]))
            minus=forward(τ.-h.*AT(dt[:,p]),ω.-h.*AT(dw[:,p]),above.-h.*AT(da[:,p]),Z.-h.*AT(dphase[:,:,:,p]))
            for k in 1:4
                @test analytic[k][:,:,:,p] ≈ (plus[k].-minus[k])./(2h) rtol=rtol atol=atol
            end
        end
    end
end
@testset "Elemental boundary derivatives CPU" begin
    check_elemental_boundaries(Float64,Array)
    check_elemental_boundaries(Float32,Array)
end
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
    @testset "Elemental boundary derivatives CUDA" begin
        check_elemental_boundaries(Float64,CUDA.CuArray)
        check_elemental_boundaries(Float32,CUDA.CuArray)
    end
end
