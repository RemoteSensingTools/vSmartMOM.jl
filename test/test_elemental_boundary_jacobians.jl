using Test, LinearAlgebra, vSmartMOM
import KernelAbstractions as KA

# Differentiate the independent forward kernels at zero albedo, zero optical
# thickness, and an exactly zero projected phase entry. These are smooth
# elemental limits; dividing a zero forward value by its parameter loses them.
function check_elemental_boundaries(FT,AT;nstokes=3)
    C=vSmartMOM.CoreRT
    n=2nstokes; ns=2; np=4
    μ=AT(repeat(FT[0.3,0.8];inner=nstokes)); μ₀=FT(0.8)
    weights=AT(fill(FT(0.25),n)); D=AT(Diagonal(repeat(FT[1,1,-1,-1][1:nstokes],2)))
    F₀=AT(vcat(FT[1.2 0.7],zeros(FT,nstokes-1,ns))); I₀=AT([one(FT);zeros(FT,nstokes-1)])
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
        # A layer that is doubled (ndoubl ≥ 1) gets the forward D row sign,
        # r⁻⁺ → D·r⁻⁺, on its reflection tangent; t⁺⁺ is unchanged.
        doubled=map(x->AT(zeros(FT,size(x)...)),(primal[1],primal[2],core[1],core[2],
                                                tangent[1],tangent[2],reverse...))
        C.get_elem_rt_fused!(backend)(doubled...,ω,τ,Z,Z,AT(dt),AT(dw),AT(dphase),AT(dphase),
            μ,weights,np,1,nstokes;ndrange=size(primal[1]))
        Dh=Array(D)
        for p in 1:np, s in 1:ns
            @test Array(doubled[5])[:,:,s,p] ≈ Dh*analytic[1][:,:,s,p] rtol=rtol atol=atol
            @test Array(doubled[6])[:,:,s,p] ≈ analytic[2][:,:,s,p] rtol=rtol atol=atol
        end
        for p in 1:np
            plus=forward(τ.+h.*AT(dt[:,p]),ω.+h.*AT(dw[:,p]),above.+h.*AT(da[:,p]),Z.+h.*AT(dphase[:,:,:,p]))
            minus=forward(τ.-h.*AT(dt[:,p]),ω.-h.*AT(dw[:,p]),above.-h.*AT(da[:,p]),Z.-h.*AT(dphase[:,:,:,p]))
            for k in 1:4
                @test analytic[k][:,:,:,p] ≈ (plus[k].-minus[k])./(2h) rtol=rtol atol=atol
            end
            # Undoubled layer (ndoubl = 0): the reverse operators are
            # r⁺⁻ = D·r⁻⁺·D and t⁻⁻ = D·t⁺⁺·D, so their tangents must be too.
            # For IQU, dphase sits in an off-block (Stokes 3 ← 1) entry, where
            # D flips sign; for IQ, D is the identity (Stokes 1 ← 2 entry);
            # for IQUV it is a same-block (Stokes 3 ← 4) entry.
            for s in 1:ns, (rev,fwd) in ((reverse[1],analytic[1]),(reverse[2],analytic[2]))
                @test Array(rev)[:,:,s,p] ≈ Dh*fwd[:,:,s,p]*Dh rtol=rtol atol=atol
            end
        end
    end
end
@testset "Elemental boundary derivatives CPU" begin
    for nstokes in (3,2,4)
        check_elemental_boundaries(Float64,Array;nstokes)
        check_elemental_boundaries(Float32,Array;nstokes)
    end
end
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
    @testset "Elemental boundary derivatives CUDA" begin
        for nstokes in (3,2,4)
            check_elemental_boundaries(Float64,CUDA.CuArray;nstokes)
            check_elemental_boundaries(Float32,CUDA.CuArray;nstokes)
        end
    end
end
