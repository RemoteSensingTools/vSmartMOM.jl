#!/usr/bin/env julia
# =============================================================================
# raman_fused_kernels.jl — VALIDATED prototype fused KernelAbstractions kernels
# replacing the per-Δn batched_mul host loops in the RRS adding-doubling solver
# (doubling_inelastic.jl / interaction_inelastic.jl). NOT yet wired into the
# package — this is the reviewed design template for integration.
#
#   CUDA_VISIBLE_DEVICES=<free> julia --project=test test/benchmarks/raman_fused_kernels.jl
#
# Motivation: the per-Δn loops launch ~13 tiny (NquadN×NquadN, NquadN=15)
# batched_mul + broadcasts per iteration, ~16% GPU-busy (launch-bound) with HBM
# round-trips for every intermediate. Fusing the whole per-(n₁,Δn) chain into one
# custom kernel (one workgroup per (n₁,Δn); 15×15 threads = one per output
# element; intermediates kept in @localmem; CSEs computed once; Raman shift
# n₀ = n₁ + i_λ₁λ₀[Δn] resolved in-kernel) recovers ~100% GPU utilization.
#
# MEASURED (L40S, f32, NquadN=15, nRaman=200):
#   doubling n-loop-2: batched_mul 234.5 ms → fused 38.4 ms (6.1×) @ nSpec=10000
#                      (22.5× @1000, 5.7× @5000; fused 99.8% GPU-busy, scales linearly)
#   interaction loop-1 tmpieT⁻⁻ (single-output kernel): correct & stable
#
# ACCURACY (vs f64 reference of the SAME inputs; curry A100 ≡ wurst L40S bit-for-bit):
#   doubling fused:  f64 2.1e-16,  f32 1.1e-7     (Neumaier f32 6e-8, +13% time)
#   interaction oT:  f64 4.65e-16, f32 1.84e-7
#   => correct to machine precision; NOT bit-identical to batched_mul (regen goldens).
#
# ⚠ BUFFER-COUNT CONSTRAINT (found the hard way): a monolithic kernel with ~18
#   @localmem (15×15) arrays CORRUPTS outputs (flaky NaN/race), even though total
#   shared (~32 KB f64) is under the 48 KB limit — likely a GPUCompiler/register
#   issue. Keep each kernel ≤ ~14 @localmem buffers. => fuse interaction as
#   PER-OUTPUT kernels (6 total), recomputing the W/Wt CSE per kernel.
#   ALWAYS validate each kernel f64-vs-batched_mul AND run ≥3× (catch races).
#
# ⚠ REFERENCE PITFALL: when writing a batched_mul reference, COPY the old slices
#   (`iet=copy(IET[:,:,n1,Δn])`) — writing IET in-place before computing IER makes
#   a view alias updated data (this masked as a 7.7% "error" until the f64 check).
#
# ⚠ KA gotchas: @localmem type must be inlined (`@localmem eltype(X) (NQ,NQ)`, not
#   a local `E=eltype`); @localmem must be at kernel top; `return` is forbidden
#   (wrap body in `if`); KA CPU backend can't run these (scalar locals don't
#   survive @synchronize splits) — use the batched_mul reference / curry f64 GPU.
# =============================================================================
ENV["CUDA_VISIBLE_DEVICES"] = get(ENV,"CUDA_VISIBLE_DEVICES","0")
using CUDA, NNlib, KernelAbstractions, LinearAlgebra, Random, Statistics, Printf
const KA = KernelAbstractions; const ⊠ = NNlib.batched_mul; const NQ = 15

# row-i · col-j dot of two NQ×NQ shared matrices; comp=true → Neumaier compensated sum
@inline function dotc(X, Y, i, j, comp)
    Tt = eltype(X); s = zero(Tt); c = zero(Tt)
    @inbounds for k in 1:NQ
        p = X[i,k]*Y[k,j]
        comp ? (t=s+p; c += ifelse(abs(s)>=abs(p),(s-t)+p,(p-t)+s); s=t) : (s += p)
    end
    return comp ? s+c : s
end

get_nn(nSpec,Δ) = (a=max(1,1-Δ); b=min(nSpec,nSpec-Δ); (a:b,(a+Δ):(b+Δ)))

# ─── DOUBLING n-loop-2 (replaces doubling_inelastic.jl:120-143) ───────────────
# Updates ier⁻⁺/iet⁺⁺ in place. Both tmp5(iet)/tmp6(ier) use OLD ier/iet (read
# once into shared), so in-place write is race-free per (n₁,Δn) block.
@kernel function dbl_nloop2!(IET, IER, @Const(TTGP), @Const(GP), @Const(R), @Const(T),
                             @Const(shifts), nSpec, comp::Bool)
    _i,_j, n1, Δn = @index(Global, NTuple); li, lj, _, _ = @index(Local, NTuple)
    sTTGP=@localmem eltype(IET) (NQ,NQ); sGP0=@localmem eltype(IET) (NQ,NQ); sR0=@localmem eltype(IET) (NQ,NQ)
    sR1=@localmem eltype(IET) (NQ,NQ); sT0=@localmem eltype(IET) (NQ,NQ); sIET=@localmem eltype(IET) (NQ,NQ); sIER=@localmem eltype(IET) (NQ,NQ)
    sQ=@localmem eltype(IET) (NQ,NQ); sGPt0=@localmem eltype(IET) (NQ,NQ); sGPr0=@localmem eltype(IET) (NQ,NQ)
    s1=@localmem eltype(IET) (NQ,NQ); s2=@localmem eltype(IET) (NQ,NQ); s3=@localmem eltype(IET) (NQ,NQ); sT5=@localmem eltype(IET) (NQ,NQ)
    Δ=@inbounds shifts[Δn]; n0=n1+Δ
    if 1<=n0<=nSpec
        @inbounds begin
            sTTGP[li,lj]=TTGP[li,lj,n1]; sGP0[li,lj]=GP[li,lj,n0]; sR0[li,lj]=R[li,lj,n0]
            sR1[li,lj]=R[li,lj,n1]; sT0[li,lj]=T[li,lj,n0]; sIET[li,lj]=IET[li,lj,n1,Δn]; sIER[li,lj]=IER[li,lj,n1,Δn]
        end
        @synchronize
        @inbounds sQ[li,lj]=dotc(sIER,sR0,li,lj,comp)+dotc(sR1,sIER,li,lj,comp)
        @inbounds sGPt0[li,lj]=dotc(sGP0,sT0,li,lj,comp); @inbounds sGPr0[li,lj]=dotc(sGP0,sR0,li,lj,comp)
        @synchronize
        @inbounds s1[li,lj]=sIET[li,lj]+dotc(sQ,sGPt0,li,lj,comp)
        @synchronize
        @inbounds sT5[li,lj]=dotc(sTTGP,s1,li,lj,comp)+dotc(sIET,sGPt0,li,lj,comp)
        @synchronize
        @inbounds s1[li,lj]=sIER[li,lj]+dotc(sQ,sGPr0,li,lj,comp); @inbounds s3[li,lj]=dotc(sIET,sGPr0,li,lj,comp)
        @synchronize
        @inbounds s2[li,lj]=dotc(sR1,sIET,li,lj,comp)+dotc(s1,sT0,li,lj,comp)
        @synchronize
        @inbounds tmp6=sIER[li,lj]+dotc(s3,sT0,li,lj,comp)+dotc(sTTGP,s2,li,lj,comp)
        @inbounds IET[li,lj,n1,Δn]=sT5[li,lj]; @inbounds IER[li,lj,n1,Δn]=tmp6
    end
end

# ─── INTERACTION SI_11 loop-1, tmpieT⁻⁻ output (interaction_inelastic.jl:413-418)
# Single-output template (12 buffers — within the safe ceiling). The other 5
# outputs (tmpieR⁻⁺, tmpieJ₀⁻ in loop1; tmpieT⁺⁺/tmpieR⁺⁻/tmpieJ₀⁺ in loop2) follow
# the same shape, each its own kernel. tmp_inv/T01_inv are precomputed (batch_inv!).
@kernel function inter1_T!(oTmm, @Const(T01), @Const(TINV), @Const(r), @Const(Rpm), @Const(tmm),
                           @Const(ier), @Const(ieRpm), @Const(ieTmm), @Const(ietmm), @Const(shifts), nSpec, comp::Bool)
    _i,_j, n1, Δn = @index(Global, NTuple); li, lj, _, _ = @index(Local, NTuple)
    st01=@localmem eltype(oTmm) (NQ,NQ); stinv=@localmem eltype(oTmm) (NQ,NQ); sr1=@localmem eltype(oTmm) (NQ,NQ)
    srpm=@localmem eltype(oTmm) (NQ,NQ); stm=@localmem eltype(oTmm) (NQ,NQ)
    sie=@localmem eltype(oTmm) (NQ,NQ); siRpm=@localmem eltype(oTmm) (NQ,NQ); siTmm=@localmem eltype(oTmm) (NQ,NQ); sitmm=@localmem eltype(oTmm) (NQ,NQ)
    sW=@localmem eltype(oTmm) (NQ,NQ); sWt=@localmem eltype(oTmm) (NQ,NQ); sA=@localmem eltype(oTmm) (NQ,NQ)
    Δ=@inbounds shifts[Δn]; n0=n1+Δ
    if 1<=n0<=nSpec
        @inbounds begin
            st01[li,lj]=T01[li,lj,n1]; stinv[li,lj]=TINV[li,lj,n0]; sr1[li,lj]=r[li,lj,n1]; srpm[li,lj]=Rpm[li,lj,n0]; stm[li,lj]=tmm[li,lj,n0]
            sie[li,lj]=ier[li,lj,n1,Δn]; siRpm[li,lj]=ieRpm[li,lj,n1,Δn]; siTmm[li,lj]=ieTmm[li,lj,n1,Δn]; sitmm[li,lj]=ietmm[li,lj,n1,Δn]
        end
        @synchronize
        @inbounds sA[li,lj]=dotc(sie,srpm,li,lj,comp)+dotc(sr1,siRpm,li,lj,comp)   # ier⊠R⁺⁻ + r⊠ieR⁺⁻
        @synchronize
        @inbounds sW[li,lj]=dotc(st01,sA,li,lj,comp)+siTmm[li,lj]                  # W = T01_inv⊠A + ieT⁻⁻
        @synchronize
        @inbounds sWt[li,lj]=dotc(sW,stinv,li,lj,comp)                            # Wt = W⊠tmp_inv
        @synchronize
        @inbounds oTmm[li,lj,n1,Δn]=dotc(st01,sitmm,li,lj,comp)+dotc(sWt,stm,li,lj,comp)  # T01_inv⊠iet⁻⁻ + Wt⊠t⁻⁻
    end
end

# ─── references (batched_mul; COPY old slices to avoid in-place view aliasing) ──
function dbl_nloop2_ref!(IET,IER,TTGP,GP,R,T,shifts)
    nSpec=size(IET,3)
    for Δn in 1:size(IET,4)
        n1,n0=get_nn(nSpec,shifts[Δn])
        @views begin
            iet=copy(IET[:,:,n1,Δn]); ier=copy(IER[:,:,n1,Δn])
            ttgp=TTGP[:,:,n1]; gp0=GP[:,:,n0]; r0=R[:,:,n0]; r1=R[:,:,n1]; t0=T[:,:,n0]
            Q=ier⊠r0 .+ r1⊠ier; GPt0=gp0⊠t0; GPr0=gp0⊠r0
            IET[:,:,n1,Δn] .= ttgp⊠(iet .+ Q⊠GPt0) .+ iet⊠GPt0
            IER[:,:,n1,Δn] .= ier .+ iet⊠GPr0⊠t0 .+ ttgp⊠(r1⊠iet .+ (ier .+ Q⊠GPr0)⊠t0)
        end
    end
end
function inter1_T_ref!(oTmm,T01,TINV,r,Rpm,tmm,ier,ieRpm,ieTmm,ietmm,shifts)
    nSpec=size(ier,3)
    for Δn in 1:size(ier,4)
        n1,n0=get_nn(nSpec,shifts[Δn])
        @views begin
            W = T01[:,:,n1]⊠(ier[:,:,n1,Δn]⊠Rpm[:,:,n0] .+ r[:,:,n1]⊠ieRpm[:,:,n1,Δn]) .+ ieTmm[:,:,n1,Δn]
            Wt = W ⊠ TINV[:,:,n0]
            oTmm[:,:,n1,Δn] .= T01[:,:,n1]⊠ietmm[:,:,n1,Δn] .+ Wt⊠tmm[:,:,n0]
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    nSpec=parse(Int,get(ENV,"NSPEC","10000")); nRaman=200; sh=collect(-100:99); shd=CuArray(sh); Random.seed!(7)
    rel(a,b)=maximum(abs.(Float64.(Array(a)).-b))/maximum(abs.(b))
    # doubling
    m3()=0.05.*rand(Float64,NQ,NQ,nSpec); m4()=0.05.*rand(Float64,NQ,NQ,nSpec,nRaman)
    TTGP,GP,R,T,IET0,IER0 = m3(),m3(),m3(),m3(),m4(),m4()
    Iref=copy(IET0);Jref=copy(IER0); dbl_nloop2_ref!(Iref,Jref,TTGP,GP,R,T,sh)
    for FT in (Float64,Float32)
        g(x)=CuArray(FT.(x)); I=g(IET0);J=g(IER0)
        dbl_nloop2!(KA.get_backend(I),(NQ,NQ,1,1))(I,J,g(TTGP),g(GP),g(R),g(T),shd,nSpec,false;ndrange=(NQ,NQ,nSpec,nRaman)); KA.synchronize(KA.get_backend(I))
        @printf("doubling   %-7s vs f64 ref: IET %.2e IER %.2e\n", FT, rel(I,Iref), rel(J,Jref))
    end
end
