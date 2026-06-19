#=
================================================================================
 fused_raman_kernels.jl — fused GPU kernels for the RRS adding–doubling solver
================================================================================

WHAT THIS FILE IS
-----------------
Custom KernelAbstractions kernels that collapse the per-Δn `batched_mul` host
loops of the rotational-Raman (RRS) adding–doubling solver into a SINGLE GPU
kernel launch each. They are drop-in replacements (same math, machine-precision
equal in f64) for two hot inner loops:

  • doubling_inelastic.jl  — the doubling "n loop 2" R/T update  (iet⁺⁺, ier⁻⁺)
  • interaction_inelastic.jl — the ScatteringInterface_11 per-Δn matrix outputs
                               (tmpieR⁻⁺, tmpieT⁻⁻, tmpieT⁺⁺, tmpieR⁺⁻)

WHY THEY EXIST (rationale)
--------------------------
The RRS solver multiplies tiny NquadN×NquadN matrices (NquadN = Nstreams·nStokes,
typically 15–48) for every spectral point (nSpec ~ 10³–10⁴) and every Raman
coupling offset (nRaman ~ 200). Two performance problems were found by profiling
(L40S, see test/benchmarks/raman_fused_kernels.jl and the project memory):

  1. The per-Δn loops issue ~13 tiny `batched_mul` + broadcast launches PER Δn,
     i.e. tens of thousands of dependent kernel launches. The GPU sat at ~16%
     utilisation — launch-latency bound, not compute bound. (An earlier, separate
     fix replaced `=` slice-assignment with `.=` to kill a 2.8 ms-per-slice
     `gpu_setindex_kernel`; that gave ~5.7×. THIS file is the next layer.)
  2. Every intermediate of the matmul chain round-tripped through HBM.

A fused kernel assigns ONE thread-block to each (n₁, Δn) pair, with NquadN×NquadN
threads (one per output-matrix element). It loads the handful of NquadN×NquadN
operand matrices into shared memory ONCE, computes the whole chain there (common
subexpressions — e.g. `Q`, `GPt0`, `GPr0`, `W`, `Wt` — computed a single time,
intermediates never leaving shared memory), resolves the Raman wavelength shift
`n₀ = n₁ + i_λ₁λ₀[Δn]` with in-kernel index arithmetic, and writes the result
once. Measured: doubling n-loop-2 6.1× @ nSpec=10000 (≈100% GPU-busy);
interaction matrix outputs validated to machine precision.

HOW THEY WORK
-------------
  • ndrange = (N, N, nSpec, nRaman); workgroup = (N, N, 1, 1) ⇒ one block / (n₁,Δn).
  • `_ie_dotc(X, Y, i, j, comp)` is an N-generic row·col dot of two shared NxN
    matrices. `comp=true` switches on Neumaier compensated summation (≈3× tighter
    f32 accuracy for ~7–13% more time; default `false`).
  • Kernels are parameterised by `Val{N}` so `@localmem (N,N)` specialises per
    NquadN. `N` is taken from `size(.,1)` of the operands at the call site.

WHAT THEY REPLACE / WHAT STAYS
------------------------------
Replaced (GPU, when they fit — see FALLBACK):
  • doubling "n loop 2"            → `dbl_nloop2_fused!`
  • interaction SI_11 tmpieR⁻⁺     → `ie_inter_Rmp!`     (loop 1)
  • interaction SI_11 tmpieT⁻⁻     → `ie_inter_Tmm!`     (loop 1)
  • interaction SI_11 tmpieT⁺⁺     → `ie_inter_Tpp!`     (loop 2)
  • interaction SI_11 tmpieR⁺⁻     → `ie_inter_Rpm2!`    (loop 2)
NOT fused (still batched_mul, intentionally): the doubling SFI source loop and
the interaction SI_11 source-vector outputs (tmpieJ₀⁺/tmpieJ₀⁻) — those are
matrix×vector (NquadN×1) and a separate follow-up; the `batch_inv!` inverses
(<0.1% of time); the elastic R/T updates; the VS and multisensor paths.

FALLBACK (this is mandatory, not optional)
------------------------------------------
The kernels use STATIC `@localmem`, hard-capped at 48 KB shared memory and 1024
threads per block. `_fused_*_fits(N, FT)` gates every launch; when it returns
`false` the caller runs the original `batched_mul` host loop unchanged. The
binding limit is SHARED MEMORY (tighter than threads):
  • doubling (14 buffers):     fused for N ≤ 20 (f64) / N ≤ 29 (f32)
  • interaction (12 buffers):  fused for N ≤ 23 (f64) / N ≤ 32 (f32)
So NquadN = 30 (e.g. the wide-window production case) and NquadN = 48 (16-stream
IQU) automatically FALL BACK. The KA CPU backend cannot run these kernels at all
(scalar locals do not survive KA's @synchronize loop-splitting) — CPU always
falls back. The fused path is therefore a GPU bonus for small NquadN with a
guaranteed, byte-for-byte-unchanged batched_mul path everywhere else.

ACCURACY / GOLDENS
------------------
Fused vs batched_mul agree to MACHINE PRECISION in f64 (≤ ~5e-16) and to ~1e-7
in f32 — but NOT bit-identically (different accumulation order). Configs that
take the fused path therefore need their RRS goldens regenerated; the Phase1b
gate uses rtol=0.02 which already absorbs the f32 difference. Validation is done
f64-vs-batched_mul on both the L40S (f32 box) and curry's A100 (real f64), which
agree bit-for-bit — see project memory `project_raman_setindex_bottleneck`.

MAINTAINER GOTCHAS (KA quirks that bit us — keep these in mind when editing)
---------------------------------------------------------------------------
  • `@localmem` element type MUST be inlined: `@localmem eltype(X) (N,N)` — a
    local `E = eltype(X)` does NOT survive macro expansion.
  • `@localmem` declarations must be at kernel top level (not inside `if`).
  • `return` is illegal in a KA kernel — guard the body with `if 1<=n0<=nSpec`.
  • DO NOT exceed ~14 `@localmem` buffers in one kernel: an 18-buffer monolith
    silently corrupted outputs (flaky NaN/race). Split into per-output kernels
    and read element-wise-only operands straight from global instead of shared.
  • When writing a batched_mul *reference*, COPY the old slices first
    (`iet = copy(IET[:,:,n1,Δn])`): writing IET in place before computing IER
    aliases a view onto updated data (this masqueraded as a 7.7% "error").
================================================================================
=#

# row-i · col-j dot of two N×N @localmem matrices (N = size(X,2), compile-time).
# comp=true ⇒ Neumaier compensated summation (tighter f32 accuracy).
@inline function _ie_dotc(X, Y, i, j, comp)
    T = eltype(X); s = zero(T); c = zero(T)
    @inbounds for k in 1:size(X, 2)
        p = X[i,k] * Y[k,j]
        comp ? (t = s + p; c += ifelse(abs(s) >= abs(p), (s-t)+p, (p-t)+s); s = t) : (s += p)
    end
    return comp ? s + c : s
end

# Static-shared (48 KB) / 1024-thread fit test for `nbuf` N×N buffers of type FT.
@inline _fused_fits(nbuf::Int, N::Int, ::Type{FT}) where {FT} =
    (nbuf * N * N * sizeof(FT) <= 48 * 1024) && (N * N <= 1024)
@inline _fused_doubling_fits(N::Int, ::Type{FT}) where {FT}    = _fused_fits(14, N, FT)
@inline _fused_interaction_fits(N::Int, ::Type{FT}) where {FT} = _fused_fits(12, N, FT)

# ───────────────────────── DOUBLING "n loop 2" ────────────────────────────────
# Replaces doubling_inelastic.jl:120-143. Writes iet⁺⁺ (tmp5) and ier⁻⁺ (tmp6)
# in place; both read the OLD iet/ier (loaded once into shared) so the per-(n₁,Δn)
# in-place write is race-free.
@kernel function dbl_nloop2_fused!(IET, IER, @Const(TTGP), @Const(GP), @Const(R), @Const(T),
                                   @Const(shifts), nSpec, ::Val{N}, comp::Bool) where {N}
    _i, _j, n1, Δn = @index(Global, NTuple); li, lj, _, _ = @index(Local, NTuple)
    sTTGP=@localmem eltype(IET) (N,N); sGP0=@localmem eltype(IET) (N,N); sR0=@localmem eltype(IET) (N,N)
    sR1=@localmem eltype(IET) (N,N); sT0=@localmem eltype(IET) (N,N); sIET=@localmem eltype(IET) (N,N); sIER=@localmem eltype(IET) (N,N)
    sQ=@localmem eltype(IET) (N,N); sGPt0=@localmem eltype(IET) (N,N); sGPr0=@localmem eltype(IET) (N,N)
    s1=@localmem eltype(IET) (N,N); s2=@localmem eltype(IET) (N,N); s3=@localmem eltype(IET) (N,N); sT5=@localmem eltype(IET) (N,N)
    Δ = @inbounds shifts[Δn]; n0 = n1 + Δ
    if 1 <= n0 <= nSpec
        @inbounds begin
            sTTGP[li,lj]=TTGP[li,lj,n1]; sGP0[li,lj]=GP[li,lj,n0]; sR0[li,lj]=R[li,lj,n0]
            sR1[li,lj]=R[li,lj,n1]; sT0[li,lj]=T[li,lj,n0]; sIET[li,lj]=IET[li,lj,n1,Δn]; sIER[li,lj]=IER[li,lj,n1,Δn]
        end
        @synchronize
        @inbounds sQ[li,lj]=_ie_dotc(sIER,sR0,li,lj,comp)+_ie_dotc(sR1,sIER,li,lj,comp)
        @inbounds sGPt0[li,lj]=_ie_dotc(sGP0,sT0,li,lj,comp); @inbounds sGPr0[li,lj]=_ie_dotc(sGP0,sR0,li,lj,comp)
        @synchronize
        @inbounds s1[li,lj]=sIET[li,lj]+_ie_dotc(sQ,sGPt0,li,lj,comp)
        @synchronize
        @inbounds sT5[li,lj]=_ie_dotc(sTTGP,s1,li,lj,comp)+_ie_dotc(sIET,sGPt0,li,lj,comp)
        @synchronize
        @inbounds s1[li,lj]=sIER[li,lj]+_ie_dotc(sQ,sGPr0,li,lj,comp); @inbounds s3[li,lj]=_ie_dotc(sIET,sGPr0,li,lj,comp)
        @synchronize
        @inbounds s2[li,lj]=_ie_dotc(sR1,sIET,li,lj,comp)+_ie_dotc(s1,sT0,li,lj,comp)
        @synchronize
        @inbounds tmp6=sIER[li,lj]+_ie_dotc(s3,sT0,li,lj,comp)+_ie_dotc(sTTGP,s2,li,lj,comp)
        @inbounds IET[li,lj,n1,Δn]=sT5[li,lj]; @inbounds IER[li,lj,n1,Δn]=tmp6
    end
end

function apply_fused_doubling_nloop2!(iet⁺⁺, ier⁻⁺, tt⁺⁺_gp_refl, gp_refl, r⁻⁺, t⁺⁺,
                                      i_λ₁λ₀_dev, architecture; comp::Bool=false)
    N, _, nSpec, nRaman = size(iet⁺⁺)
    dbl_nloop2_fused!(devi(architecture), (N, N, 1, 1))(iet⁺⁺, ier⁻⁺, tt⁺⁺_gp_refl, gp_refl, r⁻⁺, t⁺⁺,
        i_λ₁λ₀_dev, nSpec, Val(N), comp; ndrange=(N, N, nSpec, nRaman))
    return nothing
end

# ───────────────── INTERACTION SI_11 matrix outputs (per Δn) ──────────────────
# Each writes one tmpie* workspace array from composite/added ie* + elastic ops.
# Common subexpr per loop: W = T_inv⊠(ie⊠R + r⊠ieR) + ieT ; Wt = W⊠tmp_inv.
# tmp_inv / T01_inv / T21_inv are precomputed by the caller (batch_inv! + ⊠).

# loop1 tmpieR⁻⁺ = ieR⁻⁺ + T01⊠(ier⊠T⁺⁺ + r⊠ieT⁺⁺) + Wt⊠r₀⊠T⁺⁺   (ieT⁻⁻,ieR⁻⁺ read from global)
@kernel function ie_inter_Rmp!(o, @Const(T01), @Const(TINV), @Const(r), @Const(Rpm), @Const(Tpp),
                               @Const(ier), @Const(ieRpm), @Const(ieTpp), @Const(ieTmm), @Const(ieRmp),
                               @Const(shifts), nSpec, ::Val{N}, comp::Bool) where {N}
    _i,_j,n1,Δn=@index(Global,NTuple); li,lj,_,_=@index(Local,NTuple)
    sT01=@localmem eltype(o) (N,N); sTINV=@localmem eltype(o) (N,N); sr1=@localmem eltype(o) (N,N); sr0=@localmem eltype(o) (N,N)
    sRpm=@localmem eltype(o) (N,N); sTpp=@localmem eltype(o) (N,N); sie=@localmem eltype(o) (N,N); sieRpm=@localmem eltype(o) (N,N); sieTpp=@localmem eltype(o) (N,N)
    sW=@localmem eltype(o) (N,N); sWt=@localmem eltype(o) (N,N); sA=@localmem eltype(o) (N,N)
    Δ=@inbounds shifts[Δn]; n0=n1+Δ
    if 1<=n0<=nSpec
        @inbounds begin
            sT01[li,lj]=T01[li,lj,n1]; sTINV[li,lj]=TINV[li,lj,n0]; sr1[li,lj]=r[li,lj,n1]; sr0[li,lj]=r[li,lj,n0]
            sRpm[li,lj]=Rpm[li,lj,n0]; sTpp[li,lj]=Tpp[li,lj,n0]; sie[li,lj]=ier[li,lj,n1,Δn]; sieRpm[li,lj]=ieRpm[li,lj,n1,Δn]; sieTpp[li,lj]=ieTpp[li,lj,n1,Δn]
        end
        @synchronize
        @inbounds sA[li,lj]=_ie_dotc(sie,sRpm,li,lj,comp)+_ie_dotc(sr1,sieRpm,li,lj,comp)
        @synchronize
        @inbounds sW[li,lj]=_ie_dotc(sT01,sA,li,lj,comp)+ieTmm[li,lj,n1,Δn]
        @synchronize
        @inbounds sWt[li,lj]=_ie_dotc(sW,sTINV,li,lj,comp)
        @inbounds sA[li,lj]=_ie_dotc(sie,sTpp,li,lj,comp)+_ie_dotc(sr1,sieTpp,li,lj,comp)
        @synchronize
        @inbounds sW[li,lj]=_ie_dotc(sWt,sr0,li,lj,comp)
        @synchronize
        @inbounds o[li,lj,n1,Δn]=ieRmp[li,lj,n1,Δn]+_ie_dotc(sT01,sA,li,lj,comp)+_ie_dotc(sW,sTpp,li,lj,comp)
    end
end

# loop1 tmpieT⁻⁻ = T01⊠iet⁻⁻ + Wt⊠t⁻⁻
@kernel function ie_inter_Tmm!(o, @Const(T01), @Const(TINV), @Const(r), @Const(Rpm), @Const(tmm),
                               @Const(ier), @Const(ieRpm), @Const(ieTmm), @Const(ietmm),
                               @Const(shifts), nSpec, ::Val{N}, comp::Bool) where {N}
    _i,_j,n1,Δn=@index(Global,NTuple); li,lj,_,_=@index(Local,NTuple)
    sT01=@localmem eltype(o) (N,N); sTINV=@localmem eltype(o) (N,N); sr1=@localmem eltype(o) (N,N); sRpm=@localmem eltype(o) (N,N); stmm=@localmem eltype(o) (N,N)
    sie=@localmem eltype(o) (N,N); sieRpm=@localmem eltype(o) (N,N); sieTmm=@localmem eltype(o) (N,N); sietmm=@localmem eltype(o) (N,N)
    sW=@localmem eltype(o) (N,N); sWt=@localmem eltype(o) (N,N); sA=@localmem eltype(o) (N,N)
    Δ=@inbounds shifts[Δn]; n0=n1+Δ
    if 1<=n0<=nSpec
        @inbounds begin
            sT01[li,lj]=T01[li,lj,n1]; sTINV[li,lj]=TINV[li,lj,n0]; sr1[li,lj]=r[li,lj,n1]; sRpm[li,lj]=Rpm[li,lj,n0]; stmm[li,lj]=tmm[li,lj,n0]
            sie[li,lj]=ier[li,lj,n1,Δn]; sieRpm[li,lj]=ieRpm[li,lj,n1,Δn]; sieTmm[li,lj]=ieTmm[li,lj,n1,Δn]; sietmm[li,lj]=ietmm[li,lj,n1,Δn]
        end
        @synchronize
        @inbounds sA[li,lj]=_ie_dotc(sie,sRpm,li,lj,comp)+_ie_dotc(sr1,sieRpm,li,lj,comp)
        @synchronize
        @inbounds sW[li,lj]=_ie_dotc(sT01,sA,li,lj,comp)+sieTmm[li,lj]
        @synchronize
        @inbounds sWt[li,lj]=_ie_dotc(sW,sTINV,li,lj,comp)
        @synchronize
        @inbounds o[li,lj,n1,Δn]=_ie_dotc(sT01,sietmm,li,lj,comp)+_ie_dotc(sWt,stmm,li,lj,comp)
    end
end

# loop2 tmpieT⁺⁺ = T21⊠ieT⁺⁺ + W2t⊠T⁺⁺ ; W2=T21⊠(ieR⁺⁻⊠r₀+R⁺⁻⊠ier⁻⁺)+iet⁺⁺(global)
@kernel function ie_inter_Tpp!(o, @Const(T21), @Const(TINV), @Const(r), @Const(Rpm), @Const(Tpp),
                               @Const(ier), @Const(ieRpm), @Const(ieTpp), @Const(iet),
                               @Const(shifts), nSpec, ::Val{N}, comp::Bool) where {N}
    _i,_j,n1,Δn=@index(Global,NTuple); li,lj,_,_=@index(Local,NTuple)
    sT21=@localmem eltype(o) (N,N); sTINV=@localmem eltype(o) (N,N); sr0=@localmem eltype(o) (N,N); sRpm1=@localmem eltype(o) (N,N)
    sTpp=@localmem eltype(o) (N,N); sie=@localmem eltype(o) (N,N); sieRpm=@localmem eltype(o) (N,N); sieTpp=@localmem eltype(o) (N,N)
    sW=@localmem eltype(o) (N,N); sWt=@localmem eltype(o) (N,N); sA=@localmem eltype(o) (N,N)
    Δ=@inbounds shifts[Δn]; n0=n1+Δ
    if 1<=n0<=nSpec
        @inbounds begin
            sT21[li,lj]=T21[li,lj,n1]; sTINV[li,lj]=TINV[li,lj,n0]; sr0[li,lj]=r[li,lj,n0]; sRpm1[li,lj]=Rpm[li,lj,n1]
            sTpp[li,lj]=Tpp[li,lj,n0]; sie[li,lj]=ier[li,lj,n1,Δn]; sieRpm[li,lj]=ieRpm[li,lj,n1,Δn]; sieTpp[li,lj]=ieTpp[li,lj,n1,Δn]
        end
        @synchronize
        @inbounds sA[li,lj]=_ie_dotc(sieRpm,sr0,li,lj,comp)+_ie_dotc(sRpm1,sie,li,lj,comp)
        @synchronize
        @inbounds sW[li,lj]=_ie_dotc(sT21,sA,li,lj,comp)+iet[li,lj,n1,Δn]
        @synchronize
        @inbounds sWt[li,lj]=_ie_dotc(sW,sTINV,li,lj,comp)
        @synchronize
        @inbounds o[li,lj,n1,Δn]=_ie_dotc(sT21,sieTpp,li,lj,comp)+_ie_dotc(sWt,sTpp,li,lj,comp)
    end
end

# loop2 tmpieR⁺⁻ = ier⁺⁻(global) + T21⊠(ieR⁺⁻⊠t⁻⁻ + R⁺⁻⊠iet⁻⁻) + W2t⊠R⁺⁻₀⊠t⁻⁻
@kernel function ie_inter_Rpm2!(o, @Const(T21), @Const(TINV), @Const(r), @Const(Rpm), @Const(tmm),
                                @Const(ier), @Const(ieRpm), @Const(ietmm), @Const(iet), @Const(ierpm),
                                @Const(shifts), nSpec, ::Val{N}, comp::Bool) where {N}
    _i,_j,n1,Δn=@index(Global,NTuple); li,lj,_,_=@index(Local,NTuple)
    sT21=@localmem eltype(o) (N,N); sTINV=@localmem eltype(o) (N,N); sr0=@localmem eltype(o) (N,N); sRpm1=@localmem eltype(o) (N,N); sRpm0=@localmem eltype(o) (N,N)
    stmm=@localmem eltype(o) (N,N); sie=@localmem eltype(o) (N,N); sieRpm=@localmem eltype(o) (N,N); sietmm=@localmem eltype(o) (N,N)
    sW=@localmem eltype(o) (N,N); sWt=@localmem eltype(o) (N,N); sA=@localmem eltype(o) (N,N)
    Δ=@inbounds shifts[Δn]; n0=n1+Δ
    if 1<=n0<=nSpec
        @inbounds begin
            sT21[li,lj]=T21[li,lj,n1]; sTINV[li,lj]=TINV[li,lj,n0]; sr0[li,lj]=r[li,lj,n0]; sRpm1[li,lj]=Rpm[li,lj,n1]; sRpm0[li,lj]=Rpm[li,lj,n0]
            stmm[li,lj]=tmm[li,lj,n0]; sie[li,lj]=ier[li,lj,n1,Δn]; sieRpm[li,lj]=ieRpm[li,lj,n1,Δn]; sietmm[li,lj]=ietmm[li,lj,n1,Δn]
        end
        @synchronize
        @inbounds sA[li,lj]=_ie_dotc(sieRpm,sr0,li,lj,comp)+_ie_dotc(sRpm1,sie,li,lj,comp)
        @synchronize
        @inbounds sW[li,lj]=_ie_dotc(sT21,sA,li,lj,comp)+iet[li,lj,n1,Δn]
        @synchronize
        @inbounds sWt[li,lj]=_ie_dotc(sW,sTINV,li,lj,comp)
        @inbounds sA[li,lj]=_ie_dotc(sieRpm,stmm,li,lj,comp)+_ie_dotc(sRpm1,sietmm,li,lj,comp)
        @synchronize
        @inbounds sW[li,lj]=_ie_dotc(sWt,sRpm0,li,lj,comp)
        @synchronize
        @inbounds o[li,lj,n1,Δn]=ierpm[li,lj,n1,Δn]+_ie_dotc(sT21,sA,li,lj,comp)+_ie_dotc(sW,stmm,li,lj,comp)
    end
end

# launchers (caller must have checked _fused_interaction_fits(N, FT))
function apply_fused_interaction_Rmp!(o, T01, TINV, r, Rpm, Tpp, ier, ieRpm, ieTpp, ieTmm, ieRmp, shd, architecture; comp::Bool=false)
    N,_,nSpec,nRaman = size(ier)
    ie_inter_Rmp!(devi(architecture),(N,N,1,1))(o,T01,TINV,r,Rpm,Tpp,ier,ieRpm,ieTpp,ieTmm,ieRmp,shd,nSpec,Val(N),comp;ndrange=(N,N,nSpec,nRaman)); return nothing
end
function apply_fused_interaction_Tmm!(o, T01, TINV, r, Rpm, tmm, ier, ieRpm, ieTmm, ietmm, shd, architecture; comp::Bool=false)
    N,_,nSpec,nRaman = size(ier)
    ie_inter_Tmm!(devi(architecture),(N,N,1,1))(o,T01,TINV,r,Rpm,tmm,ier,ieRpm,ieTmm,ietmm,shd,nSpec,Val(N),comp;ndrange=(N,N,nSpec,nRaman)); return nothing
end
function apply_fused_interaction_Tpp!(o, T21, TINV, r, Rpm, Tpp, ier, ieRpm, ieTpp, iet, shd, architecture; comp::Bool=false)
    N,_,nSpec,nRaman = size(ier)
    ie_inter_Tpp!(devi(architecture),(N,N,1,1))(o,T21,TINV,r,Rpm,Tpp,ier,ieRpm,ieTpp,iet,shd,nSpec,Val(N),comp;ndrange=(N,N,nSpec,nRaman)); return nothing
end
function apply_fused_interaction_Rpm2!(o, T21, TINV, r, Rpm, tmm, ier, ieRpm, ietmm, iet, ierpm, shd, architecture; comp::Bool=false)
    N,_,nSpec,nRaman = size(ier)
    ie_inter_Rpm2!(devi(architecture),(N,N,1,1))(o,T21,TINV,r,Rpm,tmm,ier,ieRpm,ietmm,iet,ierpm,shd,nSpec,Val(N),comp;ndrange=(N,N,nSpec,nRaman)); return nothing
end
