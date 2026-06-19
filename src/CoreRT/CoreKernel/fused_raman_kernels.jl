#=
fused_raman_kernels.jl — fused KernelAbstractions kernels that collapse the
per-Δn batched_mul host loops in the RRS adding-doubling solver into one launch
per (n₁,Δn). One workgroup per (n₁,Δn); N×N threads (one per output element);
all per-block intermediates live in @localmem; common subexpressions are computed
once; the Raman shift n₀ = n₁ + i_λ₁λ₀[Δn] is resolved in-kernel.

Validated f64-vs-batched_mul to machine precision (2e-16; f32 ~1e-7) on L40S +
A100; ~6× faster than the batched_mul host loop at nSpec=10000 (NquadN=15).
NOT bit-identical to the batched_mul path (different accumulation order) — RRS
goldens must be regenerated for configs that take the fused path.

⚠ FALLBACK: these kernels use static @localmem, capped at 48 KB shared / 1024
threads per block. `_fused_doubling_fits(N, FT)` gates the launch; above the cap
(e.g. NquadN≥21 in f64, ≥30 in f32, or NquadN=48 IQU) the caller falls back to
the batched_mul host loop. The CPU backend cannot run these (scalar locals don't
survive KA @synchronize splits) — CPU always takes the fallback.
=#

# N-generic dot of row i of X with col j of Y; comp=true ⇒ Neumaier compensated sum.
# Loop bound comes from the statically-sized @localmem operand (compile-time N).
@inline function _ie_dotc(X, Y, i, j, comp)
    T = eltype(X); s = zero(T); c = zero(T)
    @inbounds for k in 1:size(X, 2)
        p = X[i,k] * Y[k,j]
        comp ? (t = s + p; c += ifelse(abs(s) >= abs(p), (s-t)+p, (p-t)+s); s = t) : (s += p)
    end
    return comp ? s + c : s
end

# Number of N×N @localmem buffers used by dbl_nloop2_fused! (for the fit check).
const _DBL_NLOOP2_NBUF = 14

"""
    _fused_doubling_fits(N, FT) -> Bool

True when the fused doubling n-loop-2 kernel fits the static-shared (48 KB) and
threads-per-block (1024) limits for matrix dimension `N` and float type `FT`.
"""
@inline _fused_doubling_fits(N::Int, ::Type{FT}) where {FT} =
    (_DBL_NLOOP2_NBUF * N * N * sizeof(FT) <= 48 * 1024) && (N * N <= 1024)

# Fused replacement for the RRS doubling "n loop 2" (doubling_inelastic.jl).
# Computes tmp5→iet⁺⁺ and tmp6→ier⁻⁺ for every (n₁,Δn), in place. Both outputs
# read the OLD iet/ier (loaded once into shared), so the in-place write is
# race-free per (n₁,Δn) block.
@kernel function dbl_nloop2_fused!(IET, IER, @Const(TTGP), @Const(GP), @Const(R), @Const(T),
                                   @Const(shifts), nSpec, ::Val{N}, comp::Bool) where {N}
    _i, _j, n1, Δn = @index(Global, NTuple)
    li, lj, _, _   = @index(Local, NTuple)
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

"""
    apply_fused_doubling_nloop2!(iet⁺⁺, ier⁻⁺, tt⁺⁺_gp_refl, gp_refl, r⁻⁺, t⁺⁺,
                                 i_λ₁λ₀_dev, architecture; comp=false)

Launch `dbl_nloop2_fused!` over all (n₁,Δn). `i_λ₁λ₀_dev` must already be on the
target device. Caller must have checked `_fused_doubling_fits(N, FT)`.
"""
function apply_fused_doubling_nloop2!(iet⁺⁺, ier⁻⁺, tt⁺⁺_gp_refl, gp_refl, r⁻⁺, t⁺⁺,
                                      i_λ₁λ₀_dev, architecture; comp::Bool=false)
    N, _, nSpec, nRaman = size(iet⁺⁺)
    dev = devi(architecture)
    dbl_nloop2_fused!(dev, (N, N, 1, 1))(iet⁺⁺, ier⁻⁺, tt⁺⁺_gp_refl, gp_refl, r⁻⁺, t⁺⁺,
                                         i_λ₁λ₀_dev, nSpec, Val(N), comp; ndrange=(N, N, nSpec, nRaman))
    return nothing
end
