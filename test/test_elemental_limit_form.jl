# =========================================================================
# Elemental transmission/source factors near coincident ordinates
# =========================================================================
#
# The finite-δ transmission t⁺⁺ and the direct-solar source J₀⁺ (Fell 1997,
# Eqs. 1.53–1.54; SF2023-II Eqs. 10–11) contain
#
#     μⱼ/(μᵢ - μⱼ) · (exp(-τ/μᵢ) - exp(-τ/μⱼ))                    (elastic)
#     [exp(-τ₁/μᵢ) - exp(-τ₀/μⱼ)] / (μᵢ/μⱼ - τ₁/τ₀)                (Raman)
#
# Dividing an exponential difference by the cosine difference amplifies the
# rounding of τ/μ by μ/|μᵢ - μⱼ|: up to about 5e-19/δ in Float64 and 2e-10/δ
# in Float32 at the RT output (δ = |μᵢ - μ₀|; an estimate, scene-dependent),
# outside the 4·eps window in which `rt_set_streams` merges ordinates. The
# kernels now evaluate both forms as a divided difference of exp(-x)
# (`expdiff_ratio`, midpoint form exp(-s)·sinh(x)/x for |b - a| ≤ 2), which
# never divides by the difference (decision 14e). These tests check, in
# Float64 and Float32:
#   (i)   accuracy against BigFloat for δ from 1e-12 to 1e-1, for the helpers
#         and for every elastic, linearized, canopy and Raman kernel that
#         carries the factor;
#   (ii)  continuity across the merge window (kernel factor and TOA radiance);
#   (iii) the linearized tangent against BigFloat derivatives, and the RT
#         Jacobian against finite differences, near μ₀ ≈ μᵢ;
#   (iv)  agreement with the former expression away from coincidence.
# The kernels are KernelAbstractions kernels; with
# VSMARTMOM_JACOBIAN_GPU_TEST=true the kernel checks also run on CUDA.

module ElementalLimitFormTests

using Test, LinearAlgebra, Logging, YAML
using ForwardDiff
using vSmartMOM, vSmartMOM.CoreRT
import KernelAbstractions as KA

const C = vSmartMOM.CoreRT
setprecision(BigFloat, 256)

# ---- exact references (BigFloat; inputs are the FT values, exactly) -------
bigDD(a, b) = a == b ? exp(-a) : (exp(-a) - exp(-b)) / (b - a)
ref_t(τ, μi, μj) = (τ = big(τ); μi = big(μi); μj = big(μj); (τ / μi) * bigDD(τ / μi, τ / μj))
function ref_dt(τ, μi, μj)
    τ = big(τ); μi = big(μi); μj = big(μj)
    μi == μj && return exp(-τ / μi) * (1 - τ / μi) / μi
    return μj / (μi - μj) * (exp(-τ / μj) / μj - exp(-τ / μi) / μi)
end
# Raman: [exp(-τ₁/μᵢ) - exp(-τ₀/μⱼ)] / (μᵢ/μⱼ - τ₁/τ₀) and its limit.
ref_raman(τ₁, τ₀, μi, μj) = (τ₁ = big(τ₁); τ₀ = big(τ₀); μi = big(μi); μj = big(μj);
                             (τ₀ / μi) * bigDD(τ₁ / μi, τ₀ / μj))
# Pre-v2.2 elastic expression (the cancellation-prone form).
naive_t(τ, μi, μj) = (μj / (μi - μj)) * C.expdiff_neg(τ / μi, τ / μj)
relerr(x, r) = Float64(abs(big(x) - r) / abs(r))

# Tolerance: a few ulps plus the conditioning of exp(-τ/μ) on its argument.
tol(FT, args...) = 16 * eps(FT) * (1 + maximum(Float64.(args)))

# Ordinate pairs μᵢ = μ₀ + δ (both signs) representable in FT.
function nearby(FT, μ₀, δs)
    out = FT[]
    for δ in δs, s in (1, -1)
        μ = FT(μ₀ + s * δ)
        μ != FT(μ₀) && 0 < μ <= 1 && push!(out, μ)
    end
    return out
end

const DELTAS = 10.0 .^ (-12:-1)
const TAUS = (1e-6, 1e-4, 0.01, 0.5, 2.0, 20.0)

@testset "helpers against BigFloat (i)" begin
    for FT in (Float64, Float32)
        worst = 0.0; worst_naive = 0.0; worst_d = 0.0
        for μ₀ in FT.((0.05, 0.3, 0.5, 0.95)), τ in FT.(TAUS), μi in nearby(FT, μ₀, DELTAS)
            τ / min(μi, μ₀) < 70 || continue                  # stay clear of underflow
            r = ref_t(τ, μi, μ₀)
            e = relerr(C.elem_trans_factor(τ, μi, μ₀), r)
            worst = max(worst, e / tol(FT, τ / min(μi, μ₀)))
            worst_naive = max(worst_naive, relerr(naive_t(τ, μi, μ₀), r))
            # The derivative crosses zero at τ ≈ μ, so measure its error on
            # the scale of its two terms, (exp(-B) + A·DD)/μᵢ.
            a = big(τ) / big(μi); b = big(τ) / big(μ₀)
            scale = (exp(-max(a, b)) + min(a, b) * bigDD(a, b)) / big(μi)
            ed = Float64(abs(big(C.elem_trans_factor_dτ(τ, μi, μ₀)) - ref_dt(τ, μi, μ₀)) / scale)
            worst_d = max(worst_d, ed / tol(FT, τ / min(μi, μ₀)))
        end
        @test worst <= 1
        @test worst_d <= 1
        # The former expression loses accuracy near coincidence ("before").
        @test worst_naive > (FT == Float64 ? 1e-6 : 1e-3)

        # Raman form: every separation of τ₁/μᵢ and τ₀/μⱼ, including the
        # coincidences μᵢ = μⱼ, τ₁ = τ₀ and μᵢ/μⱼ = τ₁/τ₀.
        worst_r = 0.0
        for μj in FT.((0.1, 0.6)), μi in [FT[μj, FT(0.8)]; nearby(FT, μj, DELTAS)],
            τ₀ in FT.((1e-5, 0.02, 3.0)), ρ in (0.0, 1e-10, 1e-7, 1e-4, 0.1)
            for τ₁ in (FT(τ₀ * (1 + ρ)), FT(τ₀ * (μi / μj) * (1 + ρ)))
                r = ref_raman(τ₁, τ₀, μi, μj)
                v = (τ₀ / μi) * C.expdiff_ratio(τ₁ / μi, τ₀ / μj)
                worst_r = max(worst_r, relerr(v, r) / tol(FT, τ₁ / μi, τ₀ / μj))
            end
        end
        @test worst_r <= 1
    end
end

@testset "expdiff_ratio: limit, symmetry, branch switch, Dual" begin
    for FT in (Float64, Float32)
        for a in FT.((0.0, 1e-7, 0.3, 4.0, 60.0))
            @test C.expdiff_ratio(a, a) == exp(-a)
            for b in FT.((a + 1e-3, a + 1.5, a + 2, a + 3, 2a + 1))
                @test C.expdiff_ratio(a, b) === C.expdiff_ratio(b, a)
            end
        end
        # Both sides of the |b - a| = 2 switch agree with BigFloat.
        for a in FT.((0.1, 5.0)), d in (prevfloat(FT(2)), FT(2), nextfloat(FT(2)))
            b = a + d
            @test relerr(C.expdiff_ratio(a, b), bigDD(big(a), big(b))) <= tol(FT, b)
        end
    end
    # ForwardDiff at and near a = b (no quotient to differentiate there):
    # ∂/∂a and ∂/∂b against BigFloat, ∂/∂a = ∂/∂b = -exp(-a)/2 at a = b.
    ddDD(a, b) = (a = big(a); b = big(b); a == b ? (-exp(-a) / 2, -exp(-a) / 2) :
        ((-exp(-a) * (b - a) + (exp(-a) - exp(-b))) / (b - a)^2,
         (exp(-b) * (b - a) - (exp(-a) - exp(-b))) / (b - a)^2))
    for FT in (Float64, Float32), a in FT.((0.0, 1e-5, 0.7, 3.0)),
        r in (0.0, 1e-12, 1e-6, 0.3, 1.9, 2.5)
        b = FT(a + r * (1 + a))
        g = ForwardDiff.gradient(v -> C.expdiff_ratio(v[1], v[2]), [a, b])
        ga, gb = ddDD(a, b)
        sc = exp(-big(min(a, b)))
        # The derivatives inherit the cancellation of a divided difference of
        # the derivative only through the rounding of a and b (≈ eps·max).
        @test abs(g[1] - ga) <= 64eps(FT) * (1 + b) * sc
        @test abs(g[2] - gb) <= 64eps(FT) * (1 + b) * sc
    end
    g = ForwardDiff.gradient(v -> C.elem_trans_factor(v[1], v[2], 0.5), [0.1, 0.5])
    @test g[1] ≈ Float64(ref_dt(0.1, 0.5, 0.5)) rtol = 1e-14
    # Tiny optical depth: the τ-derivative tends to 1/μᵢ (finite).
    for FT in (Float64, Float32)
        d = ForwardDiff.derivative(t -> C.elem_trans_factor(t, FT(0.3), FT(0.75)), FT(1e-25))
        @test d ≈ 1 / FT(0.3) rtol = 4eps(FT)
    end
end

@testset "continuity across the merge window (ii)" begin
    for FT in (Float64, Float32), μ₀ in FT.((0.2, 0.5, 0.9)), τ in FT.((1e-4, 0.1, 1.0))
        on = C.elem_trans_factor(τ, μ₀, μ₀)
        @test on == (τ / μ₀) * exp(-τ / μ₀)              # Fell Eq. 1.54
        μ = μ₀
        for k in 1:12                                    # 1…12 ulps either side
            μ = nextfloat(μ)
            @test abs(C.elem_trans_factor(τ, μ, μ₀) - on) <= 64eps(FT) * on
            @test abs(C.elem_trans_factor(τ, prevfloat(μ₀, k), μ₀) - on) <= 64eps(FT) * on
            @test abs(C.elem_trans_factor_dτ(τ, μ, μ₀) - C.elem_trans_factor_dτ(τ, μ₀, μ₀)) <=
                  64eps(FT) * exp(-τ / μ₀) / μ₀ * (1 + τ / μ₀)
        end
    end
end

@testset "agreement with the former form away from coincidence (iv)" begin
    for FT in (Float64, Float32), n in (3, 8, 16)
        μ, _ = vSmartMOM.Scattering.gauleg(n, zero(FT), one(FT))
        for τ in FT.((1e-4, 0.05, 1.0)), i in eachindex(μ), j in eachindex(μ)
            i == j && continue
            old = naive_t(τ, μ[i], μ[j])
            # The former form carries ≈ eps·(μᵢ+μⱼ)/|μᵢ-μⱼ| of rounding.
            bound = 8eps(FT) * (1 + τ / min(μ[i], μ[j])) * (μ[i] + μ[j]) / abs(μ[i] - μ[j])
            @test abs(C.elem_trans_factor(τ, μ[i], μ[j]) - old) <= bound * abs(old)
        end
    end
end

# ---- kernels ----------------------------------------------------------------
# Kernel outputs against the BigFloat references, on grids that contain
# near-coincident ordinates (δ small) and exactly coincident ones.
function check_kernels(FT, AT)
    backend = KA.get_backend(AT(zeros(FT, 1)))
    δ = FT == Float64 ? 1e-10 : 2e-6
    μh = FT[0.3, 0.3 + δ, 0.75, 0.75 - δ]               # one Stokes component
    μ₀ = FT(0.75)                                        # = node 3, ≈ node 4
    N = length(μh); ns = 3
    τh = FT[1e-4, 0.03, 1.5]
    ϖh = FT[0.9, 0.5, 0.99]
    wh = FT[0.2, 0.1, 0.3, 0.0]                          # node 4 is a zero-weight view
    Zh = fill(FT(0.4), N, N, 1)
    F₀h = FT[1.2 0.7 1.0]
    above = FT[0.1, 0.0, 0.5]
    D = AT(Diagonal(ones(FT, N)))
    μ = AT(μh); τ = AT(τh); ϖ = AT(ϖh); w = AT(wh); Z = AT(Zh); F₀ = AT(F₀h)
    wct02 = FT(0.5)

    # elastic forward R/T
    r = AT(zeros(FT, N, N, ns)); t = AT(zeros(FT, N, N, ns))
    C.get_elem_rt!(backend)(r, t, ϖ, τ, Z, Z, μ, w; ndrange = size(r))
    th = Array(t)
    # elastic SFI (μ₀ = node 3 for the Z column, scalar μ₀ for geometry)
    Jp = AT(zeros(FT, N, 1, ns)); Jm = AT(zeros(FT, N, 1, ns))
    C.get_elem_rt_SFI!(backend)(Jp, Jm, ϖ, τ, AT(above), Z, Z, F₀, μ, μ₀, 0, wct02,
                                1, 3, D; ndrange = size(Jp))
    Jph = Array(Jp)
    # elastic external-solar columns
    cols = ntuple(_ -> AT(zeros(FT, N, 1, ns)), 4)
    Z0 = AT(fill(FT(0.4), N, 1, 1))
    C.get_elem_rt_solar_columns!(backend)(cols..., ϖ, τ, Z0, Z0, μ, μ₀, wct02, 1,
                                          AT(FT[1]); ndrange = size(cols[1]))
    T0h = Array(cols[3])
    # linearized (fused) R/T and SFI: one tangent direction along τ
    np = 1
    dτ̇ = AT(ones(FT, ns, np)); ϖ̇ = AT(zeros(FT, ns, np)); Ż = AT(zeros(FT, N, N, 1, np))
    lin = (AT(zeros(FT, N, N, ns)), AT(zeros(FT, N, N, ns)),
           AT(zeros(FT, N, N, ns, 3)), AT(zeros(FT, N, N, ns, 3)),
           ntuple(_ -> AT(zeros(FT, N, N, ns, np)), 4)...)
    C.get_elem_rt_fused!(backend)(lin..., ϖ, τ, Z, Z, dτ̇, ϖ̇, Ż, Ż, μ, w, np, 0, 1;
                                  ndrange = (N, N, ns))
    tl = Array(lin[2]); ṫl = Array(lin[4])
    sfl = (AT(zeros(FT, N, 1, ns)), AT(zeros(FT, N, 1, ns)),
           AT(zeros(FT, N, 1, ns, 3)), AT(zeros(FT, N, 1, ns, 3)),
           AT(zeros(FT, N, 1, ns, np)), AT(zeros(FT, N, 1, ns, np)))
    C.get_elem_rt_SFI_fused!(backend)(sfl..., ϖ, τ, AT(above), AT(zeros(FT, ns, np)),
                                      Z, Z, F₀, dτ̇, ϖ̇, Ż, Ż, μ, 0, wct02, 1,
                                      AT(FT[1]), μ₀, 3, D, np; ndrange = (N, 1, ns))
    Jpl = Array(sfl[1]); J̇pl = Array(sfl[3])
    # legacy core-partials kernels (μ₀ is the ordinate of node 3 there)
    part = (AT(zeros(FT, N, N, ns)), AT(zeros(FT, N, N, ns)),
            AT(zeros(FT, N, N, ns, 3)), AT(zeros(FT, N, N, ns, 3)))
    C.get_elem_rt!(backend)(part..., ϖ, τ, Z, Z, μ, w; ndrange = (N, N, ns))
    tp = Array(part[2]); ṫp = Array(part[4])
    sfp = (AT(zeros(FT, N, 1, ns)), AT(zeros(FT, N, 1, ns)),
           AT(zeros(FT, N, 1, ns, 3)), AT(zeros(FT, N, 1, ns, 3)))
    C.get_elem_rt_SFI!(backend)(sfp..., ϖ, τ, AT(above), AT(zeros(FT, ns, np)), Z, Z, F₀,
                                μ, 0, wct02, 1, AT(FT[1]), 3, D; ndrange = (N, 1, ns))
    Jpp = Array(sfp[1]); J̇pp = Array(sfp[3])
    # external-solar column tangents
    colsl = ntuple(_ -> AT(zeros(FT, N, 1, ns, np)), 4)
    C.get_elem_rt_solar_columns_lin!(backend)(colsl..., ϖ, τ, Z0, Z0, dτ̇, ϖ̇,
        AT(zeros(FT, N, 1, 1, np)), AT(zeros(FT, N, 1, 1, np)), μ, μ₀, wct02, 1,
        AT(FT[1]); ndrange = (N, 1, ns, np))
    Ṫ0h = Array(colsl[3])

    for n in 1:ns, i in 1:N
        τn = τh[n]; zf = FT(0.4) * F₀h[1, n]
        beam = exp(-big(above[n]) / big(μ₀))
        rf = ref_t(τn, μh[i], μ₀); rdf = ref_dt(τn, μh[i], μ₀)
        tl_ = tol(FT, τn / min(μh[i], μ₀))
        @test relerr(Jph[i, 1, n], wct02 * ϖh[n] * zf * rf * beam) <= tl_
        @test relerr(Jpl[i, 1, n], wct02 * ϖh[n] * zf * rf * beam) <= tl_
        @test relerr(Jpp[i, 1, n], wct02 * ϖh[n] * zf * rf * beam) <= tl_
        @test relerr(T0h[i, 1, n], wct02 * ϖh[n] * FT(0.4) * rf) <= tl_
        # tangents: errors on the scale of the two derivative terms
        a = big(τn) / big(μh[i]); b = big(τn) / big(μ₀)
        sc = (exp(-max(a, b)) + min(a, b) * bigDD(a, b)) / big(μh[i])
        @test abs(big(J̇pl[i, 1, n, 1]) - wct02 * ϖh[n] * zf * rdf * beam) <=
              tl_ * wct02 * ϖh[n] * zf * sc * beam
        @test abs(big(J̇pp[i, 1, n, 1]) - wct02 * ϖh[n] * zf * rdf * beam) <=
              tl_ * wct02 * ϖh[n] * zf * sc * beam
        @test abs(big(Ṫ0h[i, 1, n, 1]) - wct02 * ϖh[n] * FT(0.4) * rdf) <=
              tl_ * wct02 * ϖh[n] * FT(0.4) * sc
        for j in 1:N
            (i == j || wh[j] == 0) && continue
            rt = ref_t(τn, μh[i], μh[j]); rdt = ref_dt(τn, μh[i], μh[j])
            k = ϖh[n] * FT(0.4) * wh[j]
            tlj = tol(FT, τn / min(μh[i], μh[j]))
            for tt in (th, tl, tp)
                @test relerr(tt[i, j, n], k * rt) <= tlj
            end
            a = big(τn) / big(μh[i]); b = big(τn) / big(μh[j])
            sc = k * (exp(-max(a, b)) + min(a, b) * bigDD(a, b)) / big(μh[i])
            @test abs(big(ṫl[i, j, n, 1]) - k * rdt) <= tlj * sc
            @test abs(big(ṫp[i, j, n, 1]) - k * rdt) <= tlj * sc
        end
    end

    # canopy kernels (path factors G; μᵢGⱼ ≈ μⱼGᵢ for the near pairs, and
    # μᵢGⱼ = μⱼGᵢ for every pair with horizontal leaves, G ∝ μ, where the
    # former expression was 0/0)
    for G in (FT[0.5, 0.5, 0.6, 0.6], FT(0.9) .* μh)
        t = AT(zeros(FT, N, N, ns)); r = AT(zeros(FT, N, N, ns))
        C.get_canopy_elem_rt!(backend)(r, t, ϖ, τ, AT(G), Z, Z, μ, w; ndrange = size(t))
        tc = Array(t)
        Jp = AT(zeros(FT, N, 1, ns)); Jm = AT(zeros(FT, N, 1, ns))
        C.get_canopy_elem_rt_SFI!(backend)(Jp, Jm, ϖ, τ, AT(above), AT(G), Z, Z, μ, 0, wct02,
                                           1, F₀, 3, D; ndrange = size(Jp))
        Jc = Array(Jp)
        @test all(isfinite, tc) && all(isfinite, Jc)
        for n in 1:ns, i in 1:N
            τn = big(τh[n])
            ref_c(μi, μj, Gi, Gj) = (τn / big(μi)) * bigDD(τn * big(Gi) / big(μi), τn * big(Gj) / big(μj))
            if i != 3
                v = wct02 * ϖh[n] * big(G[3]) * FT(0.4) * F₀h[1, n] *
                    ref_c(μh[i], μh[3], G[i], G[3]) * exp(-big(above[n]) / big(μh[3]))
                @test relerr(Jc[i, 1, n], v) <= tol(FT, τh[n] / min(μh[i], μh[3]))
            end
            for j in 1:N
                (μh[i] == μh[j] || wh[j] == 0) && continue
                v = ϖh[n] * big(G[j]) * FT(0.4) * wh[j] * ref_c(μh[i], μh[j], G[i], G[j])
                @test relerr(tc[i, j, n], v) <= tol(FT, τh[n] / min(μh[i], μh[j]))
            end
        end
    end

    # Raman kernels: τ₁ = dτ[n₁] (outgoing), τ₀ = dτ[n₀] (incident). The grid
    # holds exact (τ₁ = τ₀) and near (τ₁ ≈ τ₀, μᵢ/μⱼ ≈ τ₁/τ₀) coincidences.
    d = FT(0.002)
    τR = FT[d, d, d * (1 + 10eps(FT)), d * FT(0.3) / FT(0.75), d * 2]
    nR = length(τR)
    fR = FT[0.9, 0.8, 0.85, 0.7, 0.95]                   # fscattRayl per λ
    ϖR = FT[0.02, 0.03]                                  # per Raman offset
    iR = [1, -1]                                         # n₀ = n₁ + iR
    ZR = fill(FT(0.4), N, N)
    ier = AT(zeros(FT, N, N, nR, 2)); iet = AT(zeros(FT, N, N, nR, 2))
    C.get_elem_rt_RRS!(backend)(AT(fR), AT(ϖR), AT(iR), 1, ier, iet, AT(τR), AT(fR),
                                AT(ZR), AT(ZR), μ, w; ndrange = size(ier))
    ieth = Array(iet)
    ieJp = AT(zeros(FT, N, 1, nR, 2)); ieJm = AT(zeros(FT, N, 1, nR, 2))
    aboveR = AT(fill(FT(0.05), nR)); F₀R = AT(ones(FT, 1, nR))
    C.get_elem_rt_SFI_RRS!(backend)(AT(fR), AT(ϖR), AT(iR), 1, ieJp, ieJm, aboveR,
                                    AT(τR), AT(fR), AT(ZR), AT(ZR), F₀R, μ, 0, wct02, 1,
                                    AT(FT[1]), 3, D; ndrange = size(ieJp))
    ieJph = Array(ieJp)
    colsR = ntuple(_ -> AT(zeros(FT, N, 1, nR, 2)), 4)
    C.get_elem_rt_solar_columns_RRS!(backend)(AT(fR), AT(ϖR), AT(iR), colsR...,
        AT(τR), AT(fill(FT(0.4), N, 1)), AT(fill(FT(0.4), N, 1)), μ, μ₀, wct02, 1,
        AT(FT[1]); ndrange = size(colsR[1]))
    ieT0h = Array(colsR[3])
    for Δ in 1:2, n₁ in 1:nR, i in 1:N
        n₀ = n₁ + iR[Δ]
        1 <= n₀ <= nR || continue
        K = ϖR[Δ] * fR[n₀] * FT(0.4)
        tlr(μj) = tol(FT, τR[n₁] / μh[i], τR[n₀] / μj)
        @test relerr(ieJph[i, 1, n₁, Δ],
                     wct02 * K * ref_raman(τR[n₁], τR[n₀], μh[i], μh[3]) *
                     exp(-big(FT(0.05)) / big(μh[3]))) <= tlr(μh[3])
        @test relerr(ieT0h[i, 1, n₁, Δ], wct02 * K * ref_raman(τR[n₁], τR[n₀], μh[i], μ₀)) <=
              tlr(μ₀)
        for j in 1:N
            wh[j] == 0 && continue
            @test relerr(ieth[i, j, n₁, Δ], K * wh[j] * ref_raman(τR[n₁], τR[n₀], μh[i], μh[j])) <=
                  tlr(μh[j])
        end
    end

    # Vibrational Raman (`_plus` accumulation kernels and the plain VS kernel):
    # incident wavelength n₀ = 1, outgoing n₁ = i_λ₁λ₀[Δn].
    iV = [2, 3, 4]
    ϖV = FT[0.02, 0.03, 0.01]
    ier = AT(zeros(FT, N, N, nR, 1)); iet = AT(zeros(FT, N, N, nR, 1))
    C.get_elem_rt_VS!(backend)(AT(fR), AT(ϖV), AT(iV), ier, iet, AT(τR), AT(fR),
                               AT(ZR), AT(ZR), μ, w; ndrange = (N, N, length(iV)))
    ietV = Array(iet)
    ieJp = AT(zeros(FT, N, 1, nR, 1)); ieJm = AT(zeros(FT, N, 1, nR, 1))
    C.get_elem_rt_SFI_VS!(backend)(AT(fR), AT(ϖV), AT(iV), ieJp, ieJm, aboveR, AT(τR),
                                   AT(fR), AT(ZR), AT(ZR), F₀R, μ, 0, wct02, 1,
                                   AT(FT[1]), 3, D; ndrange = (N, length(iV)))
    ieJV = Array(ieJp)
    ier = AT(zeros(FT, N, N, nR, 1)); iet = AT(zeros(FT, N, N, nR, 1))
    C.get_elem_rt_VS!(backend)(AT(fR), AT(ϖV), AT(iV .- 1), 1, ier, iet, AT(τR), AT(fR),
                               AT(ZR), AT(ZR), μ, w; ndrange = (N, N, length(iV)))
    ietV0 = Array(iet)
    for (Δ, n₁) in enumerate(iV), i in 1:N
        K = ϖV[Δ] * fR[1] * FT(0.4)
        @test relerr(ieJV[i, 1, n₁, 1],
                     wct02 * K * ref_raman(τR[n₁], τR[1], μh[i], μh[3]) *
                     exp(-big(FT(0.05)) / big(μh[3]))) <= tol(FT, τR[n₁] / μh[i], τR[1] / μh[3])
        for j in 1:N
            wh[j] == 0 && continue
            v = K * wh[j] * ref_raman(τR[n₁], τR[1], μh[i], μh[j])
            @test relerr(ietV[i, j, n₁, 1], v) <= tol(FT, τR[n₁] / μh[i], τR[1] / μh[j])
            @test relerr(ietV0[i, j, n₁, 1], fR[1] * v) <= tol(FT, τR[n₁] / μh[i], τR[1] / μh[j])
        end
    end
    return nothing
end

@testset "kernels against BigFloat near coincident ordinates (i, iii)" begin
    for FT in (Float64, Float32)
        check_kernels(FT, Array)
    end
end

if get(ENV, "VSMARTMOM_JACOBIAN_GPU_TEST", "false") == "true"
    using CUDA
    CUDA.allowscalar(false)
    @testset "kernels against BigFloat near coincident ordinates (CUDA)" begin
        for FT in (Float64, Float32)
            check_kernels(FT, CUDA.CuArray)
        end
    end
end

# ---- RT level -----------------------------------------------------------------
const QS = joinpath(@__DIR__, "..", "config", "quickstart.yaml")
function qs_params(; sza, FT = "Float64", aerosol = false, τ_ref = 0.3)
    d = YAML.load_file(QS)
    rt = d["radiative_transfer"]
    rt["polarization_type"] = "Stokes_IQU()"; rt["float_type"] = FT
    if aerosol
        rt["truncation"] = "auto"; rt["greek_beta_cutoff"] = nothing
        d["scattering"] = Dict(
            "aerosols" => [Dict("τ_ref" => τ_ref, "μ" => 0.6, "σ" => 1.6, "nᵣ" => 1.45,
                                "nᵢ" => 0.001, "p₀" => 500.0, "σp" => 100.0)],
            "r_max" => 20.0, "nquad_radius" => 50, "λ_ref" => 0.770,
            "decomp_type" => "NAI2()")
    end
    g = d["geometry"]
    g["sza"] = sza; g["vza"] = [60.0, 30.0, 0.0]; g["vaz"] = [180.0, 90.0, 0.0]
    p = read_parameters(d)
    p.sza = p.float_type(sza)
    return p
end
quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Warn))

@testset "TOA continuity across the merge window (ii)" begin
    # nstreams = 3 Gauss has a weighted node at μ = 0.5 (SZA 60° on-node).
    # Within 4 eps the SZA is merged onto the node; just outside it the
    # radiance must continue smoothly (with the former expression it differed
    # from the on-node result by ≈ 8e-5 at 5 eps in Float64, max-normalized).
    for FT in (Float64, Float32)
        e = eps(FT)
        model(μ₀) = quiet(() -> model_from_parameters(
            qs_params(; sza = acosd(Float64(μ₀)), FT = string(FT))))
        on = rt_run(model(FT(0.5))).toa
        for k in (2, 4, 4.5, 6, 10, 20), s in (1, -1)
            m = model(FT(0.5 + s * k * e))               # merged iff |δ| ≤ 4 eps
            # Merged: the node 0.5 was moved onto μ₀. Outside: node and μ₀ coexist.
            k <= 2 && @test !(FT(0.5) in Array(m.quad_points.qp_μ))
            k >= 6 && @test FT(0.5) in Array(m.quad_points.qp_μ)
            @test rt_run(m).toa ≈ on rtol = (FT == Float64 ? 1e-12 : 2e-5)
        end
    end
end

@testset "Jacobian vs finite difference near μ₀ ≈ μᵢ (iii)" begin
    # Aerosol τ_ref column (layout column 2) of the analytic Jacobian against
    # a central difference, with μ₀ 1e-10 from the weighted Gauss node 0.5.
    sza = acosd(0.5 + 1e-10)
    model, lin = quiet(() -> model_from_parameters(LinMode(), qs_params(; sza, aerosol = true)))
    @test model.quad_points.μ₀ != 0.5
    R, T, dR, dT = quiet(() -> rt_run(model, lin, 1, size(lin.τ̇_abs[1], 1), 1))
    h = 1e-4
    fwd(τ) = quiet(() -> rt_run(model_from_parameters(qs_params(; sza, aerosol = true,
                                                               τ_ref = τ))).toa)
    fd = (fwd(0.3 + h) .- fwd(0.3 - h)) ./ (2h)
    @test dR[:, :, :, 2] ≈ fd rtol = 1e-6
end

end # module
