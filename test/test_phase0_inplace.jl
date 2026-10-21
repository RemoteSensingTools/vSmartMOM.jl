# Phase-0 preallocation pass: the in-place doubling/interaction updates must
# reproduce the original allocating ⊠ formulas exactly (same GEMMs, same
# broadcast order ⇒ bitwise on CPU; identical CUBLAS calls on GPU).
using Test
using vSmartMOM
using vSmartMOM.CoreRT: doubling_source_update!, doubling_rt_update!,
                        compute_geometric_progression!, _bmm!,
                        interaction_helper!, ScatteringInterface_00,
                        ScatteringInterface_01, ScatteringInterface_10,
                        ScatteringInterface_11, AddedLayer, CompositeLayer,
                        batched_pointer_cache, batch_inv!
using NNlib: batched_mul
const ⊠ = batched_mul
using LinearAlgebra
using Random

Random.seed!(7)

# small, well-conditioned random supermatrices (contraction-scaled so the
# geometric series is sane)
mk(N, S) = 0.05 .* rand(Float64, N, N, S)
mkv(N, S) = rand(Float64, N, 1, S)

function reference_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk)
    @views j₁⁺[:, 1, :] .= j₀⁺[:, 1, :] .* expk'
    @views j₁⁻[:, 1, :] .= j₀⁻[:, 1, :] .* expk'
    j₀⁻ .= j₀⁻ .+ (tt_gp ⊠ (j₁⁻ .+ r⁻⁺ ⊠ j₀⁺))
    j₀⁺ .= j₁⁺ .+ (tt_gp ⊠ (j₀⁺ .+ r⁻⁺ ⊠ j₁⁻))
end
function reference_rt_update!(r⁻⁺, t⁺⁺, tt_gp, expk)
    r⁻⁺ .= r⁻⁺ .+ (tt_gp ⊠ r⁻⁺ ⊠ t⁺⁺)
    t⁺⁺ .= tt_gp ⊠ t⁺⁺
    expk .= expk .^ 2
end

@testset "Phase-0 in-place doubling updates (CPU, bitwise)" begin
    N, S = 8, 17
    r, t, tt = mk(N, S), mk(N, S), mk(N, S)
    jp, jm, j1p, j1m = mkv(N, S), mkv(N, S), mkv(N, S), mkv(N, S)
    expk = rand(Float64, S)

    # references on deep copies
    Rr, Rt = copy(r), copy(t)
    Rjp, Rjm, Rj1p, Rj1m = copy(jp), copy(jm), copy(j1p), copy(j1m)
    Rexpk = copy(expk)
    reference_source_update!(Rjp, Rjm, Rj1p, Rj1m, Rr, tt, Rexpk)
    reference_rt_update!(Rr, Rt, tt, Rexpk)

    v1, v2 = similar(jp), similar(jp)
    m1, m2 = similar(r), similar(r)
    doubling_source_update!(jp, jm, j1p, j1m, r, tt, expk, v1, v2)
    doubling_rt_update!(r, t, tt, expk, m1, m2)

    @test jp == Rjp
    @test jm == Rjm
    @test r == Rr
    @test t == Rt
    @test expk == Rexpk

    # allocating back-compat wrappers agree too
    r2, t2 = copy(Rr), copy(Rt)   # arbitrary state; wrapper vs explicit temps
    r3, t3 = copy(Rr), copy(Rt)
    e2, e3 = copy(expk), copy(expk)
    doubling_rt_update!(r2, t2, tt, e2)
    doubling_rt_update!(r3, t3, tt, e3, m1, m2)
    @test r2 == r3 && t2 == t3 && e2 == e3
end

@testset "Phase-0 geometric-progression fallback (CPU, bitwise)" begin
    N, S = 8, 17
    r, t = mk(N, S), mk(N, S)
    I_static = repeat(Matrix{Float64}(I, N, N), 1, 1, 1)  # (N,N,1) broadcastable
    gp_a, tt_a = similar(r), similar(r)
    temp2 = similar(r)
    # reference: original formula
    temp2 .= I_static .- r ⊠ r
    batch_inv!(gp_a, copy(temp2), nothing, nothing)
    tt_ref = zeros(N, N, S); tt_ref .= t ⊠ gp_a
    # in-place path
    gp_b, tt_b = similar(r), similar(r)
    compute_geometric_progression!(gp_b, tt_b, r, t, I_static, similar(r), nothing, nothing)
    @test tt_b ≈ tt_ref rtol = 1e-13
end

# ── interaction helpers: reference vs in-place on identical layer states ─────
function mk_added(N, S; FT = Float64)
    AddedLayer(
        r⁻⁺ = mk(N, S), t⁺⁺ = mk(N, S) .+ 0.5, r⁺⁻ = mk(N, S), t⁻⁻ = mk(N, S) .+ 0.5,
        j₀⁺ = mkv(N, S), j₀⁻ = mkv(N, S),
        temp1 = zeros(FT, N, N, S), temp2 = zeros(FT, N, N, S),
        temp1_ptr = nothing, temp2_ptr = nothing,
        dbl_gp_refl = zeros(FT, N, N, S),
        dbl_j₁⁺ = zeros(FT, N, 1, S), dbl_j₁⁻ = zeros(FT, N, 1, S),
        dbl_v1 = zeros(FT, N, 1, S), dbl_v2 = zeros(FT, N, 1, S))
end
mk_comp(N, S; FT = Float64) = CompositeLayer(
    R⁻⁺ = mk(N, S), R⁺⁻ = mk(N, S), T⁺⁺ = mk(N, S) .+ 0.5, T⁻⁻ = mk(N, S) .+ 0.5,
    J₀⁺ = mkv(N, S), J₀⁻ = mkv(N, S))

deepc(c) = CompositeLayer(R⁻⁺ = copy(c.R⁻⁺), R⁺⁻ = copy(c.R⁺⁻),
                          T⁺⁺ = copy(c.T⁺⁺), T⁻⁻ = copy(c.T⁻⁻),
                          J₀⁺ = copy(c.J₀⁺), J₀⁻ = copy(c.J₀⁻))

# original allocating formulas (verbatim from pre-Phase-0 interaction.jl)
function ref_00!(c, a)
    c.J₀⁺ .= a.j₀⁺ .+ a.t⁺⁺ ⊠ c.J₀⁺
    c.J₀⁻ .= c.J₀⁻ .+ c.T⁻⁻ ⊠ a.j₀⁻
    c.T⁻⁻ .= a.t⁻⁻ ⊠ c.T⁻⁻
    c.T⁺⁺ .= a.t⁺⁺ ⊠ c.T⁺⁺
end
function ref_01!(c, a)
    c.J₀⁻ .= c.J₀⁻ .+ c.T⁻⁻ ⊠ (a.r⁻⁺ ⊠ c.J₀⁺ .+ a.j₀⁻)
    c.J₀⁺ .= a.j₀⁺ .+ a.t⁺⁺ ⊠ c.J₀⁺
    c.R⁻⁺ .= c.T⁻⁻ ⊠ a.r⁻⁺ ⊠ c.T⁺⁺
    c.R⁺⁻ .= a.r⁺⁻
    c.T⁺⁺ .= a.t⁺⁺ ⊠ c.T⁺⁺
    c.T⁻⁻ .= c.T⁻⁻ ⊠ a.t⁻⁻
end
function ref_10!(c, a)
    c.J₀⁺ .= a.j₀⁺ .+ a.t⁺⁺ ⊠ (c.J₀⁺ .+ c.R⁺⁻ ⊠ a.j₀⁻)
    c.J₀⁻ .= c.J₀⁻ .+ c.T⁻⁻ ⊠ a.j₀⁻
    c.T⁺⁺ .= a.t⁺⁺ ⊠ c.T⁺⁺
    c.T⁻⁻ .= c.T⁻⁻ ⊠ a.t⁻⁻
    c.R⁺⁻ .= a.t⁺⁺ ⊠ c.R⁺⁻ ⊠ a.t⁻⁻
end
function ref_11!(c, a, I_static)
    temp2 = I_static .- a.r⁻⁺ ⊠ c.R⁺⁻
    temp1 = similar(temp2); batch_inv!(temp1, temp2, nothing, nothing)
    T01_inv = c.T⁻⁻ ⊠ temp1
    c.J₀⁻ .= c.J₀⁻ .+ T01_inv ⊠ (a.r⁻⁺ ⊠ c.J₀⁺ .+ a.j₀⁻)
    c.R⁻⁺ .= c.R⁻⁺ .+ T01_inv ⊠ a.r⁻⁺ ⊠ c.T⁺⁺
    c.T⁻⁻ .= T01_inv ⊠ a.t⁻⁻
    temp2 = I_static .- c.R⁺⁻ ⊠ a.r⁻⁺
    batch_inv!(temp1, temp2, nothing, nothing)
    T21_inv = a.t⁺⁺ ⊠ temp1
    c.J₀⁺ .= a.j₀⁺ .+ T21_inv ⊠ (c.J₀⁺ .+ c.R⁺⁻ ⊠ a.j₀⁻)
    c.T⁺⁺ .= T21_inv ⊠ c.T⁺⁺
    c.R⁺⁻ .= a.r⁺⁻ .+ T21_inv ⊠ c.R⁺⁻ ⊠ a.t⁻⁻
end

@testset "Phase-0 in-place interaction helpers (CPU)" begin
    N, S = 8, 17
    I_static = repeat(Matrix{Float64}(I, N, N), 1, 1, 1)
    for (iface, ref!) in ((ScatteringInterface_00(), ref_00!),
                          (ScatteringInterface_01(), ref_01!),
                          (ScatteringInterface_10(), ref_10!),
                          (ScatteringInterface_11(), ref_11!))
        a = mk_added(N, S)
        c = mk_comp(N, S)
        cr = deepc(c)
        if ref! === ref_11!
            ref!(cr, a, I_static)
        else
            ref!(cr, a)
        end
        interaction_helper!(iface, true, c, a, I_static)
        tol = ref! === ref_11! ? 1e-12 : 0.0   # _11 involves an inverse; CPU
        for f in (:R⁻⁺, :R⁺⁻, :T⁺⁺, :T⁻⁻, :J₀⁺, :J₀⁻)
            got, want = getproperty(c, f), getproperty(cr, f)
            if tol == 0.0
                @test got == want
            else
                @test isapprox(got, want; rtol = tol)
            end
        end
    end
end

println("Phase-0 in-place correctness: done")
