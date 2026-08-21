# Bitwise A/B validation of the tabulated Z-moment path (ZMomentTables +
# shared per-m Π lists) against the plain compute_Z_moments.
#
# The tabulated path performs THE SAME sums in THE SAME order — P/R/T come
# from one recurrence sweep instead of per-call sweeps, and the l-accumulation
# runs on static blocks instead of allocating slices — so the contract is
# bitwise equality (0.0), not a tolerance.
using Test
using vSmartMOM
using vSmartMOM.Scattering
using vSmartMOM.Scattering: compute_Z_moments, ZMomentTables, make_Π_lists, GreekCoefs
using vSmartMOM.CoreRT: Stokes_I, Stokes_IQ, Stokes_IQU, Stokes_IQUV

# Deterministic synthetic Greek coefficients of length l_max (values are
# irrelevant for the equality contract; use a smooth decaying sequence).
function synthetic_greeks(l_max)
    l = collect(1.0:l_max)
    GreekCoefs(0.3 ./ l,          # α
               1.0 ./ l,          # β
               0.05 ./ l,         # γ
               0.8 ./ l,          # δ
               0.01 ./ l,         # ϵ
               0.7 ./ l)          # ζ
end

@testset "tabulated Z moments == plain (bitwise)" begin
    μ = collect(range(0.05, 0.995; length = 10))
    for pol in (Stokes_I{Float64}(), Stokes_IQ{Float64}(), Stokes_IQU{Float64}(), Stokes_IQUV{Float64}())
        for l_max in (3, 17, 130)          # Rayleigh-like, δ-M-like, Mie-like
            gr = synthetic_greeks(l_max)
            tables = ZMomentTables(μ, l_max)
            for m in (0, 1, 2, l_max - 1, l_max, l_max + 3)  # incl. m ≥ l_max ⇒ Z ≡ 0
                Πm = make_Π_lists(pol, tables, m)
                Z⁺⁺_ref, Z⁻⁺_ref = compute_Z_moments(pol, μ, gr, m)
                Z⁺⁺_tab, Z⁻⁺_tab = compute_Z_moments(pol, μ, gr, m;
                                                     tables = tables, Π_pair = Πm)
                @test Z⁺⁺_tab == Z⁺⁺_ref
                @test Z⁻⁺_tab == Z⁻⁺_ref
            end
        end
    end
end

@testset "tables at global l_max serve shorter expansions" begin
    # A mode whose Greek expansion is SHORTER than the table's l_max must give
    # identical Z — the recurrences are upward, so P/R/T slices don't depend
    # on how far the sweep continued.
    μ = collect(range(0.05, 0.995; length = 8))
    pol = Stokes_IQU{Float64}()
    gr_short = synthetic_greeks(5)
    tables_global = ZMomentTables(μ, 200)
    for m in (0, 1, 3)
        Πm = make_Π_lists(pol, tables_global, m)
        Z⁺⁺_ref, Z⁻⁺_ref = compute_Z_moments(pol, μ, gr_short, m)
        Z⁺⁺_tab, Z⁻⁺_tab = compute_Z_moments(pol, μ, gr_short, m;
                                             tables = tables_global, Π_pair = Πm)
        @test Z⁺⁺_tab == Z⁺⁺_ref
        @test Z⁻⁺_tab == Z⁻⁺_ref
    end
end

@testset "tables without Π_pair auto-build (silent-zero regression)" begin
    # Passing `tables` while omitting `Π_pair` must NOT silently skip the
    # accumulation and return Z ≡ 0 — the method builds the Π lists itself.
    μ = collect(range(0.05, 0.995; length = 8))
    pol = Stokes_IQU{Float64}()
    gr = synthetic_greeks(17)
    tables = ZMomentTables(μ, 17)
    for m in (0, 1, 4)
        Z⁺⁺_ref, Z⁻⁺_ref = compute_Z_moments(pol, μ, gr, m)
        Z⁺⁺_tab, Z⁻⁺_tab = compute_Z_moments(pol, μ, gr, m; tables = tables)
        @test Z⁺⁺_tab == Z⁺⁺_ref            # not zeros — bitwise equal
        @test Z⁻⁺_tab == Z⁻⁺_ref
        @test any(!iszero, Z⁺⁺_tab)
    end
end

@testset "tables built for a different μ are rejected" begin
    μ  = collect(range(0.05, 0.995; length = 8))
    μ2 = collect(range(0.10, 0.90;  length = 8))   # same length, different nodes
    μ3 = collect(range(0.05, 0.995; length = 6))   # different length
    pol = Stokes_IQU{Float64}()
    gr = synthetic_greeks(17)
    tables = ZMomentTables(μ, 17)
    @test_throws ArgumentError compute_Z_moments(pol, μ2, gr, 1; tables = tables)
    @test_throws ArgumentError compute_Z_moments(pol, μ3, gr, 1; tables = tables)
end

@testset "kill-switch default" begin
    # _Z_TABLES_ENABLED=false ⇒ build_m_invariant_cache stores z_tables===nothing
    # ⇒ every call falls back to the plain method. The end-to-end radiance A/B
    # (tables on vs off, R_on == R_off bitwise) runs in the GPU benchmark
    # harness; here we pin the default so CI catches an accidental flip.
    @test vSmartMOM.CoreRT._Z_TABLES_ENABLED[] == true
end

println("test_z_tables: done")
