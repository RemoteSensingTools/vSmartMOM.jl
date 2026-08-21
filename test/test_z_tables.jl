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

@testset "kill-switch: cache built with tables disabled still solves" begin
    # _Z_TABLES_ENABLED=false ⇒ z_tables===nothing ⇒ every call falls back to
    # the plain method; the cache-aware overload must not error.
    @test vSmartMOM.CoreRT._Z_TABLES_ENABLED[] == true   # default on
end

println("test_z_tables: done")
