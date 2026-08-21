# Regression for the sparse (layer-resolved) aerosol combine in
# constructCoreOpticalProperties — the _SPARSE_AERO_COMBINE_ENABLED path.
#
# Two regimes:
#   dense  — a mode spans >1 layer ⇒ _detect_mode_layers returns `nothing`
#            and the sparse switch must be a no-op (identical code path).
#   sparse — every mode occupies exactly one layer (the satellite campaign's
#            layer-resolved diagonal τ) ⇒ mode i is combined ONLY into layer
#            mode_layers[i], skipping the O(nAero×nZ) zero-τ adds.
#
# Contract (kill-switch docstring): the sparse path skips the dense path's
# zero-τ round-trips `(ϖ·τ + 0)/τ`, making it tolerance-equivalent and
# marginally MORE exact; in every measured case the results are bitwise
# equal. We assert bitwise here and would want to know if that ever drifts.
using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.CoreRT: constructCoreOpticalProperties, build_m_invariant_cache,
                        _SPARSE_AERO_COMBINE_ENABLED, _detect_mode_layers
using vSmartMOM.InelasticScattering: noRS

params = vSmartMOM.parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
model  = model_from_parameters(params)
FT     = vSmartMOM.CoreRT.float_type(model)
RS     = noRS()   # zero-arg pure-elastic default (ϖ_Cabannes = 1 per band)

@testset "mode-layer detection" begin
    τ = copy(model.τ_aer[1])                    # [iAer, nSpec, iz]
    @test _detect_mode_layers(τ) === nothing || _detect_mode_layers(τ) isa Vector{Int}
    # force the layer-resolved diagonal shape: mode 1 only in layer 3
    τ .= 0; τ[1, :, 3] .= FT(0.01)
    @test _detect_mode_layers(τ) == [3]
    # a mode spanning two layers ⇒ dense combine required
    τ[1, :, 5] .= FT(0.01)
    @test _detect_mode_layers(τ) === nothing
    # an all-zero mode is flagged 0, not treated as spanning
    τ .= 0
    @test _detect_mode_layers(τ) == [0]
end

"Run the full m-loop's constructCoreOpticalProperties A/B and compare bitwise."
function ab_layers(model, m)
    cache_on  = (_SPARSE_AERO_COMBINE_ENABLED[] = true;  build_m_invariant_cache(RS, [1], model))
    lp_on,  _ = constructCoreOpticalProperties(RS, [1], m, model, cache_on)
    cache_off = (_SPARSE_AERO_COMBINE_ENABLED[] = false; build_m_invariant_cache(RS, [1], model))
    lp_off, _ = constructCoreOpticalProperties(RS, [1], m, model, cache_off)
    _SPARSE_AERO_COMBINE_ENABLED[] = true
    same = true
    for (a, b) in zip(lp_on, lp_off)
        same &= Array(a.τ) == Array(b.τ) && Array(a.ϖ) == Array(b.ϖ) &&
                Array(a.Z⁺⁺) == Array(b.Z⁺⁺) && Array(a.Z⁻⁺) == Array(b.Z⁻⁺)
    end
    return same
end

@testset "sparse == dense (bitwise) — dense regime (mode spans all layers)" begin
    # JacobianTestFast's aerosol profile spans many layers ⇒ ml === nothing on
    # both sides; the switch must change nothing at all.
    for m in (0, 1)
        @test ab_layers(model, m)
    end
end

@testset "sparse == dense (bitwise) — layer-resolved diagonal regime" begin
    # Rewrite τ_aer into the satellite campaign's shape: the single mode
    # occupies exactly one layer, so the sparse path activates.
    τ = model.τ_aer[1]
    τ .= 0
    τ[1, :, 4] .= FT(0.02)
    @test _detect_mode_layers(τ) == [4]
    for m in (0, 1, 2)
        @test ab_layers(model, m)
    end
end

println("test_sparse_combine: done")
