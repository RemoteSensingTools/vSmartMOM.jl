# Rotational Raman (RRS) line weights on the band grid.
#
# A band that holds no rotational line is rejected (decision 4 of the v2.2
# decision log): the elastic part would still be reduced to the Cabannes
# fraction, so the Raman fraction 1 − ϖ_Cabannes (≈ 3 %) would be lost.

module RamanRRSLinesTests

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Test

const IS = InelasticScattering

function _rrs(::Type{FT}, ν̃, T = 250.0) where {FT}
    n2, o2 = IS.getRamanAtmoConstants(ν̃, T)
    return IS.RRS(n2 = n2, o2 = o2,
        greek_raman = IS.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl = [FT(1)], ϖ_Cabannes = [FT(1)],
        ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = 1, n_Raman = 0,
        F₀ = zeros(FT, 1, 1), SIF₀ = zeros(FT, 1, 1))
end

# The two RRS calls only read `ϖ_Cabannes` and the band grid from the model.
_band_model(ν, ϖ_C) = (ϖ_Cabannes = [ϖ_C], atmosphere = (spec_bands = [ν],))

@testset "RRS band without rotational lines is rejected" begin
    for FT in (Float32, Float64)
        # 11 points around 25000 cm⁻¹: ±2.5 cm⁻¹, inside the smallest shifts
        # (≈ ±11.9 cm⁻¹ for N₂, ±14.3 cm⁻¹ for O₂).
        ν = collect(FT, (25000 - 2.5):0.5:(25000 + 2.5))
        rs = _rrs(FT, 25000.0)
        err = try
            CoreRT.getRamanSSProp!(rs, 1e7 / 25000, ν)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("no N₂/O₂ rotational Raman line", msg)
        @test occursin("noRS()", msg)
        @test occursin("±11.9 cm⁻¹ (N₂)", msg) && occursin("±14.3 cm⁻¹ (O₂)", msg)

        # Weights that sum to zero (e.g. the placeholder vector of an RRS
        # object that never went through getRamanSSProp!) are rejected by the
        # normalization every RRS solver applies before its first moment.
        @test_throws ArgumentError IS.normalize_raman_weights!(
            _rrs(FT, 25000.0), _band_model(ν, FT(0.968)), 1)
    end

    # ±25 cm⁻¹ holds the lowest N₂/O₂ lines: no error, and the normalized
    # weights carry exactly the Raman fraction.
    ν = collect((25000 - 25.0):0.5:(25000 + 25.0))
    rs = _rrs(Float64, 25000.0)
    CoreRT.getRamanSSProp!(rs, 1e7 / 25000, ν)
    @test rs.n_Raman > 0
    ϖ_C = rs.ϖ_Cabannes[1]
    IS.normalize_raman_weights!(rs, _band_model(ν, ϖ_C), 1)
    @test sum(rs.ϖ_λ₁λ₀) ≈ 1 - ϖ_C rtol = 1e-12
end

end # module
