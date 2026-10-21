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

@testset "apply_gridlines! keeps line strength and position" begin
    λ₀ = 1e7 / 13000
    ν₀ = 1e7 / λ₀
    S(σ, Δ) = σ * (ν₀ + Δ)^4
    for FT in (Float32, Float64), grid in (collect(FT, 12990:0.5:13010),
                                          reverse(collect(FT, 12990:0.5:13010)))
        σ_out = zeros(FT, length(grid))
        # A line 0.3 of a step above a grid point.
        Δ = 3.15
        IS.apply_gridlines!([Δ], [1e-45], λ₀, grid, σ_out)
        @test count(>(0), σ_out) == 2
        @test sum(σ_out) ≈ S(1e-45, Δ) rtol = 1e-6
        @test sum(σ_out .* (Float64.(grid) .- ν₀)) / sum(σ_out) ≈ Δ atol = 2e-3
        k = findfirst(>(0), σ_out)
        i_low = grid[1] < grid[end] ? k : k + 1
        @test σ_out[i_low] / sum(σ_out) ≈ 0.7 rtol = 1e-3
        # A line on a grid point goes to that point only.
        IS.apply_gridlines!([-4.0], [1e-45], λ₀, grid, σ_out)
        @test count(>(0), σ_out) == 1
        @test σ_out[findfirst(==(FT(ν₀ - 4)), grid)] ≈ S(1e-45, -4.0) rtol = 1e-6
        # Lines on the end points are kept on them; lines outside are dropped.
        IS.apply_gridlines!([-10.0, 10.0], [1e-45, 1e-45], λ₀, grid, σ_out)
        @test count(>(0), σ_out) == 2
        @test σ_out[findfirst(==(FT(12990)), grid)] ≈ S(1e-45, -10.0) rtol = 1e-6
        @test σ_out[findfirst(==(FT(13010)), grid)] ≈ S(1e-45, 10.0) rtol = 1e-6
        IS.apply_gridlines!([-12.0, 15.0], [1e-45, 1e-45], λ₀, grid, σ_out)
        @test all(iszero, σ_out)
    end
end

@testset "RRS offsets are taken from the reference grid point" begin
    ν̃ = 25000.0
    n2, o2 = IS.getRamanAtmoConstants(ν̃, 250.0)
    σ_Rayl = IS.compute_optical_Rayl(1e7 / ν̃, n2, o2)
    h = 0.5
    odd = collect((ν̃ - 30):h:(ν̃ + 30))   # the band mean is a grid point
    even = odd[2:end]                     # the band mean lies between two
    # Reference wavenumbers on and between grid points; ν̃ is the nearest
    # grid point to both.
    results = map(((odd, ν̃), (even, ν̃ + 0.2))) do (ν, ν_ref)
        rs = _rrs(Float64, ν̃)
        CoreRT.getRamanSSProp!(rs, 1e7 / ν_ref, ν)
        (copy(rs.i_λ₁λ₀), copy(rs.ϖ_λ₁λ₀))
    end
    # The same offsets and relative weights: an off-grid reference no longer
    # moves the lines (by up to h/2). The absolute scale differs by the
    # Rayleigh cross section at the given λ, which rt_run normalizes away.
    @test results[1][1] == results[2][1]
    @test results[1][2] ./ sum(results[1][2]) ≈ results[2][2] ./ sum(results[2][2]) rtol = 1e-12

    i, w = results[1]
    @test issorted(i)
    # Every source offset has its mirror image (Stokes and anti-Stokes lines
    # of the same pair of levels have the same |shift|).
    @test sort(-i) == i
    # i > 0: the source lies above the output (Stokes, wavenumber decreases).
    # Stokes lines are stronger than their anti-Stokes partners.
    for (Δn, iΔ) in enumerate(i)
        iΔ > 0 || continue
        @test w[Δn] > w[findfirst(==(-iΔ), i)]
    end

    # The weight-averaged Stokes shift equals the strength-averaged shift of
    # the Stokes lines inside the band (linear interpolation keeps positions).
    lines = Tuple{Float64, Float64}[]
    for mol in (n2, o2), J in 0:28
        c = mol.effCoeff
        Δ, σ = c.Δν̃_RoRaman_coeff_JtoJp2[J], c.σ_RoRaman_coeff_JtoJp2[J]
        (σ > 0 && abs(Δ) < 30) && push!(lines, (Δ, mol.vmr * σ * (ν̃ + Δ)^4))
    end
    S_stokes = sum(last, lines)
    @test sum(w[i .> 0]) ≈ S_stokes / σ_Rayl rtol = 1e-10
    @test sum(w[i .> 0] .* i[i .> 0]) * h / sum(w[i .> 0]) ≈
          -sum(prod, lines) / S_stokes rtol = 1e-10

    # A reference wavelength outside the band is rejected, one on its edge is not.
    @test_throws ArgumentError CoreRT.getRamanSSProp!(_rrs(Float64, ν̃), 1e7 / (ν̃ + 40), odd)
    @test_throws ArgumentError CoreRT.getRamanSSProp!(_rrs(Float64, ν̃), 1e7 / (last(odd) + 0.2), odd)
    rs_edge = _rrs(Float64, ν̃)
    CoreRT.getRamanSSProp!(rs_edge, 1e7 / last(odd), odd)
    @test rs_edge.n_Raman > 0
    # Also with a Float32 grid and a Float32 λ built from the band edge.
    odd32 = Float32.(odd)
    for ν_edge in (first(odd32), last(odd32))
        rs32 = _rrs(Float32, ν̃)
        CoreRT.getRamanSSProp!(rs32, Float32(1e7 / ν_edge), odd32)
        @test rs32.n_Raman > 0
    end
end

@testset "Unused getRamanSSProp! overloads are removed" begin
    # R6: the 4-argument RRS method reversed both ϖ_λ₁λ₀ and i_λ₁λ₀, which
    # swapped Stokes and anti-Stokes weights. R8: the VS_0to1/VS_1to0 method
    # assigned fields that the immutable VS structs do not have. Neither was
    # called; the RRS entry point is getRamanSSProp!(rs, λ, ν).
    @test !hasmethod(CoreRT.getRamanSSProp!, Tuple{IS.RRS, Any, Any, Any})
    @test !hasmethod(CoreRT.getRamanSSProp!, Tuple{IS.VS_0to1, Any, Any, Any})
    @test !hasmethod(CoreRT.getRamanSSProp!, Tuple{IS.VS_1to0, Any, Any, Any})
    @test hasmethod(CoreRT.getRamanSSProp!, Tuple{IS.RRS, Any, Any})
end

end # module
