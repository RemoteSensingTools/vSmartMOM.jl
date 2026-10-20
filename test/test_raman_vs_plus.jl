# Vibrational Raman over concatenated bands (VS_0to1_plus, VS_1to0_plus).
#
# One incident wavelength (band 1) and its N₂/O₂ vibrational target bands, pure
# Rayleigh, one layer, flat F₀. Covers the forward run, the source-band
# depletion (decision 7, R5) and the model builder (decision 8).

module RamanVSPlusTests

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Test

const IS = InelasticScattering
const FT = Float64
const λ_INC = 395.0   # nm

function _vs_params()
    p = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "O2ABandParamsVS.yaml"))
    p.absorption_params = nothing
    p.scattering_params = nothing
    p.vza = [30.0]
    p.vaz = [0.0]
    return p
end

_g1() = IS.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])

_vs(::Type{IS.VS_0to1_plus}, n2, o2) = IS.VS_0to1_plus{FT}(n2 = n2, o2 = o2)
_vs(::Type{IS.VS_1to0_plus}, n2, o2) = IS.VS_1to0_plus{FT}(
    bandSpecLim = UnitRange{Int}[], iBand = Int[], grid_in = StepRangeLen{FT}[],
    n2 = n2, o2 = o2,
    greek_raman = _g1(), greek_raman_VS_n2 = _g1(), greek_raman_VS_o2 = _g1(),
    fscattRayl = zeros(FT, 1), ϖ_Cabannes = zeros(FT, 1),
    ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
    ϖ_λ₁λ₀_VS_n2 = zeros(FT, 1), i_λ₁λ₀_VS_n2 = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀_VS_n2 = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀_VS_n2 = zeros(FT, 1, 1),
    ϖ_λ₁λ₀_VS_o2 = zeros(FT, 1), i_λ₁λ₀_VS_o2 = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀_VS_o2 = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀_VS_o2 = zeros(FT, 1, 1),
    i_λ₁λ₀_all = zeros(Int, 1), i_ref = 1, n_Raman = 1,
    F₀ = zeros(FT, 1, 1), SIF₀ = zeros(FT, 1, 1))

"Build a VS plus model; return `(rs, model, params)`."
function _vs_model(::Type{T}; params = _vs_params()) where {T}
    n2, o2 = IS.getRamanAtmoConstants(1e7 / λ_INC, 300.0)
    rs = _vs(T, n2, o2)
    model = model_from_parameters(rs, λ_INC, params)
    return rs, model, params
end

function _flat_beam!(rs, model)
    nspec = sum(length, rs.grid_in)
    rs.F₀ = zeros(FT, CoreRT.polarization_type(model).n, nspec)
    rs.F₀[1, :] .= 1
    return nspec
end

@testset "VS_0to1_plus forward run" begin
    rs, model, _ = _vs_model(IS.VS_0to1_plus)
    nspec = _flat_beam!(rs, model)
    R, T, ieR, ieT = CoreRT.rt_run(rs, model, rs.iBand)
    @test size(R) == (1, 3, nspec)
    @test all(isfinite, R) && all(isfinite, ieR)
    @test all(isfinite, T) && all(isfinite, ieT)
    @test R[1, 1, 1] > 0
    # Only the target wavelengths of the Raman lines receive inelastic light.
    active = sort(unique(filter(>(0), rs.i_λ₁λ₀_all)))
    @test all(>(0), ieR[1, 1, active])
    @test all(iszero, ieR[1, 1, setdiff(1:nspec, active)])
end

@testset "VS plus: the incident band is the Raman source" begin
    # Both directions take the vibrational Raman light from band 1 (the
    # incident wavelength), so band 1 is depleted by the VRS + RVRS fraction
    # and the target bands keep ϖ_Cabannes = 1. The depletion covers both
    # directions (0→1 and 1→0); the weights of a VS_0to1_plus run carry the
    # 0→1 part, i.e. all but ~6e-5 of it at the profile's effective
    # temperature (the builder recomputes the molecular constants there).
    for T in (IS.VS_0to1_plus, IS.VS_1to0_plus)
        rs, _, _ = _vs_model(T)
        @test rs.ϖ_Cabannes[1] < 1
        @test rs.ϖ_Cabannes[2:end] == ones(FT, length(rs.ϖ_Cabannes) - 1)
    end
    rs, _, _ = _vs_model(IS.VS_0to1_plus)
    depletion = 1 - rs.ϖ_Cabannes[1]
    raman = sum(rs.ϖ_λ₁λ₀) + sum(rs.ϖ_λ₁λ₀_VS_n2) + sum(rs.ϖ_λ₁λ₀_VS_o2)
    @test 0 < depletion - raman < 1e-4 * depletion
end

@testset "VS plus builder: resolved depolarization, caller's params untouched" begin
    # Decision 8: the builder uses the depolarization resolved as in the
    # standard model_from_parameters (depol: -1 → molecular constants) for
    # τ_Rayleigh and the Rayleigh/Cabannes Greek coefficients, and works on a
    # copy of the parameters.
    for depol in (-1.0, 0.0)
        p = _vs_params()
        p.spec_bands = p.spec_bands[1:1]
        p.brdf = p.brdf[1:1]
        p.depol = depol
        spec_before = deepcopy(p.spec_bands)
        brdf_before = p.brdf
        rs, model, _ = _vs_model(IS.VS_0to1_plus; params = p)
        @test p.spec_bands == spec_before
        @test p.brdf === brdf_before && length(p.brdf) == 1
        @test length(CoreRT.get_spec_bands(model)) == 3

        # Band 1 is the incident wavelength alone; the standard builder on that
        # one point gives the reference Rayleigh optics.
        q = _vs_params()
        q.spec_bands = [collect(rs.grid_in[1])]
        q.brdf = q.brdf[1:1]
        q.depol = depol
        ref = model_from_parameters(q)
        @test model.τ_rayl[1] ≈ ref.τ_rayl[1] rtol = 1e-12
        @test model.greek_rayleigh[1].β ≈ ref.greek_rayleigh[1].β rtol = 1e-12
        @test model.greek_rayleigh[1].δ ≈ ref.greek_rayleigh[1].δ rtol = 1e-12
        @test model.greek_cabannes[1].β ≈ ref.greek_cabannes[1].β rtol = 1e-12
        @test model.ϖ_Cabannes[1] ≈ ref.ϖ_Cabannes[1] rtol = 1e-12
    end

    # Float32 model (builder only).
    p = _vs_params()
    p.float_type = Float32
    n2, o2 = IS.getRamanAtmoConstants(1e7 / λ_INC, 300.0)
    rs32 = IS.VS_0to1_plus{Float32}(n2 = n2, o2 = o2)
    model32 = model_from_parameters(rs32, λ_INC, p)
    @test eltype(model32.τ_rayl[1]) == Float32
    @test model32.greek_rayleigh[1] isa vSmartMOM.Scattering.GreekCoefs{Float32}
    @test model32.greek_cabannes[1] isa vSmartMOM.Scattering.GreekCoefs{Float32}
    @test model32.τ_rayl[1] ≈ model_from_parameters(
        IS.VS_0to1_plus{FT}(n2 = n2, o2 = o2), λ_INC, _vs_params()).τ_rayl[1] rtol = 1e-5

    # With gas absorption (O₂ line-by-line from the HITRAN artifact): the
    # builder reads the lines in the model float type.
    p = _vs_params()
    p_abs = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "O2ABandParamsVS.yaml"))
    p.absorption_params = p_abs.absorption_params
    n2, o2 = IS.getRamanAtmoConstants(1e7 / 765.0, 300.0)
    rs = IS.VS_0to1_plus{FT}(n2 = n2, o2 = o2)
    model = model_from_parameters(rs, 765.0, p)
    @test all(isfinite, model.τ_abs[1]) && sum(model.τ_abs[1]) > 0
end

end # module
