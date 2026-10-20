# Vibrational Raman over concatenated bands (VS_0to1_plus, VS_1to0_plus).
#
# One incident wavelength (band 1) and its N₂/O₂ vibrational target bands, pure
# Rayleigh, one layer, flat F₀. Covers the forward run itself.

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

end # module
