# Run from test/, using --project=../docs. These execute the maintained Literate
# sources, including plots; the documentation renderer alone does not execute them.
using Test, vSmartMOM

const tutorial_names = isempty(ARGS) ?
    ["QuickStart", "IO", "CoreRT", "Jacobians", "HybridAD",
     "Absorption", "Scattering", "Surfaces", "Canopy", "MieDeepDive"] : ARGS
const tutorial_dir = joinpath(pkgdir(vSmartMOM), "docs", "src", "pages", "tutorials")

@testset "Executable tutorials" begin
    for name in tutorial_names
        @testset "$name" begin
            path = joinpath(tutorial_dir, "Tutorial_$(name).jl")
            example = Module(gensym(Symbol(name)))
            Base.include(example, path)
            # Validate the outputs that these examples promise, not just parsing.
            for field in (:R, :T, :dR, :dT, :R_cpu, :dR_cpu, :R_canopy)
                if isdefined(example, field)
                    @test all(isfinite, Base.invokelatest(getfield, example, field))
                end
            end
            if isdefined(example, :result)
                result = Base.invokelatest(getfield, example, :result)
                @test size(result.toa_jacobian, 4) == vSmartMOM.CoreRT.n_total(result.layout)
            end
            @test true # completion, including tutorials without RT output arrays
        end
    end
end
