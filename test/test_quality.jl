using Aqua
using JSON
using TOML
using Test
using vSmartMOM

@testset "Taplo schema config" begin
    root = normpath(joinpath(@__DIR__, ".."))
    taplo = TOML.parsefile(joinpath(root, ".taplo.toml"))

    @test haskey(taplo, "rule")

    schema_dir = joinpath(root, "schemas")
    schema_names = [
        "julia-artifacts.schema.json",
        "julia-project.schema.json",
        "vsmartmom-parameters.schema.json",
    ]
    schema_paths = Set(joinpath(schema_dir, name) for name in schema_names)

    for schema_path in schema_paths
        schema = JSON.parsefile(schema_path)
        @test schema["\$schema"] == "http://json-schema.org/draft-04/schema#"
        @test schema["type"] == "object"
        @test !isempty(schema["title"])
        @test haskey(schema, "x-taplo")
    end

    for rule in taplo["rule"]
        schema = get(rule, "schema", Dict{String,Any}())
        schema_path = get(schema, "path", nothing)
        schema_path isa String || continue
        startswith(schema_path, "taplo://") && continue

        local_schema_path = normpath(joinpath(root, schema_path))
        @test local_schema_path in schema_paths
        @test isfile(local_schema_path)
    end

    for toml_path in ["Project.toml", "test/Project.toml", "docs/Project.toml", "Artifacts.toml", ".taplo.toml"]
        @test TOML.parsefile(joinpath(root, toml_path)) isa Dict
    end

    parameter_schema = JSON.parsefile(joinpath(schema_dir, "vsmartmom-parameters.schema.json"))
    radiative_transfer = parameter_schema["definitions"]["radiative_transfer"]

    @test Set(parameter_schema["required"]) == Set([
        "radiative_transfer",
        "geometry",
        "atmospheric_profile",
    ])
    @test "Δ_angle" in radiative_transfer["required"]
    @test haskey(radiative_transfer["properties"], "Δ_angle")
end

@testset "Aqua" begin
    Aqua.test_all(vSmartMOM;
        # Keep method-ambiguity rollout separate; current branch still has
        # legacy broad signatures that make this too noisy for the first gate.
        ambiguities = false,
        # CUDA and lazy artifact precompilation can leave background tasks in
        # the test environment. Keep this check separate from the package
        # hygiene gate until those initialization paths are isolated.
        persistent_tasks = false,
    )
end
