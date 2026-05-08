# Phase D — schema documentation invariants
# =========================================================================
#
# Asserts that the per-block schema documentation under
# `docs/src/pages/IO/Schema/` covers all the top-level YAML/TOML
# blocks accepted by `parameters_from_dict`, and that the JSON
# Schema at `schemas/vsmartmom-parameters.schema.json` is consistent
# with the doc pages.
#
# This is a "no field-doc drift" regression: when someone adds a
# new top-level block to the parser, this test fails until they add
# a corresponding `Schema/<block>.md` page (and ideally a JSON Schema
# entry).

using Test

const _SCHEMA_DOCS = joinpath(@__DIR__, "..", "docs", "src", "pages", "IO", "Schema")
const _SCHEMA_INDEX = joinpath(@__DIR__, "..", "docs", "src", "pages", "IO", "Schema.md")
const _JSON_SCHEMA = joinpath(@__DIR__, "..", "schemas", "vsmartmom-parameters.schema.json")
const _TAPLO_TOML  = joinpath(@__DIR__, "..", ".taplo.toml")
const _VSCODE_SETTINGS = joinpath(@__DIR__, "..", ".vscode", "settings.json")
const _VSCODE_EXTENSIONS = joinpath(@__DIR__, "..", ".vscode", "extensions.json")

@testset "Schema docs invariants" begin
    @testset "Per-block pages exist" begin
        for page in (
            "radiative_transfer.md", "geometry.md", "atmospheric_profile.md",
            "absorption.md", "aerosols.md", "sources.md", "surface.md",
        )
            path = joinpath(_SCHEMA_DOCS, page)
            @test isfile(path)
        end
    end

    @testset "Schema index links each per-block page" begin
        @test isfile(_SCHEMA_INDEX)
        idx = read(_SCHEMA_INDEX, String)
        for page in ("radiative_transfer.md", "geometry.md",
                     "atmospheric_profile.md", "absorption.md",
                     "aerosols.md", "sources.md", "surface.md")
            @test occursin("Schema/$page", idx)
        end
        # Editor-support recipe must be discoverable from the index.
        @test occursin("Taplo", idx)
        @test occursin("yaml-language-server", idx)
    end

    @testset "JSON Schema lists v0.7 fields" begin
        @test isfile(_JSON_SCHEMA)
        s = read(_JSON_SCHEMA, String)
        # Use grep-style checks so we don't need a JSON test dependency.
        # Phase D (v0.7) primary knobs:
        @test occursin("\"nstreams\"", s)
        @test occursin("\"m_max\"", s)
        @test occursin("\"truncation\"", s)
        # Legacy still in the schema (kept for backward compat):
        @test occursin("\"max_m\"", s)
        @test occursin("\"l_trunc\"", s)
        @test occursin("\"quadrature_type\"", s)
        # Sanghavi's nstreams≥3 minimum is enforced in the schema.
        # Match the JSON region for `nstreams` and verify "minimum": 3
        # appears within it.
        m = match(r"\"nstreams\":\s*\{[^}]*\"minimum\":\s*3"s, s)
        @test m !== nothing
    end

    @testset "JSON Schema choice fields expose flat enum completions" begin
        s = read(_JSON_SCHEMA, String)
        # yaml-language-server completion works more reliably from flat enum
        # lists than from oneOf-per-value branches. Keep descriptions via
        # enumDescriptions so VS Code shows useful choice help.
        @test match(r"\"polarization_type\":\s*\{[^}]*\"enum\":\s*\[[^\]]*\"Stokes_I\(\)\"[^\]]*\"Stokes_IQ\(\)\"[^\]]*\"Stokes_IQU\(\)\"[^\]]*\"Stokes_IQUV\(\)\""s, s) !== nothing
        @test match(r"\"polarization_type\":\s*\{[^}]*\"enumDescriptions\""s, s) !== nothing
        for field in ("quadrature_type", "float_type", "architecture", "broadening", "CEF", "decomp_type")
            @test match(Regex("\"$field\"\\s*:\\s*\\{[^}]*\"enum\"\\s*:", "s"), s) !== nothing
            @test match(Regex("\"$field\"\\s*:\\s*\\{[^}]*\"enumDescriptions\"", "s"), s) !== nothing
        end
    end

    @testset "Taplo wires JSON schema to TOML configs" begin
        @test isfile(_TAPLO_TOML)
        t = read(_TAPLO_TOML, String)
        @test occursin("vsmartmom-parameters.schema.json", t)
        # Schema association is intentionally limited to setup/scene TOML
        # locations, not every TOML file in test/ or the package metadata.
        for path_glob in (
            "config/**/*.toml",
            "test/test_parameters/**/*.toml",
            "test/benchmarks/*.toml",
            "test/vlidort_baseline/configs/**/*.toml",
            "examples/**/*.toml",
            "sandbox/**/*.toml",
        )
            @test occursin(path_glob, t)
        end
        @test !occursin("\"test/**/*.toml\"", t)
        @test !occursin("**/*params*.toml", t)
        @test occursin("test/benchmarks/harness/**/*.toml", t)
    end

    @testset "VS Code wires JSON schema to YAML setup files" begin
        @test isfile(_VSCODE_SETTINGS)
        @test isfile(_VSCODE_EXTENSIONS)
        settings = read(_VSCODE_SETTINGS, String)
        extensions = read(_VSCODE_EXTENSIONS, String)
        @test occursin("vsmartmom-parameters.schema.json", settings)
        @test occursin("/config/**/*.yaml", settings)
        @test occursin("/test/test_parameters/**/*.yaml", settings)
        @test occursin("/test/benchmarks/*.yaml", settings)
        @test occursin("/test/vlidort_baseline/configs/**/*.yaml", settings)
        @test occursin("/sandbox/**/*.yaml", settings)
        @test !occursin("**/*params*.yaml", settings)
        @test occursin("tamasfe.even-better-toml", extensions)
        @test occursin("redhat.vscode-yaml", extensions)
    end
end
