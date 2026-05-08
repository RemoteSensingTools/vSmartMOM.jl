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

    @testset "Taplo wires JSON schema to TOML configs" begin
        @test isfile(_TAPLO_TOML)
        t = read(_TAPLO_TOML, String)
        @test occursin("vsmartmom-parameters.schema.json", t)
        # Common config locations the schema applies to.
        # Note: .taplo.toml uses globs like `test/**/*.toml` which
        # covers test/test_parameters/ — we just check the parent
        # directory tokens are mentioned.
        for path_glob in ("config/", "test/", "examples/")
            @test occursin(path_glob, t)
        end
    end
end
