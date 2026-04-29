using Documenter
using DocumenterVitepress
using Literate
using vSmartMOM
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using CairoMakie

function copy_landing_page_icons()
    icon_src = joinpath(@__DIR__, "src", "assets", "icons")
    isdir(icon_src) || return nothing

    destinations = [joinpath(@__DIR__, "build", ".documenter", "public", "assets", "icons")]
    bases_file = joinpath(@__DIR__, "build", "bases.txt")
    if isfile(bases_file)
        bases = readlines(bases_file)
        append!(destinations,
                [joinpath(@__DIR__, "build", string(i), "assets", "icons")
                 for i in eachindex(bases)])
    end

    for dest in destinations
        mkpath(dest)
        for icon in readdir(icon_src)
            endswith(icon, ".svg") || continue
            cp(joinpath(icon_src, icon), joinpath(dest, icon); force = true)
        end
    end
    return nothing
end

function build()

    tutorials = ["Tutorial_QuickStart.jl", "Tutorial_Absorption.jl", "Tutorial_Scattering.jl", "Tutorial_MieDeepDive.jl", "Tutorial_IO.jl", "Tutorial_CoreRT.jl", "Tutorial_Surfaces.jl", "Tutorial_Canopy.jl", "Tutorial_Jacobians.jl", "Tutorial_GPU.jl", "Tutorial_HybridAD.jl"]
    tutorials_paths = [joinpath(@__DIR__, "src", "pages", "tutorials", tutorial) for tutorial in tutorials]

    for tutorial in tutorials_paths
        Literate.markdown(tutorial, joinpath(@__DIR__, "src", "pages", "tutorials");
                          codefence = "```julia" => "```")
    end

    tutorials_md = [joinpath("pages", "tutorials", tutorial[1:end-3]) * ".md" for tutorial in tutorials]

    pages = Any[
        "Home"                  => "index.md",
        "Manual"                => Any[
                                    "Quick Start (5 min)" => "pages/quickstart.md",
                                    "Configure a Scene" => "pages/IO/Overview.md",
                                    "Compute Jacobians" => "pages/jacobians.md",
                                    "Run on GPU" => "pages/gpu.md",
                                   ],
        "Concepts"              => Any[
                                    "Core RT Theory" => "pages/vSmartMOM/CoreRTTheory.md",
                                    "Architecture & Design" => "design.md",
                                    "Absorption" => "pages/Absorption/Overview.md",
                                    "Scattering" => "pages/Scattering/Overview.md",
                                    "Surfaces" => "pages/Surfaces/Overview.md",
                                    "Inelastic Scattering" => "pages/Inelastic/Overview.md",
                                    "Aerosols" => "pages/Aerosols/Overview.md",
                                    "Solar Model" => "pages/SolarModel/Overview.md",
                                   ],
        "Developer Guides"      => Any[
                                    "Add a Surface BRDF" => "pages/extending/surfaces.md",
                                    "Add a Raman Mode" => "pages/extending/raman.md",
                                    "GEOS-Chem Integration" => "pages/geoschem_integration.md",
                                   ],
        "Tutorials"             => tutorials_md,
        "Library"               => Any[
                                    "Overview" => "pages/api_reference.md",
                                    "Top-Level API" => "pages/api/top_level.md",
                                    "CoreRT" => "pages/api/core_rt.md",
                                    "IO" => "pages/api/io.md",
                                    "Absorption" => "pages/api/absorption.md",
                                    "Scattering" => "pages/api/scattering.md",
                                    "Surfaces" => "pages/api/surfaces.md",
                                    "Inelastic Scattering" => "pages/api/inelastic.md",
                                    "Aerosols" => "pages/api/aerosols.md",
                                    "SolarModel" => "pages/api/solar_model.md",
                                    "Experimental Helpers" => "pages/api/experimental.md",
                                    "Developer Coverage" => "pages/internal_api_coverage.md",
                                    "Function Index" => "pages/api/function_index.md",
                                   ],
        "Resources"             => Any[
                                    "Input Schema" => "pages/IO/Schema.md",
                                    "IO Examples" => "pages/IO/Examples.md",
                                    "HITRAN Data Management" => "pages/Absorption/HITRAN_Data.md",
                                    "References" => "pages/vSmartMOM/References.md",
                                    "Release Notes" => "pages/release_notes.md",
                                   ],
    ]

    ref_name = get(ENV, "GITHUB_REF_NAME", "")
    deploy_devbranch = ref_name == "unified-vsmartmom" ? "unified-vsmartmom" : "main"
    deploy_devurl = ref_name == "unified-vsmartmom" ? "unified-vsmartmom" : "dev"

    format = MarkdownVitepress(
        repo = "github.com/RemoteSensingTools/vSmartMOM.jl",
        devbranch = deploy_devbranch,
        devurl = deploy_devurl,
        deploy_decision = get(ENV, "CI", "false") == "true" ?
            nothing : Documenter.DeployDecision(all_ok = false),
        deploy_url = "https://remotesensingtools.github.io/vSmartMOM.jl",
        description = "Polarized atmospheric radiative transfer and remote sensing tools in Julia.",
        assets = [
            "assets/favicon.ico",
            "assets/logo.png",
            "assets/icons/logo.svg",
            "assets/icons/scattering.svg",
            "assets/icons/absorption.svg",
            "assets/icons/radiative_transfer.svg",
        ],
        sidebar_drawer = true,
    )

    # This way it shows warnings of functions that have not been documented
    makedocs(
        sitename = "vSmartMOM",
        format = format,
        clean = true,
        modules = [vSmartMOM],
        pages = pages,
        checkdocs = :exports,
    )

    copy_landing_page_icons()

end

build()

# Deploy only in CI contexts; local docs builds should not attempt git deployment.
if get(ENV, "CI", "false") == "true"
    # Keep main as the canonical dev docs, but allow branch-specific dev deployment
    # for unified-vsmartmom so pushes to that branch publish docs as well.
    ref_name = get(ENV, "GITHUB_REF_NAME", "")
    deploy_devbranch = ref_name == "unified-vsmartmom" ? "unified-vsmartmom" : "main"
    deploy_devurl = ref_name == "unified-vsmartmom" ? "unified-vsmartmom" : "dev"

    DocumenterVitepress.deploydocs(
        root = @__DIR__,
        repo = "github.com/RemoteSensingTools/vSmartMOM.jl.git",
        target = "build",
        push_preview = true,
        devbranch = deploy_devbranch,
        devurl = deploy_devurl,
    )
end
