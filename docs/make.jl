using Documenter
using Literate
using vSmartMOM
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using CairoMakie

function build()

    tutorials = ["Tutorial_QuickStart.jl", "Tutorial_Absorption.jl", "Tutorial_Scattering.jl", "Tutorial_MieDeepDive.jl", "Tutorial_IO.jl", "Tutorial_CoreRT.jl", "Tutorial_Surfaces.jl", "Tutorial_Canopy.jl", "Tutorial_Jacobians.jl", "Tutorial_GPU.jl", "Tutorial_HybridAD.jl"]
    tutorials_paths = [joinpath(@__DIR__, "src", "pages", "tutorials", tutorial) for tutorial in tutorials]

    for tutorial in tutorials_paths
        Literate.markdown(tutorial, joinpath(@__DIR__, "src", "pages", "tutorials");
                          codefence = "```julia" => "```")
    end

    tutorials_md = [joinpath("pages", "tutorials", tutorial[1:end-3]) * ".md" for tutorial in tutorials]

    pages = Any[
        "Getting Started"       => "index.md",
        "Quick Start (5 min)"   => "pages/quickstart.md",
        "Configure a Scene"     => "pages/IO/Overview.md",
        "Compute Jacobians"     => "pages/jacobians.md",
        "Run on GPU"            => "pages/gpu.md",
        "Extend vSmartMOM"      => Any[
                                    "Add a Surface BRDF" => "pages/extending/surfaces.md",
                                    "Add a Raman Mode" => "pages/extending/raman.md",
                                   ],
        "Architecture & Design" => "design.md",
        "Core RT Theory"        => "pages/vSmartMOM/CoreRTTheory.md",
        "API Reference"         => Any[
                                    "Public API" => "pages/api_reference.md",
                                    "Developer API Coverage" => "pages/internal_api_coverage.md",
                                   ],
        "Release Notes"         => "pages/release_notes.md",
        "IO Reference"          => Any[
                                    "Schema" => "pages/IO/Schema.md",
                                    "Examples" => "pages/IO/Examples.md",
                                   ],
        "vSmartMOM"             => Any[
                                    "Overview" => "pages/vSmartMOM/Overview.md",
                                    "Example" => "pages/vSmartMOM/Example.md",
                                    "User-Defined RT Parameters" => "pages/vSmartMOM/InputParametersGuide.md",
                                    "Methods & Types" => "pages/vSmartMOM/Types.md",
                                    "References" => "pages/vSmartMOM/References.md"
                                    ],
        "Absorption"            => [
                                    "Overview" => "pages/Absorption/Overview.md",
                                    "HITRAN Data Management" => "pages/Absorption/HITRAN_Data.md",
                                    "Example" => "pages/Absorption/Example.md",
                                    "Methods & Types" => "pages/Absorption/Types.md",
                                    "References" => "pages/Absorption/References.md"
                                    ],
        "Scattering"            => [
                                    "Overview" => "pages/Scattering/Overview.md",
                                    "Example" => "pages/Scattering/Example.md",
                                    "Methods & Types" => "pages/Scattering/Types.md",
                                    "References" => "pages/Scattering/References.md"
                                    ],
        "Surfaces"              => "pages/Surfaces/Overview.md",
        "Inelastic Scattering"  => "pages/Inelastic/Overview.md",
        "Aerosols"              => "pages/Aerosols/Overview.md",
        "SolarModel"            => "pages/SolarModel/Overview.md",
        "GEOSChem Integration"  => "pages/geoschem_integration.md",
        "Tutorials"             => tutorials_md
    ]

    mathengine = MathJax(Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ))

    # The format will make other pages in parallel with the index page
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = mathengine,
        collapselevel = 1,
        assets = ["assets/favicon.ico"],
        size_threshold_warn = 250 * 1024,
        size_threshold = 500 * 1024,
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

end

build()

# Deploy only in CI contexts; local docs builds should not attempt git deployment.
if get(ENV, "CI", "false") == "true"
    # Keep main as the canonical dev docs, but allow branch-specific dev deployment
    # for unified-vsmartmom so pushes to that branch publish docs as well.
    ref_name = get(ENV, "GITHUB_REF_NAME", "")
    deploy_devbranch = ref_name == "unified-vsmartmom" ? "unified-vsmartmom" : "main"
    deploy_devurl = ref_name == "unified-vsmartmom" ? "unified-vsmartmom" : "dev"

    deploydocs(
        repo = "github.com/RemoteSensingTools/vSmartMOM.jl.git",
        target = "build",
        push_preview = true,
        devbranch = deploy_devbranch,
        devurl = deploy_devurl,
    )
end
