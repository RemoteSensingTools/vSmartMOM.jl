using Documenter
using Literate
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using Plots 

ENV["GKSwstype"] = "nul"

function build()

    tutorials = ["Tutorial_Absorption.jl", "Tutorial_Scattering.jl"] # , 
    tutorials_paths = [joinpath(@__DIR__, "src", "pages", "tutorials", tutorial) for tutorial in tutorials]

    for tutorial in tutorials_paths
        Literate.markdown(tutorial, joinpath(@__DIR__, "src", "pages", "tutorials"))
    end

    tutorials_md = [joinpath("pages", "tutorials", tutorial[1:end-3]) * ".md" for tutorial in tutorials]

    pages = Any[
        "Getting Started"       => "index.md",
        "RadiativeTransfer"     => Any[
                                    "Overview" => "pages/RadiativeTransfer/Overview.md", 
                                    "Example" => "pages/RadiativeTransfer/Example.md",
                                    "User-Defined RT Parameters" => "pages/RadiativeTransfer/InputParametersGuide.md",
                                    "Methods & Types" => "pages/RadiativeTransfer/Types.md",
                                    "References" => "pages/RadiativeTransfer/References.md"
                                    ],
        "Absorption"            => [
                                    "Overview" => "pages/Absorption/Overview.md",
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
        assets = ["assets/favicon.ico"]
    )

    # This way it shows warnings of functions that have not been documented
    makedocs(
        sitename = "Radiative Transfer",
        format = format,
        clean = false,
        modules = [RadiativeTransfer],
        pages = pages,
    )

end

build()

deploydocs(
    repo = "github.com/RadiativeTransfer/RadiativeTransfer.jl.git",
    target = "build",
    push_preview = true,
)

