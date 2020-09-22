using Documenter, RadiativeTransfer, RadiativeTransfer.CrossSection

# generated_dir = joinpath(@__DIR__, "src", "generated") # generated files directory
# rm(generated_dir, force = true, recursive = true)
# mkpath(generated_dir)

pages = Any[
    "Home"           => "index.md",
    "CrossSection"   => "pages/CrossSection.md"
    ]

# mathengine = MathJax(Dict(
#     :TeX => Dict(
#         :equationNumbers => Dict(:autoNumber => "AMS"),
#         :Macros => Dict(),
#     ),
# ))

# format = Documenter.HTML(
#     prettyurls = get(ENV, "CI", nothing) == "true",
#     mathengine = mathengine,
#     collapselevel = 1,
# )

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
)

# This way it shows warnings of functions that have not been documented
makedocs(
    sitename = "Radiative Transfer",
    format = format,
    clean = false,
    modules = [RadiativeTransfer],
    pages = pages,
)

deploydocs(
    repo = "github.com/RupeshJey/RadiativeTransfer.jl.git",
    target = "build",
    push_preview = true,
)
