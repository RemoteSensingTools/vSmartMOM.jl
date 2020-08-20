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

makedocs(sitename="Radiative Transfer", 
         pages=pages)
