include("./CrossSection.jl")
using .CrossSection.HITRAN
using .CrossSection.VoigtShape

df = readHITRAN("/Users/cfranken/Desktop/02_hit04.par", 2, 1, 0, 30000)

computeVoigt(df, collect(1:20000), 950, 298, 1, 100)
