include("./CrossSection.jl")
using .CrossSection.HITRAN
using .CrossSection.DopplerShape
using .CrossSection.VoigtShape

df = readHITRAN("/net/fluo/data1/projects/RTM/Mie/hitran/par/02_hit04.par", 2, -1, 6300, 6380)

# computeVoigt(df, collect(1:20000), 950, 298, 1, 100)

computeVoigt(df, 6300, 6380, 0.001, 1, 290, 1, 100)
