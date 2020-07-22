include("./CrossSection.jl")
using .CrossSection.HITRAN
using .CrossSection.VoigtShape

df = readHITRAN("/net/fluo/data1/projects/RTM/Mie/hitran/par/02_hit04.par", 2, 1, 0, 30000)

computeVoigt(df, collect(1:20000), 950, 298, 1, 100)
